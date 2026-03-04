#include "atmosphere_step_solver.h"
#include "fluids/OceanAlbedoModel.h"

#include "atmospheric_thermodynamics.h"
#include "surface_atmosphere_coupler.h"
#include "atmosphere/EvaporationModel.h"
#include "fluids/PhaseModel.h"

#include <QtConcurrent/QtConcurrentMap>
#include <QtCore/QLoggingCategory>
#include <QtCore/QString>
#include <QtCore/QtMath>

#include <cmath>
#include <numeric>

Q_LOGGING_CATEGORY(atmosphereProfileLog, "solar.atmosphere.profile")
Q_LOGGING_CATEGORY(atmosphereMassLog, "solar.atmosphere.mass")

namespace {
// Базовая скорость массообмена для испарения (м/с).
// Bulk transfer: C_E ≈ 1.5e-3, типичный ветер ~5 м/с → C_E * U ≈ 7.5e-3.
// Минимальное значение при штиле (~0.5 м/с ветра).
constexpr double kEvaporationBaseTransferVelocityMps = 7.5e-4;
// Коэффициент Дальтона (безразмерный): типичный C_E ~ 1.2e-3 для океана.
constexpr double kDaltonCoefficient = 1.2e-3;
constexpr double kEvaporationTempMinK = 260.0;
constexpr double kEvaporationTempRangeK = 20.0;
constexpr double kWaterMassTolerance = 1.0e-6;
// Удельная теплота парообразования воды, Дж/кг.
constexpr double kLatentHeatVaporization = 2.501e6;
// Поглощение H₂O в коротковолновом ближнем ИК (слабее, чем длинноволновое).
constexpr double kH2OMassAbsorptionSw = 0.01;
// Доля сублимации: ледяной океан обеспечивает ~10% влаги от жидкого.
// В реальности сублимация льда при -10°C ≈ 10–15% от испарения воды при 0°C.
constexpr double kIceSublimationFraction = 0.1;

// ── Многополосная спектральная модель парникового эффекта ──────────────
//
// 4 ИК-полосы с долями Планка при ~280K:
//   Band 0  Far-IR   (>17 мкм)  — вращательная полоса H₂O, CIA CO₂ при высоком P
//   Band 1  Window   (8–13 мкм) — «окно прозрачности»: SF₆, NF₃, NH₃
//   Band 2  CO₂      (13–17 мкм) — основная полоса CO₂
//   Band 3  H₂O+CH₄  (5–8 мкм) — колебательная полоса H₂O, метан
//
// Каждый газ имеет собственные коэффициенты поглощения в каждой полосе.
// Солвер выполняет независимый двухпоточный перенос для каждой полосы,
// корректно моделируя независимое насыщение:
// SF₆ эффективен даже когда CO₂+H₂O полностью насыщены (поглощает в окне).
constexpr int kNumBands = kRadiationBandCount; // из layered_radiation_solver.h

// H₂O: линейное массовое поглощение по полосам (м²/кг).
// Взвешенное среднее ≈ 0.1 м²/кг (совместимо со старой однополосной моделью).
constexpr double kH2OBandLw[kNumBands] = {0.17, 0.01, 0.08, 0.22};

// Логарифмическая формула для газов: τ = coeff × log1p(scale × column × pressure).
// Массив logScale по полосам; коэффициент logCoeff одинаковый (0.16).
constexpr double kLogCoeff = 0.16;

// CO₂: основное поглощение в полосе 15 мкм (band 2), слабое в остальных.
// Band 0 (Far-IR): при высоком давлении (>10 атм) collision-induced absorption
// и давленческое уширение крыльев 15-мкм полосы дают заметное поглощение.
// Значение 0.3 обеспечивает ~τ=2 на слой для Венеры (92 атм, 96% CO₂).
constexpr double kCO2BandScale[kNumBands] = {0.3, 0.5, 10.0, 0.5};
// SF₆: поглощение 10.5 мкм — почти целиком в окне (band 1). GWP ≈ 23500.
constexpr double kSF6BandScale[kNumBands] = {0.0, 200000.0, 0.0, 0.0};
// NF₃: поглощение ~10 мкм, аналогично SF₆. GWP ≈ 17200.
constexpr double kNF3BandScale[kNumBands] = {0.0, 150000.0, 0.0, 0.0};
// CH₄: 7.7 мкм и 3.3 мкм — в основном band 3, слабо в окне. GWP ≈ 30.
constexpr double kCH4BandScale[kNumBands] = {0.0, 0.1, 0.0, 0.5};
// NH₃: 10–11 мкм (окно) + широкие вращательные полосы. GWP не определён (короткоживущий).
constexpr double kNH3BandScale[kNumBands] = {0.1, 500.0, 0.1, 0.1};
// H₂: collision-induced absorption (CIA). Значим только для массивных H₂-атмосфер.
constexpr double kH2BandScale[kNumBands] = {0.001, 0.001, 0.001, 0.001};

// Вычисляет оптические толщины по полосам для одного слоя.
// baseTau — базовая τ из инициализации (серая, добавляется к каждой полосе).
// h2oKg — водяной пар (кг/м²); gasColumns[6] — CO₂, SF₆, NF₃, CH₄, NH₃, H₂ (кг/м²);
// pressure — давление (атм); outTau[kNumBands] — результат.
inline void computeBandTaus(double baseTau,
                            double h2oKg,
                            const double gasColumns[6],
                            double pressure,
                            double outTau[kNumBands]) {
    const double *bandScales[6] = {
        kCO2BandScale, kSF6BandScale, kNF3BandScale,
        kCH4BandScale, kNH3BandScale, kH2BandScale};
    for (int b = 0; b < kNumBands; ++b) {
        double tau = baseTau + kH2OBandLw[b] * h2oKg;
        for (int g = 0; g < 6; ++g) {
            const double scale = bandScales[g][b];
            if (scale > 0.0 && gasColumns[g] > 0.0) {
                tau += kLogCoeff *
                    std::log1p(scale * gasColumns[g] * pressure);
            }
        }
        outTau[b] = tau;
    }
}

// ── Время жизни газов (распад/фотолиз, в секундах) ────────────────────
// CH₄: ~12 лет (реакция с OH-радикалом, УФ-фотолиз).
constexpr double kCH4LifetimeSeconds = 12.0 * 365.25 * 86400.0;
// NH₃: ~14 дней (быстрое осаждение, реакция с кислотами).
constexpr double kNH3LifetimeSeconds = 14.0 * 86400.0;
// SF₆: ~3200 лет (крайне стабильный).
constexpr double kSF6LifetimeSeconds = 3200.0 * 365.25 * 86400.0;
// NF₃: ~550 лет.
constexpr double kNF3LifetimeSeconds = 550.0 * 365.25 * 86400.0;
// H₂: ~2 года (реакция с OH, утечка в космос).
constexpr double kH2LifetimeSeconds = 2.0 * 365.25 * 86400.0;

// ── Атмосферное убегание (Jeans escape) ───────────────────────────────
constexpr double kBoltzmannConstant = 1.380649e-23;   // Дж/К
constexpr double kGravConst = 6.67430e-11;             // м³/(кг·с²)
constexpr double kEarthMass = 5.9722e24;               // кг
constexpr double kAmuToKg = 1.66054e-27;               // кг/а.е.м.
// Молярные массы (г/моль ≈ а.е.м.).
constexpr double kMolarMassH2 = 2.016;
constexpr double kMolarMassHe = 4.003;
constexpr double kMolarMassCH4 = 16.04;
constexpr double kMolarMassNH3 = 17.03;
// Время сглаживания для интенсивности осадков (EMA) в секундах.
constexpr double kPrecipitationEmaTimeSeconds = 3600.0;

double potentialEvaporationKgPerM2(const SurfacePointState &surface,
                                   const AtmosphericLayerState &layer,
                                   double windSpeedMps,
                                   double dtSeconds) {
    if (dtSeconds <= 0.0) {
        return 0.0;
    }

    const double surfaceTemperature = surface.temperatureKelvin();
    if (surfaceTemperature <= 0.0) {
        return 0.0;
    }

    // Температурный множитель: испарение растёт при переходе через диапазон 260–280 K.
    const double temperatureFactor = qBound(0.0,
                                            (surfaceTemperature - kEvaporationTempMinK) /
                                                kEvaporationTempRangeK,
                                            1.0);
    if (temperatureFactor <= 0.0) {
        return 0.0;
    }

    const double saturationDensity =
        EvaporationModel::saturationVaporDensityKgPerM3(surfaceTemperature);
    const double layerThickness = layer.thicknessMeters();
    const double vaporDensity = (layerThickness > 0.0)
        ? qMax(0.0, layer.waterVaporKgPerM2()) / layerThickness
        : 0.0;
    const double deficit = qMax(0.0, saturationDensity - vaporDensity);
    // Bulk aerodynamic formula: E = C_E * U * (rho_sat - rho),
    // где C_E — коэффициент Дальтона (~1.2e-3), U — скорость ветра.
    // При штиле используем минимальную базовую скорость переноса.
    const double effectiveWind = qMax(0.5, std::abs(windSpeedMps));
    const double transferVelocity =
        qMax(kEvaporationBaseTransferVelocityMps, kDaltonCoefficient * effectiveWind);
    const double evaporationRateKgPerM2S = transferVelocity * deficit * temperatureFactor;
    return evaporationRateKgPerM2S * dtSeconds;
}

double columnWaterMassKgPerM2(const AtmosphericColumn &column) {
    double total = 0.0;
    const auto &layers = column.layers();
    for (const auto &layer : layers) {
        total += qMax(0.0, layer.waterVaporKgPerM2());
        total += qMax(0.0, layer.liquidWaterKgPerM2());
        total += qMax(0.0, layer.iceWaterKgPerM2());
    }
    return total;
}

// Вертикальная турбулентная диффузия температуры между слоями.
// Дополняет конвективную коррекцию: конвекция исправляет сверхадиабатические
// градиенты (нижний слой горячее), а диффузия рассасывает устойчивые инверсии
// (верхний слой горячее). Без этого горячие верхние слои остывают только
// радиационно, что при малой оптической толщине занимает сотни шагов.
//
// ВАЖНО: используем плотность × Cp × h (теплоёмкость слоя) для определения
// скорости изменения температуры. Без учёта массы плотный нижний слой (1 атм)
// терял бы тепло к разрежённому верхнему (0.01 атм) в ~100 раз быстрее
// реального → искусственное переохлаждение поверхности.
void applyVerticalTemperatureDiffusion(QVector<AtmosphericLayerState> &layers,
                                       double kz, double dtSeconds) {
    const int n = layers.size();
    if (n < 2 || kz <= 0.0 || dtSeconds <= 0.0) {
        return;
    }
    constexpr double kMinDz = 1.0e-3;
    constexpr double kMinHeatCapacity = 1.0;
    constexpr double kMaxCourant = 0.5;

    QVector<double> temps(n);
    QVector<double> thickness(n);
    QVector<double> heatCap(n);    // C = ρ × Cp × h [J/(m²·K)]
    QVector<double> density(n);    // ρ [kg/m³]
    double minDz = 1.0e30;
    for (int i = 0; i < n; ++i) {
        temps[i] = layers[i].temperatureKelvin();
        thickness[i] = qMax(layers[i].thicknessMeters(), kMinDz);
        heatCap[i] = qMax(layers[i].heatCapacityJPerM2K(), kMinHeatCapacity);
        density[i] = qMax(layers[i].densityKgPerM3(), 1.0e-6);
        minDz = qMin(minDz, thickness[i]);
    }
    minDz = qMax(minDz, kMinDz);

    // Ограничение диффузионного числа Куранта: dt <= C_vol * dz² / (ρ × Cp × Kz).
    // Для нижнего слоя (максимальная плотность) Courant наибольший.
    const double maxStableDt = kMaxCourant * minDz * minDz / kz;
    if (maxStableDt <= 0.0) {
        return;
    }

    const int substeps = qMax(1, static_cast<int>(qCeil(dtSeconds / maxStableDt)));
    const double stepDt = dtSeconds / static_cast<double>(substeps);
    const int ifaceCount = n - 1;
    QVector<double> flux(ifaceCount);      // Энергетический поток [W/m²]
    QVector<double> nextTemps(n);

    for (int step = 0; step < substeps; ++step) {
        for (int i = 0; i < ifaceCount; ++i) {
            const double dz = qMax(0.5 * (thickness[i] + thickness[i + 1]), kMinDz);
            // Плотность на границе: среднее из соседних слоёв.
            const double rhoInterface = 0.5 * (density[i] + density[i + 1]);
            // Энергетический поток: F = -ρ × Cp × Kz × dT/dz [W/m²].
            // Cp ≈ 1004 J/(kg·K) для сухого воздуха.
            // Реальный поток масштабируется плотностью на границе, чтобы
            // обмен между плотным нижним и разрежённым верхним слоями
            // соответствовал реальному турбулентному переносу тепла.
            constexpr double kAirCp = 1004.0;
            flux[i] = -rhoInterface * kAirCp * kz *
                (temps[i + 1] - temps[i]) / dz;
        }

        for (int k = 0; k < n; ++k) {
            // Нижняя граница: нулевой поток (поверхность обменивается через coupler).
            const double fluxBelow = (k > 0) ? flux[k - 1] : 0.0;
            // Верхняя граница: нулевой поток (нет перемешивания через верх атмосферы).
            const double fluxAbove = (k < ifaceCount) ? flux[k] : 0.0;
            // ΔT = ΔF × dt / C_layer, где C_layer = ρ × Cp × h [J/(m²·K)].
            nextTemps[k] = qMax(0.0,
                temps[k] + (fluxBelow - fluxAbove) * stepDt / heatCap[k]);
        }

        temps.swap(nextTemps);
    }

    for (int i = 0; i < n; ++i) {
        layers[i].setTemperatureKelvin(temps[i]);
    }
}

double gasShareById(const AtmosphereComposition &composition, const QString &gasId) {
    const auto fractions = composition.fractions();
    for (const auto &fraction : fractions) {
        if (fraction.id == gasId) {
            return qMax(0.0, fraction.share);
        }
    }
    return 0.0;
}
} // namespace

AtmosphereStepSolver::AtmosphereStepSolver(const AtmosphereComposition &composition,
                                           double gravityMps2,
                                           double timeStepSeconds,
                                           double dayLengthSeconds,
                                           bool isRetrograde,
                                           double planetMassEarths,
                                           double planetRadiusKm)
    : radiationSolver_(timeStepSeconds)
    , convectiveSolver_(composition, gravityMps2)
    , gravityMps2_(gravityMps2)
    , rSpecific_(AtmosphericThermodynamics::specificGasConstant(composition))
    , specificHeatCp_(AtmosphericThermodynamics::specificHeatCp(composition))
    , timeStepSeconds_(timeStepSeconds)
    , dayLengthSeconds_(dayLengthSeconds)
    , isRetrograde_(isRetrograde)
    , co2Share_(gasShareById(composition, QStringLiteral("co2")))
    , sf6Share_(gasShareById(composition, QStringLiteral("sf6")))
    , nf3Share_(gasShareById(composition, QStringLiteral("nf3")))
    , ch4Share_(gasShareById(composition, QStringLiteral("ch4")))
    , nh3Share_(gasShareById(composition, QStringLiteral("nh3")))
    , h2Share_(gasShareById(composition, QStringLiteral("h2")))
{
    Q_UNUSED(planetMassEarths);
    Q_UNUSED(planetRadiusKm);
}

const SurfaceMaterial &AtmosphereStepSolver::materialForPoint(
    const QHash<QString, SurfaceMaterial> &materialsById,
    const SurfaceMaterial &defaultMaterial,
    const QString &materialId) const {
    const auto it = materialsById.constFind(materialId);
    return it != materialsById.cend() ? *it : defaultMaterial;
}

void AtmosphereStepSolver::updateLayerPressures(double surfacePressureAtm,
                                                QVector<AtmosphericLayerState> &layers) const {
    if (surfacePressureAtm <= 0.0 || rSpecific_ <= 0.0 || gravityMps2_ <= 0.0) {
        return;
    }

    double pressureAtmAtBottom = surfacePressureAtm;
    for (int layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
        const double layerTemperatureKelvin = layers.at(layerIndex).temperatureKelvin();
        const double layerThicknessMeters = layers.at(layerIndex).thicknessMeters();
        double layerPressureAtm = 0.0;
        if (pressureAtmAtBottom > 0.0 && layerTemperatureKelvin > 0.0 &&
            layerThicknessMeters > 0.0) {
            // Масштабная высота: H = R_specific * T / g, где
            // R_specific [Дж/(кг·К)], T [К], g [м/с²], H [м].
            const double scaleHeightMeters =
                (rSpecific_ * layerTemperatureKelvin) / gravityMps2_;
            if (scaleHeightMeters > 0.0) {
                // Гидростатика: dP/dz = -P/H ⇒ P(z+dz) = P(z) * exp(-dz/H).
                // Давление в середине слоя: P_mid = P_bottom * exp(-0.5 * dz / H) [атм].
                layerPressureAtm =
                    pressureAtmAtBottom * qExp(-0.5 * layerThicknessMeters / scaleHeightMeters);
                pressureAtmAtBottom =
                    pressureAtmAtBottom * qExp(-layerThicknessMeters / scaleHeightMeters);
            }
        }
        layers[layerIndex].setPressureAtm(layerPressureAtm);

        // Согласуем P–T–rho–C (давление–температура–плотность–теплоёмкость),
        // чтобы термодинамика оставалась устойчивой при изменении профиля.
        double densityKgPerM3 = 0.0;
        if (layerPressureAtm > 0.0 && layerTemperatureKelvin > 0.0 && rSpecific_ > 0.0) {
            const double pressurePa = layerPressureAtm * 101325.0;
            // Уравнение состояния идеального газа: rho = P / (R_specific * T).
            densityKgPerM3 = pressurePa / (rSpecific_ * layerTemperatureKelvin);
        }
        layers[layerIndex].setDensityKgPerM3(densityKgPerM3);

        const double heatCapacityJPerM2K =
            (densityKgPerM3 > 0.0 && specificHeatCp_ > 0.0 && layerThicknessMeters > 0.0)
                ? densityKgPerM3 * specificHeatCp_ * layerThicknessMeters
                : 0.0;
        layers[layerIndex].setHeatCapacityJPerM2K(heatCapacityJPerM2K);
    }
}

void AtmosphereStepSolver::runLayeredStep(const LayeredStepInput &input) {
    const int pointCount = input.surfaceGrid.points().size();
    const int columnCount = input.atmosphereGrid.columnCount();
    const int processedCount = qMin(pointCount, columnCount);
    if (processedCount <= 0) {
        return;
    }
    verticalWindMixingSolver_.setMixingCoefficient(input.verticalWindMixingCoefficientKz);
    verticalMoistureMixingSolver_.setMixingCoefficient(input.verticalMoistureMixingCoefficientKz);
    int logPointIndex = input.logPointIndex;
    if (logPointIndex < 0 || logPointIndex >= processedCount) {
        // Если индекс не задан, пишем лог для первой ячейки: так проще сравнивать шаги.
        logPointIndex = 0;
    }

    const double minTopPressureAtm = input.atmosphereGrid.minTopPressureAtm();

    // Структура для отложенного изменения числа слоёв в колонке.
    // Собираем запросы параллельно, применяем последовательно,
    // чтобы не вызывать updateColumnLayerCountFixedThickness из нескольких потоков.
    struct LayerResizeRequest {
        int newLayerCount = 0;
        double layerThicknessMeters = 0.0;
    };
    QVector<LayerResizeRequest> resizeRequests(processedCount);

    // Индексы точек для параллельной обработки.
    QVector<int> indices(processedCount);
    std::iota(indices.begin(), indices.end(), 0);

    // ── Фаза 1: поколоночная физика (параллельно) ──────────────────────
    // Каждая точка/колонка обрабатывается независимо: испарение, перемешивание
    // влаги, осадки, радиация, конвекция, теплообмен с поверхностью.
    // Все солверы используют const-методы, данные разделены по индексу i.
    //
    // ВАЖНО: QVector использует implicit sharing (copy-on-write).
    // Неконстантный operator[] вызывает detach(), который НЕ потокобезопасен.
    // Принудительно отделяем данные ДО параллельной секции через data().
    SurfacePoint *pointsPtr = input.surfaceGrid.points().data();
    AtmosphericColumn *columnsPtr = input.atmosphereGrid.columns().data();

    QtConcurrent::blockingMap(indices, [&, pointsPtr, columnsPtr](int i) {
        auto &point = pointsPtr[i];
        AtmosphericColumn &column = columnsPtr[i];
        const double localInsolation =
            (i < input.localInsolations.size()) ? input.localInsolations.at(i) : 0.0;
        const SurfaceMaterial &material =
            materialForPoint(input.materialsById, input.defaultMaterial, point.materialId);
        double surfaceAlbedo = qBound(0.0, material.albedo, 1.0);
        if (point.materialId == QStringLiteral("ocean")) {
            const double cosZenith =
                (i < input.localCosZeniths.size()) ? input.localCosZeniths.at(i) : 0.0;
            surfaceAlbedo = OceanAlbedoModel::albedoForPhase(point.waterPhase, cosZenith);
        }

        auto &layers = column.layers();
        if (!layers.isEmpty()) {
            // Океан — неограниченный источник воды для испарения.
            // Жидкий океан даёт полный запас влаги, замёрзший — сублимацию
            // (≈10% от испарения). Без сублимации атмосфера полностью
            // теряет водяной пар → исчезает парниковый эффект H₂O →
            // необратимый снежок (ice-albedo death spiral).
            if (point.materialId == QLatin1String("ocean")) {
                const double maxWater =
                    point.surfaceMoisture.settings().maxStorageKgPerM2;
                if (point.waterPhase == PhaseModel::Phase::Liquid) {
                    point.surfaceMoisture.setWaterKgPerM2(maxWater);
                } else {
                    // Сублимация: лёд → пар напрямую, без промежуточной
                    // жидкой фазы. Скорость ~10–15% от испарения при 0°C.
                    point.surfaceMoisture.setWaterKgPerM2(
                        maxWater * kIceSublimationFraction);
                }
            }

            // Испарение из поверхностного слоя: переносим влагу в нижний слой атмосферы.
            const double windSpeed =
                std::hypot(layers.first().windUMps(), layers.first().windVMps());
            const double potentialEvaporation =
                potentialEvaporationKgPerM2(point.state, layers.first(),
                                            windSpeed, timeStepSeconds_);
            const double evaporationKgPerM2 =
                point.surfaceMoisture.applyEvaporation(potentialEvaporation);
            if (evaporationKgPerM2 > 0.0) {
                layers[0].setWaterVaporKgPerM2(
                    layers[0].waterVaporKgPerM2() + evaporationKgPerM2);
                // Испарительное охлаждение поверхности: Q_evap = L_v * E.
                // Эта энергия уходит из поверхности в скрытую теплоту водяного пара.
                // На Земле ~80 Вт/м² (глобальное среднее).
                const double evaporativeCoolingWPerM2 =
                    kLatentHeatVaporization * evaporationKgPerM2 / qMax(1.0, timeStepSeconds_);
                point.state.applySurfaceFlux(-evaporativeCoolingWPerM2, timeStepSeconds_);
            }
        }

        // Фазовый баланс влаги обновляем до радиации, чтобы конденсация влияла
        // на альбедо облаков уже в текущем шаге.
        // Вертикальное перемешивание влаги выполняем перед фазовым балансом,
        // чтобы конденсация пересчитывалась уже после переноса влаги вверх/вниз.
        verticalMoistureMixingSolver_.mix(column, timeStepSeconds_, minTopPressureAtm);

        const double columnWaterBefore = columnWaterMassKgPerM2(column);
        const double precipitationKgPerM2 =
            evaporationModel_.updateColumnWithPrecipitation(column, timeStepSeconds_);
        const double columnWaterAfter =
            columnWaterMassKgPerM2(column) + precipitationKgPerM2;
        const double massDelta = columnWaterAfter - columnWaterBefore;
        const double massScale = qMax(1.0, columnWaterBefore);
        if (std::abs(massDelta) > kWaterMassTolerance * massScale) {
            qCWarning(atmosphereMassLog)
                << "Mass balance warning"
                << "index=" << i
                << "before=" << columnWaterBefore
                << "after=" << columnWaterAfter
                << "delta=" << massDelta;
        }
        point.precipitationRateKgPerM2 = precipitationKgPerM2;
        const double precipitationAlpha =
            (timeStepSeconds_ > 0.0 && kPrecipitationEmaTimeSeconds > 0.0)
                ? (1.0 - std::exp(-timeStepSeconds_ / kPrecipitationEmaTimeSeconds))
                : 1.0;
        // Храним EMA интенсивности осадков (кг/м² за интервал), а не интеграл по времени.
        point.precipitationKgPerM2 =
            (1.0 - precipitationAlpha) * point.precipitationKgPerM2 +
            precipitationAlpha * precipitationKgPerM2;
        if (precipitationKgPerM2 > 0.0) {
            const bool isOcean = (point.materialId == QLatin1String("ocean"));
            const bool isLiquidOcean = isOcean &&
                (point.waterPhase == PhaseModel::Phase::Liquid);
            const double phaseTemperatureK =
                (point.temperatureK > 0.0) ? point.temperatureK : point.airTemperatureK;
            // Упрощённый фазовый порог: ниже 273.15 K осадки считаем снегом.
            // Снег накапливается на суше и замёрзшем океане (морской лёд).
            const bool isSnow =
                (phaseTemperatureK > 0.0 && phaseTemperatureK <= 273.15) &&
                !isLiquidOcean;
            if (isSnow) {
                point.snowKgPerM2 += precipitationKgPerM2;
            } else if (!isLiquidOcean) {
                // На жидком океане дождь просто возвращается в воду —
                // surfaceMoisture не накапливается.
                point.surfaceMoisture.addPrecipitation(precipitationKgPerM2);
            }
        }

        const double condensationAlbedo =
            evaporationModel_.cloudAlbedoFromCondensation(column);
        // Используем индикатор конденсации как визуальную непрозрачность облаков.
        point.cloudOpacity = qBound(0.0, condensationAlbedo, 1.0);
        const double cloudShortwaveTransmission =
            qBound(0.0,
                   input.cloudShortwaveTransmission * (1.0 - condensationAlbedo),
                   1.0);

        // Многополосный парниковый эффект: вычисляем τ по полосам для каждого слоя,
        // затем солвер делает независимый двухпоточный перенос для каждой полосы.
        // Это корректно обрабатывает насыщение: SF₆ в окне работает независимо от CO₂/H₂O.
        QVector<LayeredRadiationSolver::BandOpticalDepths> bandTaus(layers.size());
        for (int li = 0; li < layers.size(); ++li) {
            const double h2oKg = qMax(0.0, layers.at(li).waterVaporKgPerM2());
            const double layerDensity = qMax(0.0, layers.at(li).densityKgPerM3());
            const double layerThickness = qMax(0.0, layers.at(li).thicknessMeters());
            const double pressure = qMax(0.0, layers.at(li).pressureAtm());
            const double gasColumns[6] = {
                co2Share_ * layerDensity * layerThickness,
                sf6Share_ * layerDensity * layerThickness,
                nf3Share_ * layerDensity * layerThickness,
                ch4Share_ * layerDensity * layerThickness,
                nh3Share_ * layerDensity * layerThickness,
                h2Share_  * layerDensity * layerThickness,
            };
            // Базовая τ (из инициализации) распределяется как серая (во все полосы).
            const double baseTau = qMax(0.0, layers.at(li).opticalDepthLongwave());
            computeBandTaus(baseTau, h2oKg, gasColumns, pressure,
                            bandTaus[li].tauLw);
            // Коротковолновая поправка H₂O (одна полоса).
            layers[li].setOpticalDepthShortwave(
                layers.at(li).opticalDepthShortwave() + kH2OMassAbsorptionSw * h2oKg);
        }

        // Радиационный расчёт: независимый двухпоточный перенос по ИК-полосам.
        LayeredRadiationSolver::SurfaceRadiativeFluxes surfaceRadFluxes;
        const QVector<double> layerDeltas =
            radiationSolver_.solveMultiband(column,
                                            bandTaus,
                                            localInsolation,
                                            surfaceAlbedo,
                                            cloudShortwaveTransmission,
                                            point.state.temperatureKelvin(),
                                            surfaceRadFluxes);
        const int layerCount = qMin(layers.size(), layerDeltas.size());
        for (int layerIndex = 0; layerIndex < layerCount; ++layerIndex) {
            const double updatedTemperature =
                qMax(0.0, layers.at(layerIndex).temperatureKelvin() + layerDeltas.at(layerIndex));
            layers[layerIndex].setTemperatureKelvin(updatedTemperature);
        }

        // Облачное поглощение: энергия, поглощённая плотными облаками (Венера и т.п.),
        // депозитируется в верхнюю треть атмосферных слоёв.
        // Без этого 20% инсоляции Венеры (~130 Вт/м²) просто исчезает.
        if (input.cloudAbsorptionFraction > 0.0 && localInsolation > 0.0 &&
            layerCount >= 2) {
            const double cloudAbsorbedWPerM2 =
                localInsolation * input.cloudAbsorptionFraction;
            // Депозитируем в верхнюю треть слоёв (≥1 слой) с весами ∝ exp.
            const int depositLayers = qMax(1, layerCount / 3);
            const int startLayer = layerCount - depositLayers;
            // Экспоненциальные веса: верхний слой получает больше всего.
            double weightSum = 0.0;
            QVector<double> weights(depositLayers);
            for (int k = 0; k < depositLayers; ++k) {
                weights[k] = std::exp(static_cast<double>(k));
                weightSum += weights[k];
            }
            if (weightSum > 0.0) {
                for (int k = 0; k < depositLayers; ++k) {
                    const int li = startLayer + k;
                    const double hc = layers.at(li).heatCapacityJPerM2K();
                    if (hc > 0.0) {
                        const double fraction = weights[k] / weightSum;
                        const double dT = cloudAbsorbedWPerM2 * fraction *
                                          timeStepSeconds_ / hc;
                        layers[li].setTemperatureKelvin(
                            qMax(0.0, layers.at(li).temperatureKelvin() +
                                      qBound(-20.0, dT, 20.0)));
                    }
                }
            }
        }

        // Восстанавливаем коротковолновую τ (базовая LW не менялась).
        for (int li = 0; li < layers.size(); ++li) {
            const double h2oKg = qMax(0.0, layers.at(li).waterVaporKgPerM2());
            layers[li].setOpticalDepthShortwave(
                layers.at(li).opticalDepthShortwave() - kH2OMassAbsorptionSw * h2oKg);
        }

        // Применяем радиационный баланс к поверхности:
        // Net = SW_absorbed + LW_down(backradiation) - ε σ T_surface^4.
        // Раньше поверхность не получала ни SW от атмосферы, ни backradiation.
        {
            const double surfaceNetRadiativeFlux =
                surfaceRadFluxes.shortwaveAbsorbedWPerM2 +
                surfaceRadFluxes.longwaveDownWPerM2;
            // Сохраняем стабилизированный поток в грунт ДО updateTemperature,
            // чтобы зафиксировать значение при текущей T_surface (аналогично
            // нелойэрному пути в solar_display.cpp, где используется прямое
            // присвоение, а не +=).
            point.subsurfaceFluxWPerM2 =
                point.state.stabilizedRadiativeFlux(surfaceNetRadiativeFlux, timeStepSeconds_);
            point.state.updateTemperature(surfaceNetRadiativeFlux, 0.0, timeStepSeconds_);
            point.temperatureK = point.state.temperatureKelvin();
        }
        // Сохраняем потоки для диагностики.
        point.shortwaveSurfaceWPerM2 = surfaceRadFluxes.shortwaveAbsorbedWPerM2;
        point.longwaveDownWPerM2 = surfaceRadFluxes.longwaveDownWPerM2;
        point.longwaveUpWPerM2 = point.state.surfaceEmittedFlux();

        convectiveSolver_.adjust(column);

        // Вертикальная диффузия температуры — турбулентное перемешивание сверх
        // конвективной коррекции: рассасывает устойчивые инверсии (горячий слой
        // над холодным), которые конвекция не трогает.
        applyVerticalTemperatureDiffusion(
            layers, input.verticalMoistureMixingCoefficientKz, timeStepSeconds_);

        const double initialAirTemperature =
            (point.airTemperatureK > 0.0) ? point.airTemperatureK : point.temperatureK;
        if (layers.isEmpty()) {
            point.airTemperatureK = initialAirTemperature;
            return;
        }

        point.airTemperatureK = layers.first().temperatureKelvin();

        // Снежная изоляция: термическое сопротивление R = d / k_snow.
        // Снег (плотность ~300 кг/м³, теплопроводность ~0.2 Вт/(м·К))
        // замедляет теплообмен поверхности с воздухом.
        double snowResistance = 0.0;
        if (point.snowKgPerM2 > 0.0) {
            constexpr double kSnowConductivity = 0.2;    // Вт/(м·К)
            constexpr double kSnowDensityKgPerM3 = 300.0; // кг/м³
            const double snowDepthM = point.snowKgPerM2 / kSnowDensityKgPerM3;
            snowResistance = snowDepthM / kSnowConductivity;
        }

        SurfaceAtmosphereCoupler coupler(input.heatTransferCoefficientWPerM2K);
        double surfaceAirFluxWPerM2 = 0.0;
        coupler.exchangeHeat(point.state,
                             layers[0],
                             material.roughnessLengthMeters,
                             timeStepSeconds_,
                             &surfaceAirFluxWPerM2,
                             snowResistance);
        point.airTemperatureK = layers[0].temperatureKelvin();
        point.temperatureK = point.state.temperatureKelvin();
        point.surfaceAirFluxWPerM2 = surfaceAirFluxWPerM2;
        // Поток в грунт учитывает радиацию и обмен с воздухом.
        point.subsurfaceFluxWPerM2 += -surfaceAirFluxWPerM2;

        updateLayerPressures(point.pressureAtm, layers);

        // Собираем запросы на изменение числа слоёв (применим после параллельной фазы).
        if (minTopPressureAtm > 0.0 && !layers.isEmpty()) {
            const AtmosphericLayerState &topLayer = layers.last();
            const double fixedLayerThicknessMeters = topLayer.thicknessMeters();
            if (fixedLayerThicknessMeters > 0.0 &&
                topLayer.temperatureKelvin() > 0.0 &&
                topLayer.pressureAtm() > 0.0 &&
                gravityMps2_ > 0.0 &&
                rSpecific_ > 0.0) {
                // Масштабная высота: H = R_specific * T / g, тогда
                // P(z + dz) = P(z) * exp(-dz / H). Для нового слоя берём dz = толщину слоя.
                // Толщины фиксируем, чтобы новые слои добавлялись без перерасчёта геометрии.
                const double scaleHeightMeters =
                    (rSpecific_ * topLayer.temperatureKelvin()) / gravityMps2_;
                if (scaleHeightMeters > 0.0) {
                    const double nextLayerPressureAtm =
                        topLayer.pressureAtm() *
                        qExp(-fixedLayerThicknessMeters / scaleHeightMeters);
                    if (nextLayerPressureAtm > minTopPressureAtm) {
                        resizeRequests[i] = {layers.size() + 1, fixedLayerThicknessMeters};
                    } else if (topLayer.pressureAtm() < minTopPressureAtm &&
                               layers.size() > 1) {
                        resizeRequests[i] = {layers.size() - 1, fixedLayerThicknessMeters};
                    }
                }
            }
        }
    }); // конец параллельной фазы 1

    // ── Фаза 1.5: отложенное обновление числа слоёв (последовательно) ──
    // updateColumnLayerCountFixedThickness обновляет общий layerCount_,
    // поэтому вызываем строго последовательно.
    for (int i = 0; i < processedCount; ++i) {
        const auto &request = resizeRequests.at(i);
        if (request.newLayerCount > 0) {
            input.atmosphereGrid.updateColumnLayerCountFixedThickness(
                i, request.newLayerCount, request.layerThicknessMeters);
        }
    }

    // Логирование профиля для выбранной точки.
    {
        const auto &logLayers = input.atmosphereGrid.columns()[logPointIndex].layers();
        qCInfo(atmosphereProfileLog) << "Atmosphere profile (layered step)"
                                     << "index=" << logPointIndex
                                     << "layerCount=" << logLayers.size();
        for (int layerIndex = 0; layerIndex < logLayers.size(); ++layerIndex) {
            const AtmosphericLayerState &layer = logLayers.at(layerIndex);
            qCInfo(atmosphereProfileLog) << "Layer"
                                         << layerIndex
                                         << "heightKm=" << layer.heightMeters() / 1000.0
                                         << "temperatureK=" << layer.temperatureKelvin()
                                         << "pressureAtm=" << layer.pressureAtm();
        }
    }

    // ── Фаза 2: растительность (последовательно — диффузия по соседям) ──
    vegetationModel_.update(input.surfaceGrid.points(), timeStepSeconds_, co2Share_);

    // ── Фаза 2.5: океанские течения и горизонтальный перенос тепла ──────
    oceanCurrentSolver_.step(input.surfaceGrid,
                              dayLengthSeconds_,
                              isRetrograde_,
                              timeStepSeconds_);

    // ── Фаза 2.7: сток поверхностных вод (последовательно — соседи) ─────
    // Вода стекает от каждой точки к наиболее низкому соседу. Если сосед —
    // жидкий океан, вода удаляется (возвращается в океан). Это создаёт
    // условия для жизни в долинах и низинах, куда собирается вода.
    {
        constexpr double kDrainageFraction = 0.15;     // доля стока за шаг
        constexpr double kFieldCapacityKgPerM2 = 200.0; // полевая влагоёмкость
        auto &points = input.surfaceGrid.points();
        const int n = points.size();
        // Буфер изменений: delta[i] — сколько воды добавить/убрать из точки i.
        QVector<double> drainageDelta(n, 0.0);
        for (int i = 0; i < n; ++i) {
            auto &pt = points[i];
            // Жидкий океан не содержит «поверхностной влаги» — пропускаем.
            if (pt.materialId == QLatin1String("ocean") &&
                pt.waterPhase == PhaseModel::Phase::Liquid) {
                continue;
            }
            const double water = pt.surfaceMoisture.waterKgPerM2();
            // Стекает только излишек сверх полевой влагоёмкости.
            const double excess = water - kFieldCapacityKgPerM2;
            if (excess <= 0.0) {
                continue;
            }
            // Найти самого низкого соседа.
            int lowestIdx = -1;
            double lowestHeight = pt.heightKm;
            for (int ni : pt.neighborIndices) {
                if (ni >= 0 && ni < n && points.at(ni).heightKm < lowestHeight) {
                    lowestHeight = points.at(ni).heightKm;
                    lowestIdx = ni;
                }
            }
            if (lowestIdx < 0) {
                continue; // нет соседей ниже — вода стоит
            }
            const double drainAmount = excess * kDrainageFraction;
            drainageDelta[i] -= drainAmount;
            // Если сосед — жидкий океан, вода просто исчезает (в океан).
            const auto &neighbor = points.at(lowestIdx);
            if (neighbor.materialId == QLatin1String("ocean") &&
                neighbor.waterPhase == PhaseModel::Phase::Liquid) {
                // вода уходит в океан — не добавляем соседу
            } else {
                drainageDelta[lowestIdx] += drainAmount;
            }
        }
        // Применяем все изменения одним проходом.
        for (int i = 0; i < n; ++i) {
            if (drainageDelta.at(i) != 0.0) {
                const double newWater = qMax(0.0,
                    points[i].surfaceMoisture.waterKgPerM2() + drainageDelta.at(i));
                points[i].surfaceMoisture.setWaterKgPerM2(newWater);
            }
        }
    }

    // ── Фаза 3: динамика ветра (соседние градиенты) ─────────────────────
    dynamicsSolver_.updateLayerWinds(input.surfaceGrid,
                                    input.atmosphereGrid,
                                    dayLengthSeconds_,
                                    isRetrograde_,
                                    timeStepSeconds_,
                                    1);

    // ── Фаза 4: вертикальное перемешивание ветра (параллельно) ──────────
    // Каждая колонка обрабатывается независимо, mix() — const-метод.
    // Принудительный detach: вектор колонок мог перестроиться в фазе 1.5.
    AtmosphericColumn *columnsForMixing = input.atmosphereGrid.columns().data();
    QtConcurrent::blockingMap(indices, [&, columnsForMixing](int i) {
        verticalWindMixingSolver_.mix(columnsForMixing[i], timeStepSeconds_);
    });

    // ── Фаза 5: горизонтальная адвекция ─────────────────────────────────
    advectionSolver_.advectLayerWinds(input.surfaceGrid,
                                      input.atmosphereGrid,
                                      dayLengthSeconds_,
                                      timeStepSeconds_,
                                      1);
    advectionSolver_.advectLayerMoisture(input.surfaceGrid,
                                         input.atmosphereGrid,
                                         timeStepSeconds_);
}

// ── Разложение газов и атмосферное убегание ────────────────────────────
bool AtmosphereStepSolver::applyGasChemistryAndEscape(
    AtmosphereComposition &composition,
    double timeStepSeconds,
    double planetMassEarths,
    double planetRadiusKm,
    double exosphericTemperatureK) {
    if (timeStepSeconds <= 0.0) return false;

    bool changed = false;

    // ── 1. Химическое разложение (экспоненциальный распад) ─────────────
    struct GasDecay {
        const char *id;
        double lifetimeSeconds;
        const char *productId; // nullptr = удаляется, "co2" = конвертируется
        double productMassRatio; // масса продукта / масса исходного
    };
    const GasDecay decays[] = {
        // CH₄ + 2O₂ → CO₂ + 2H₂O: 16 г → 44 г CO₂ (коэффициент 44/16 ≈ 2.75).
        {"ch4", kCH4LifetimeSeconds, "co2", 44.0 / 16.0},
        // NH₃ → осаждение (удаляется из атмосферы). 4NH₃+3O₂ → 2N₂+6H₂O.
        {"nh3", kNH3LifetimeSeconds, nullptr, 0.0},
        // SF₆ → крайне медленный распад.
        {"sf6", kSF6LifetimeSeconds, nullptr, 0.0},
        // NF₃ → медленный распад.
        {"nf3", kNF3LifetimeSeconds, nullptr, 0.0},
        // H₂ → реакция с OH (без значимых продуктов в модели).
        {"h2",  kH2LifetimeSeconds,  nullptr, 0.0},
    };
    for (const auto &decay : decays) {
        const QString gasId = QLatin1String(decay.id);
        const double currentMass = composition.massGigatons(gasId);
        if (currentMass <= 1.0e-12) continue;

        const double decayFactor = std::exp(-timeStepSeconds / decay.lifetimeSeconds);
        const double newMass = currentMass * decayFactor;
        const double decayedMass = currentMass - newMass;
        if (decayedMass < 1.0e-15) continue;

        composition.setMassGigatons(gasId, qMax(0.0, newMass));
        changed = true;

        // Продукт разложения (если есть).
        if (decay.productId) {
            const QString productId = QLatin1String(decay.productId);
            const double productMass = decayedMass * decay.productMassRatio;
            composition.setMassGigatons(
                productId, composition.massGigatons(productId) + productMass);
        }
    }

    // ── 2. Атмосферное убегание (Jeans escape) ────────────────────────
    if (planetMassEarths <= 0.0 || planetRadiusKm <= 0.0 ||
        exosphericTemperatureK <= 0.0) {
        return changed;
    }

    const double radiusM = planetRadiusKm * 1000.0;
    const double planetMassKg = planetMassEarths * kEarthMass;
    // Экзобаза ≈ планетарный радиус + 500 км (упрощение).
    const double exobaseRadiusM = radiusM + 500000.0;
    const double surfaceAreaM2 = 4.0 * M_PI * exobaseRadiusM * exobaseRadiusM;
    const double totalMassKg = composition.totalMassKg();
    if (totalMassKg <= 0.0) return changed;

    struct EscapeGas {
        const char *id;
        double molarMass; // г/моль ≈ а.е.м.
    };
    // Только лёгкие газы существенно убегают.
    const EscapeGas escapeGases[] = {
        {"h2",  kMolarMassH2},
        {"he",  kMolarMassHe},
        {"ch4", kMolarMassCH4},
        {"nh3", kMolarMassNH3},
    };
    for (const auto &gas : escapeGases) {
        const QString gasId = QLatin1String(gas.id);
        const double massFractionGt = composition.massGigatons(gasId);
        if (massFractionGt <= 1.0e-12) continue;

        const double molecularMassKg = gas.molarMass * kAmuToKg;
        // λ = G·M·m / (r_exo · k_B · T_exo) — параметр убегания.
        const double lambda = kGravConst * planetMassKg * molecularMassKg /
                              (exobaseRadiusM * kBoltzmannConstant * exosphericTemperatureK);
        // Для λ > 40 убегание пренебрежимо мало (exp(-40) ≈ 4e-18).
        if (lambda > 40.0) continue;

        // Тепловая скорость: v_th = sqrt(2 k T / m).
        const double vThermal = std::sqrt(
            2.0 * kBoltzmannConstant * exosphericTemperatureK / molecularMassKg);
        // Число плотность на экзобазе: n ≈ P / (k T) ≈ (m_gas·g)/(A·kT).
        // Упрощение: используем долю газа в общей массе.
        const double gasMassKg = massFractionGt * 1.0e12; // Гт → кг
        // Средняя плотность на экзобазе (грубая оценка из баро-формулы
        // с одной шкалой высот; точность не критична — важен порядок).
        const double scaleHeight = kBoltzmannConstant * exosphericTemperatureK /
                                   (molecularMassKg * kGravConst * planetMassKg /
                                    (exobaseRadiusM * exobaseRadiusM));
        const double nExobase = gasMassKg / (molecularMassKg * surfaceAreaM2 * scaleHeight);

        // Поток Джинса: Φ = n · v_th/(2√π) · (1 + λ) · exp(−λ).
        const double flux = nExobase * vThermal / (2.0 * std::sqrt(M_PI)) *
                            (1.0 + lambda) * std::exp(-lambda);
        // Потеря массы за шаг (кг).
        const double massLossKg = flux * molecularMassKg * surfaceAreaM2 * timeStepSeconds;
        const double massLossGt = massLossKg * 1.0e-12;
        if (massLossGt < 1.0e-15) continue;

        const double newMass = qMax(0.0, massFractionGt - massLossGt);
        composition.setMassGigatons(gasId, newMass);
        changed = true;
    }

    return changed;
}
