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
// Массовый коэффициент поглощения водяного пара в длинноволновом (ИК) диапазоне.
// Типичные значения для полосно-осреднённого H₂O: 0.08–0.3 м²/кг.
// При 25 кг/м² (средний столбик) даёт τ ≈ 2.5 — реалистичный парниковый эффект.
constexpr double kH2OMassAbsorptionLw = 0.1;
// Поглощение H₂O в коротковолновом ближнем ИК (слабее, чем длинноволновое).
constexpr double kH2OMassAbsorptionSw = 0.01;
// Доля сублимации: ледяной океан обеспечивает ~10% влаги от жидкого.
// В реальности сублимация льда при -10°C ≈ 10–15% от испарения воды при 0°C.
constexpr double kIceSublimationFraction = 0.1;
// Массовый коэффициент поглощения CO₂ в длинноволновом ИК (15-мкм полоса).
// При 1 атм уширение линий давлением включено: k_eff = k_base × P.
// Для Земли (1 атм): k × P × co2_column ≈ 0.35 × 1 × 3 ≈ 1.05 → реалистично.
constexpr double kCO2MassAbsorptionLw = 0.35;
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
                                           bool isRetrograde)
    : radiationSolver_(timeStepSeconds)
    , convectiveSolver_(composition, gravityMps2)
    , gravityMps2_(gravityMps2)
    , rSpecific_(AtmosphericThermodynamics::specificGasConstant(composition))
    , specificHeatCp_(AtmosphericThermodynamics::specificHeatCp(composition))
    , timeStepSeconds_(timeStepSeconds)
    , dayLengthSeconds_(dayLengthSeconds)
    , isRetrograde_(isRetrograde)
    , co2Share_(gasShareById(composition, QStringLiteral("co2"))) {}

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
            const double phaseTemperatureK =
                (point.temperatureK > 0.0) ? point.temperatureK : point.airTemperatureK;
            // Упрощённый фазовый порог: ниже 273.15 K осадки считаем снегом.
            const bool isSnow =
                (phaseTemperatureK > 0.0 && phaseTemperatureK <= 273.15) &&
                point.materialId != QLatin1String("ocean");
            if (isSnow) {
                point.snowKgPerM2 += precipitationKgPerM2;
            } else {
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

        // Динамический парниковый эффект H₂O и CO₂: добавляем вклад обоих
        // в оптическую толщину каждого слоя перед расчётом радиации.
        //
        // H₂O — переменный: больше влаги → сильнее парник → теплее поверхность
        // → больше испарения → положительная обратная связь.
        //
        // CO₂ — фиксированный, но с уширением линий давлением (pressure broadening).
        // Базовая τ при инициализации слишком мала для земных условий (0.06% CO₂
        // при 1 атм → τ ≈ 0.05). Реальный CO₂ даёт τ ≈ 1.0 из-за уширения.
        // Добавляем поправку: k_co2 × co2_fraction × density × thickness × P.
        for (int li = 0; li < layers.size(); ++li) {
            const double h2oKg = qMax(0.0, layers.at(li).waterVaporKgPerM2());
            // CO₂ column mass in layer (kg/m²) × pressure broadening factor.
            const double co2ColumnKg = co2Share_ *
                qMax(0.0, layers.at(li).densityKgPerM3()) *
                qMax(0.0, layers.at(li).thicknessMeters());
            const double pressureBroadening =
                qMax(0.0, layers.at(li).pressureAtm());
            const double co2TauLw =
                kCO2MassAbsorptionLw * co2ColumnKg * pressureBroadening;
            layers[li].setOpticalDepthLongwave(
                layers.at(li).opticalDepthLongwave() +
                kH2OMassAbsorptionLw * h2oKg + co2TauLw);
            layers[li].setOpticalDepthShortwave(
                layers.at(li).opticalDepthShortwave() + kH2OMassAbsorptionSw * h2oKg);
        }

        // Поверхность и нижний слой — разные сущности: передаём температуру поверхности явно,
        // чтобы не подменять её температурой слоя и не получать самоподогрев атмосферы.
        LayeredRadiationSolver::SurfaceRadiativeFluxes surfaceRadFluxes;
        const QVector<double> layerDeltas =
            radiationSolver_.solve(column,
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

        // Восстанавливаем базовую оптическую толщину (без H₂O и CO₂-поправки),
        // чтобы на следующем шаге не накапливать вклад повторно.
        for (int li = 0; li < layers.size(); ++li) {
            const double h2oKg = qMax(0.0, layers.at(li).waterVaporKgPerM2());
            const double co2ColumnKg = co2Share_ *
                qMax(0.0, layers.at(li).densityKgPerM3()) *
                qMax(0.0, layers.at(li).thicknessMeters());
            const double pressureBroadening =
                qMax(0.0, layers.at(li).pressureAtm());
            const double co2TauLw =
                kCO2MassAbsorptionLw * co2ColumnKg * pressureBroadening;
            layers[li].setOpticalDepthLongwave(
                layers.at(li).opticalDepthLongwave() -
                kH2OMassAbsorptionLw * h2oKg - co2TauLw);
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

        SurfaceAtmosphereCoupler coupler(input.heatTransferCoefficientWPerM2K);
        double surfaceAirFluxWPerM2 = 0.0;
        coupler.exchangeHeat(point.state,
                             layers[0],
                             material.roughnessLengthMeters,
                             timeStepSeconds_,
                             &surfaceAirFluxWPerM2);
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
