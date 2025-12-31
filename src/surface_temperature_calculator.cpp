#include "surface_temperature_calculator.h"

#include "atmospheric_pressure_model.h"
#include "surface_height_model.h"
#include "surface_temperature_state.h"
#include "surface_atmosphere_coupler.h"

#include <QtCore/QThread>
#include <QtCore/QtMath>

#include <algorithm>
#include <array>
#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kKelvinOffset = 273.15;
constexpr double kPi = 3.14159265358979323846;
constexpr double kEarthRadiusKm = 6371.0;
constexpr double kEarthAreaKm2 = 510072000.0;
constexpr double kEarthAtmosphereMassGt = 5140000.0;
constexpr double kDefaultSurfaceRoughness = 20.0;
constexpr double kDefaultBasinShape = 3.5;
constexpr double kEarthWaterGigatons = 1.4e9;
constexpr int kDailyTimeSteps = 96;
constexpr int kSpinUpDays = 6;
constexpr double kDryAirSpecificHeatJPerKgK = 1004.0;
constexpr double kEarthGravityMPerS2 = 9.80665;
constexpr double kStandardPressurePa = 101325.0;
constexpr double kDefaultHeatTransferWPerM2K = 8.0;

struct TraceGasSpec {
    const char *id;
    double greenhousePower;
};

constexpr std::array<TraceGasSpec, 4> kTraceGases = {{
    {"ch4", 0.5},
    {"nh3", 1.5},
    {"sf6", 300.0},
    {"nf3", 250.0},
}};


double gasMassGigatons(const AtmosphereComposition &atmosphere, const QString &gasId) {
    return qMax(0.0, atmosphere.massGigatons(gasId));
}

double estimateSurfaceWaterGigatons(const SurfaceMaterial &material) {
    // В текущем UI нет явного управления гидросферой, поэтому применяем мягкую эвристику:
    // океаническую поверхность считаем водной, остальные материалы — сухими.
    if (material.id == QLatin1String("ocean")) {
        return kEarthWaterGigatons;
    }
    return 0.0;
}
}  // namespace

SurfaceTemperatureCalculator::SurfaceTemperatureCalculator(double solarConstant,
                                                           const SurfaceMaterial &material,
                                                           double dayLengthDays,
                                                           RotationMode rotationMode,
                                                           const AtmosphereComposition &atmosphere,
                                                           double greenhouseOpacity,
                                                           double presetCloudAlbedo,
                                                           double atmospherePressureAtm,
                                                           double surfaceGravity,
                                                           double planetRadiusKm,
                                                           bool useAtmosphericModel,
                                                           int meridionalTransportSteps,
                                                           HeightSourceType heightSourceType,
                                                           const QString &heightmapPath,
                                                           double heightmapScaleKm,
                                                           quint32 heightSeed,
                                                           bool useContinentsHeight,
                                                           const SubsurfaceModelSettings &subsurfaceSettings)
    : solarConstant_(solarConstant),
      material_(material),
      dayLengthDays_(dayLengthDays),
      rotationMode_(rotationMode),
      atmosphere_(atmosphere),
      greenhouseOpacity_(qBound(0.0, greenhouseOpacity, 0.999)),
      presetCloudAlbedo_(qBound(0.0, presetCloudAlbedo, 1.0)),
      atmospherePressureAtm_(atmospherePressureAtm),
      surfaceGravity_(surfaceGravity),
      planetRadiusKm_(planetRadiusKm),
      useAtmosphericModel_(useAtmosphericModel),
      meridionalTransportSteps_(qMax(1, meridionalTransportSteps)),
      heightSourceType_(heightSourceType),
      heightmapPath_(heightmapPath),
      heightmapScaleKm_(heightmapScaleKm),
      heightSeed_(heightSeed),
      useContinentsHeight_(useContinentsHeight),
      subsurfaceSettings_(subsurfaceSettings) {}

void SurfaceTemperatureCalculator::setMeridionalTransportSteps(int steps) {
    meridionalTransportSteps_ = qMax(1, steps);
}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesByLatitude(
    int latitudePoints) const {
    return temperatureRangesByLatitude(latitudePoints, ProgressCallback{}, nullptr);
}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesByLatitude(
    int latitudePoints,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag) const {
    const int totalLatitudes = latitudePoints > 1 ? latitudePoints : 0;
    return temperatureRangesByLatitudeForSegment(latitudePoints, solarConstant_, 0.0,
                                                 progressCallback, cancelFlag, nullptr, 0,
                                                 totalLatitudes);
}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesForOrbitSegment(
    const OrbitSegment &segment,
    double referenceDistanceAU,
    double obliquityDegrees,
    double perihelionArgumentDegrees,
    int latitudePoints,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag,
    const std::atomic_bool *pauseFlag) const {
    if (cancelFlag && cancelFlag->load()) {
        return {};
    }
    if (latitudePoints <= 1) {
        return {};
    }

    const double obliquityRadians = qDegreesToRadians(obliquityDegrees);
    const double perihelionArgumentRadians = qDegreesToRadians(perihelionArgumentDegrees);
    // Сезонная деклинация: угол между лучами звезды и экватором планеты.
    // δ = asin(sin(наклон оси) * sin(истинная долгота звезды)).
    const double solarLongitude = segment.trueAnomalyRadians + perihelionArgumentRadians;
    const double declinationDegrees = qRadiansToDegrees(
        std::asin(std::sin(obliquityRadians) * std::sin(solarLongitude)));
    // Инсоляция меняется с расстоянием как 1 / r^2 относительно опорной дистанции.
    const double segmentSolarConstant =
        solarConstant_ * std::pow(referenceDistanceAU / segment.distanceAU, 2.0);
    const int totalLatitudes = latitudePoints;

    return temperatureRangesByLatitudeForSegment(latitudePoints,
                                                 segmentSolarConstant,
                                                 declinationDegrees,
                                                 progressCallback,
                                                 cancelFlag,
                                                 pauseFlag,
                                                 0,
                                                 totalLatitudes);
}

QVector<QVector<TemperatureRangePoint>> SurfaceTemperatureCalculator::temperatureRangesByOrbitSegments(
    const QVector<OrbitSegment> &segments,
    double referenceDistanceAU,
    double obliquityDegrees,
    double perihelionArgumentDegrees,
    int latitudePoints,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag) const {
    QVector<QVector<TemperatureRangePoint>> results;
    if (cancelFlag && cancelFlag->load()) {
        return results;
    }
    if (segments.isEmpty() || latitudePoints <= 1) {
        return results;
    }

    const double obliquityRadians = qDegreesToRadians(obliquityDegrees);
    const double perihelionArgumentRadians = qDegreesToRadians(perihelionArgumentDegrees);
    const int latitudesCount = latitudePoints;
    const int totalProgress = latitudesCount * segments.size();

    results.reserve(segments.size());
    for (int i = 0; i < segments.size(); ++i) {
        if (cancelFlag && cancelFlag->load()) {
            return {};
        }

        const auto &segment = segments.at(i);
        // Сезонная деклинация: угол между лучами звезды и экватором планеты.
        // δ = asin(sin(наклон оси) * sin(истинная долгота звезды)).
        const double solarLongitude = segment.trueAnomalyRadians + perihelionArgumentRadians;
        const double declinationDegrees = qRadiansToDegrees(
            std::asin(std::sin(obliquityRadians) * std::sin(solarLongitude)));
        // Инсоляция меняется с расстоянием как 1 / r^2 относительно опорной дистанции.
        const double segmentSolarConstant =
            solarConstant_ * std::pow(referenceDistanceAU / segment.distanceAU, 2.0);
        const int progressOffset = i * latitudesCount;

        results.push_back(temperatureRangesByLatitudeForSegment(latitudePoints,
                                                                segmentSolarConstant,
                                                                declinationDegrees,
                                                                progressCallback,
                                                                cancelFlag,
                                                                nullptr,
                                                                progressOffset,
                                                                totalProgress));
    }

    return results;
}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesByLatitudeForSegment(
    int latitudePoints,
    double segmentSolarConstant,
    double declinationDegrees,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag,
    const std::atomic_bool *pauseFlag,
    int progressOffset,
    int totalProgress) const {
    QVector<TemperatureRangePoint> points =
        radiativeBalanceByLatitudeForSegment(latitudePoints,
                                             segmentSolarConstant,
                                             declinationDegrees,
                                             progressCallback,
                                             cancelFlag,
                                             pauseFlag,
                                             progressOffset,
                                             totalProgress);
    // Текущий расчёт атмосферной поправки (аэродинамика, адиабатический градиент)
    // закомментирован по запросу: ниже применяется модель из React-кода.
#if 0
    if (!useAtmosphericModel_ || points.isEmpty()) {
        // Для безатмосферного режима возвращаем чисто радиационно-кондуктивный баланс.
        return points;
    }

    // Вертикальная поправка: корректируем температуру по адиабатическому градиенту
    // сразу после радиационного баланса, перед горизонтальной циркуляцией.
    AtmosphericLapseRateModel lapseRateModel(atmospherePressureAtm_, atmosphere_,
                                             surfaceGravity_);
    QVector<TemperatureRangePoint> adjusted;
    adjusted.reserve(points.size());
    for (const auto &point : points) {
        adjusted.push_back(lapseRateModel.applyLapseRate(point));
    }

    AtmosphericCirculationModel circulationModel(dayLengthDays_, rotationMode_,
                                                  atmospherePressureAtm_);
    circulationModel.setAtmosphereMassKg(atmosphere_.totalMassKg());
    circulationModel.setMeridionalTransportSteps(meridionalTransportSteps_);
    return circulationModel.applyHeatTransport(adjusted);
#endif
    return points;
}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::radiativeBalanceByLatitudeForSegment(
    int latitudePoints,
    double segmentSolarConstant,
    double declinationDegrees,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag,
    const std::atomic_bool *pauseFlag,
    int progressOffset,
    int totalProgress) const {
    QVector<TemperatureRangePoint> points;
    if (latitudePoints <= 1) {
        return points;
    }
    if (cancelFlag && cancelFlag->load()) {
        return {};
    }

    points.reserve(latitudePoints);

    int processedLatitudes = 0;
    const double declinationRadians = qDegreesToRadians(declinationDegrees);
    const double stepDegrees = 180.0 / static_cast<double>(latitudePoints - 1);
    const double safeRadiusKm = qMax(0.1, planetRadiusKm_);
    const double areaScale = std::pow(safeRadiusKm / kEarthRadiusKm, 2.0);
    const double surfaceGravity = qMax(0.0, surfaceGravity_);
    const double totalGas = atmosphere_.totalMassGigatons();
    const double seaLevelPressureAtm =
        (totalGas > 0.0)
            ? (totalGas / kEarthAtmosphereMassGt) * (surfaceGravity / 9.8) / areaScale
            : 0.0;
    const double co2Mass = gasMassGigatons(atmosphere_, QStringLiteral("co2"));
    const double waterGigatons = estimateSurfaceWaterGigatons(material_);
    const double planetAreaKm2 = kEarthAreaKm2 * areaScale;
    const double avgDepth = (planetAreaKm2 > 0.0) ? waterGigatons / planetAreaKm2 : 0.0;
    double potentialCoverage = 0.0;
    if (waterGigatons > 0.0) {
        const double fillFactor = (avgDepth * kDefaultBasinShape) / kDefaultSurfaceRoughness;
        potentialCoverage = 1.0 - std::exp(-fillFactor * 3.0);
    }
    potentialCoverage = qBound(0.0, potentialCoverage, 1.0);

    const double colDensity = 1.0 / qMax(1.0, areaScale);
    const double albedo = qBound(0.0, material_.albedo, 1.0);
    const double baseRadiativeTemp =
        std::pow((segmentSolarConstant * (1.0 - albedo)) / (4.0 * kStefanBoltzmannConstant),
                 0.25);
    const SurfaceHeightModel heightModel(heightSourceType_, heightmapPath_, heightmapScaleKm_,
                                         heightSeed_, useContinentsHeight_);

    const bool isTidallyLocked = rotationMode_ == RotationMode::TidalLocked;
    // Подбираем число шагов по длительности суток, чтобы шаг по времени не разрастался
    // на медленно вращающихся планетах (например, Меркурий).
    const int scaledStepsPerDay = qRound(dayLengthDays_ * 24.0);
    const int stepsPerDay = qMax(kDailyTimeSteps, scaledStepsPerDay);
    const int spinUpDays = isTidallyLocked ? 2 : kSpinUpDays;
    const double dayLengthSeconds = qMax(0.01, dayLengthDays_) * 86400.0;
    const double timeStepSeconds =
        (stepsPerDay > 0) ? (dayLengthSeconds / static_cast<double>(stepsPerDay)) : 0.0;
    // Глобальный средний поток перед альбедо: отражение применяется в SurfaceTemperatureState,
    // чтобы избежать повторного учета в суточной итерации.
    const double globalAverageInsolation = segmentSolarConstant / 4.0;

    for (int i = 0; i < latitudePoints; ++i) {
        if (cancelFlag && cancelFlag->load()) {
            return {};
        }
        if (pauseFlag) {
            // Пауза из интерфейса: удерживаем вычисления, но проверяем отмену.
            while (pauseFlag->load()) {
                if (cancelFlag && cancelFlag->load()) {
                    return {};
                }
                QThread::msleep(50);
            }
        }

        const double axisDegrees = (rotationMode_ == RotationMode::Normal)
                                       ? (-90.0 + stepDegrees * static_cast<double>(i))
                                       : (stepDegrees * static_cast<double>(i));
        const double latitudeRadians = qDegreesToRadians(axisDegrees);
        const double tanVal = -std::tan(latitudeRadians) * std::tan(declinationRadians);
        double hourAngleLimit = 0.0;
        if (tanVal >= 1.0) {
            hourAngleLimit = 0.0;
        } else if (tanVal <= -1.0) {
            hourAngleLimit = kPi;
        } else {
            hourAngleLimit = std::acos(tanVal);
        }
        // Средний дневной поток повторяет реактовскую формулу через H0
        // (угол захода/восхода) и учитывает полярные дни/ночи.
        const double dailyFactor =
            (hourAngleLimit * std::sin(latitudeRadians) * std::sin(declinationRadians) +
             std::cos(latitudeRadians) * std::cos(declinationRadians) *
                 std::sin(hourAngleLimit)) /
            kPi;
        const bool hasInsolation = dailyFactor > 0.0;

        // Высоту берём из карты рельефа на меридиане 0°, переводя километры в метры.
        const double heightMeters = heightModel.heightKmAt(axisDegrees, 0.0) * 1000.0;
        // Давление по барометрической формуле рассчитываем от уровня моря с учётом высоты.
        const double pressureAtm = AtmosphericPressureModel::pressureAtHeightAtm(
            seaLevelPressureAtm,
            heightMeters,
            qMax(1.0, baseRadiativeTemp),
            atmosphere_,
            surfaceGravity);

        // Модель оптической толщины повторяет формулы из React-кода:
        // tau_CO2 = ln(1 + M_CO2 * colDensity * power * 0.001) * broadening,
        // где broadening учитывает расширение линий при росте давления.
        const double broadening = std::pow(qMax(0.5, pressureAtm * 2.0), 0.32);
        const double baseTau =
            std::log(1.0 + co2Mass * colDensity * 0.02 * 0.001) * broadening;
        double traceTau = 0.0;
        for (const auto &spec : kTraceGases) {
            const double traceMass = gasMassGigatons(atmosphere_, QLatin1String(spec.id));
            traceTau += (traceMass / qMax(1.0, areaScale)) * spec.greenhousePower * 1e-6 *
                        broadening;
        }

        double pressureClouds =
            pressureAtm > 0.05 ? 0.25 * (1.0 - std::exp(-pressureAtm)) : 0.0;
        const double surfAlbedoPre =
            (1.0 - potentialCoverage) * albedo + potentialCoverage * 0.06;
        const double tEffPre =
            std::pow((segmentSolarConstant * (1.0 - qMax(surfAlbedoPre, pressureClouds))) /
                         (4.0 * kStefanBoltzmannConstant),
                     0.25);
        const double tBasePre =
            tEffPre * std::pow(1.0 + 0.75 * (baseTau + traceTau), 0.25);

        double evaporation = 0.0;
        if (potentialCoverage > 0.0 && tBasePre > 263.0) {
            evaporation = potentialCoverage * std::exp((tBasePre - 280.0) / 15.0);
        }
        const double waterClouds = qMin(0.5, evaporation * 0.28);
        // Сернокислотные облака (H₂SO₄) Венеры отражают больше, чем водные облака,
        // поэтому учитываем минимум по пресету, даже если испарение воды отсутствует.
        const double dynamicCloudAlbedo = qMin(0.88, pressureClouds + waterClouds);
        const double presetCloudAlbedo = qBound(0.0, presetCloudAlbedo_, 0.88);
        const double cloudAlbedo = qMax(dynamicCloudAlbedo, presetCloudAlbedo);
        const double waterTau = qMin(8.0, evaporation * 1.5);
        double totalTau = baseTau + traceTau + waterTau;
        if (!useAtmosphericModel_ && greenhouseOpacity_ > 0.0) {
            // Дополнительная непрозрачность, когда атмосферная модель выключена,
            // но задана парниковая "шторка" вручную.
            totalTau += -std::log(qMax(1e-6, 1.0 - greenhouseOpacity_));
        }
        const double ghMult = std::pow(1.0 + 0.75 * totalTau, 0.25);
        // Приводим оптическую толщину к коэффициенту парникового эффекта для модели
        // SurfaceTemperatureState: в стационаре T^4 ∝ 1 / (1 - G).
        const double greenhouseOpacity =
            (totalTau > 0.0) ? (1.0 - 1.0 / (1.0 + 0.75 * totalTau)) : 0.0;

        const double transport =
            (pressureAtm > 50.0)
                ? 0.99
                : (pressureAtm > 0.001
                       ? qMin(1.0, 0.15 * std::log(pressureAtm * 100.0 + 1.0))
                       : 0.0);
        const double rotBlock =
            (dayLengthDays_ < 2.0 && pressureAtm < 10.0) ? 0.65 : 1.0;
        const double meridionalTransport = transport * rotBlock;

        double dynamicSurfAlbedo = albedo * (1.0 - potentialCoverage);
        if (tBasePre < 260.0) {
            dynamicSurfAlbedo += potentialCoverage * 0.70;
        } else {
            dynamicSurfAlbedo += potentialCoverage * 0.06;
        }
        const double finalAlbedo = qMax(dynamicSurfAlbedo, cloudAlbedo);
        const double tEff =
            std::pow((segmentSolarConstant * (1.0 - finalAlbedo)) /
                         (4.0 * kStefanBoltzmannConstant),
                     0.25);
        const double tGlobalAvg = tEff * ghMult;

        const double tLatRad =
            std::pow(qMax(0.1,
                          segmentSolarConstant * (1.0 - finalAlbedo) * dailyFactor) /
                         kStefanBoltzmannConstant,
                     0.25) *
            ghMult;
        const double tBase =
            tLatRad * (1.0 - meridionalTransport) + tGlobalAvg * meridionalTransport;

        // В базовом радиационном балансе не ограничиваем температуру атмосферным "полом":
        // слабые звезды должны давать холодные поверхности вплоть до физического минимума ~3 K.
        const double absFloor = 3.0;

        // Шаговое обновление температуры на протяжении суток:
        // S_inst = S0 * cos(zenith) при освещении, иначе 0.
        // Поток смешивается с глобальным средним, чтобы имитировать меридиональный перенос.
        SurfaceTemperatureState state(tBase,
                                      finalAlbedo,
                                      greenhouseOpacity,
                                      absFloor,
                                      material_,
                                      subsurfaceSettings_);
        const double gravity = (surfaceGravity_ > 0.0) ? surfaceGravity_ : kEarthGravityMPerS2;
        const double columnMassKgPerM2 =
            (pressureAtm > 0.0) ? (pressureAtm * kStandardPressurePa / gravity) : 0.0;
        const double airHeatCapacity = columnMassKgPerM2 * kDryAirSpecificHeatJPerKgK;
        AtmosphericCellState airState(tBase, airHeatCapacity);
        // Коэффициент теплообмена h_c в Вт/(м^2·К) связывает поверхность и воздух.
        // В реальности он зависит от ветра, давления и турбулентности; здесь масштабируем
        // от давления, чтобы тонкая атмосфера слабее влияла на поверхность.
        const double couplingScale = qBound(0.0, pressureAtm, 1.0);
        SurfaceAtmosphereCoupler coupler(kDefaultHeatTransferWPerM2K * couplingScale);

        const int totalSteps = stepsPerDay * (spinUpDays + 1);
        double tMin = state.temperatureKelvin();
        double tMax = state.temperatureKelvin();
        double tSum = 0.0;
        double tDaySum = 0.0;
        double tNightSum = 0.0;
        int dayCount = 0;
        int nightCount = 0;
        for (int step = 0; step < totalSteps; ++step) {
            const double phase =
                (stepsPerDay > 0) ? (2.0 * kPi * (static_cast<double>(step % stepsPerDay) +
                                                 0.5) /
                                        static_cast<double>(stepsPerDay))
                                  : 0.0;
            const double hourAngle = isTidallyLocked ? 0.0 : (phase - kPi);
            const double cosZenith =
                std::sin(latitudeRadians) * std::sin(declinationRadians) +
                std::cos(latitudeRadians) * std::cos(declinationRadians) *
                    std::cos(hourAngle);
            // Суточная инсоляция до учета альбедо.
            const double localInsolation =
                segmentSolarConstant * qMax(0.0, cosZenith);
            const double blendedInsolation =
                localInsolation * (1.0 - meridionalTransport) +
                globalAverageInsolation * meridionalTransport;

            const double absorbedFlux = state.absorbedFlux(blendedInsolation);
            const double emittedFlux = state.emittedFlux();
            state.updateTemperature(absorbedFlux, emittedFlux, timeStepSeconds);
            coupler.exchangeSensibleHeat(state, airState, timeStepSeconds);

            if (step >= stepsPerDay * spinUpDays) {
                const double temp = state.temperatureKelvin();
                tMin = qMin(tMin, temp);
                tMax = qMax(tMax, temp);
                tSum += temp;
                if (cosZenith > 0.0) {
                    tDaySum += temp;
                    ++dayCount;
                } else {
                    tNightSum += temp;
                    ++nightCount;
                }
            }
        }

        const double tMean = (stepsPerDay > 0) ? (tSum / stepsPerDay) : state.temperatureKelvin();
        const double tDay = (dayCount > 0) ? (tDaySum / dayCount) : tMean;
        const double tNight = (nightCount > 0) ? (tNightSum / nightCount) : tMean;

        TemperatureRangePoint point;
        point.latitudeDegrees = axisDegrees;
        point.hasInsolation = hasInsolation;
        point.minimumKelvin = tMin;
        point.maximumKelvin = tMax;
        point.meanDailyKelvin = tMean;
        point.meanDayKelvin = tDay;
        point.meanNightKelvin = tNight;
        point.minimumCelsius = tMin - kKelvinOffset;
        point.maximumCelsius = tMax - kKelvinOffset;
        point.meanDailyCelsius = tMean - kKelvinOffset;
        point.meanDayCelsius = tDay - kKelvinOffset;
        point.meanNightCelsius = tNight - kKelvinOffset;
        points.push_back(point);

        ++processedLatitudes;
        if (progressCallback) {
            progressCallback(progressOffset + processedLatitudes, totalProgress);
        }
    }

    return points;
}
