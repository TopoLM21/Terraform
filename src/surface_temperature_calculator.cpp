#include "surface_temperature_calculator.h"

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
                                                           double atmospherePressureAtm,
                                                           double surfaceGravity,
                                                           double planetRadiusKm,
                                                           bool useAtmosphericModel,
                                                           int meridionalTransportSteps)
    : solarConstant_(solarConstant),
      material_(material),
      dayLengthDays_(dayLengthDays),
      rotationMode_(rotationMode),
      atmosphere_(atmosphere),
      greenhouseOpacity_(qBound(0.0, greenhouseOpacity, 0.999)),
      atmospherePressureAtm_(atmospherePressureAtm),
      surfaceGravity_(surfaceGravity),
      planetRadiusKm_(planetRadiusKm),
      useAtmosphericModel_(useAtmosphericModel),
      meridionalTransportSteps_(qMax(1, meridionalTransportSteps)) {}

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
                                                 progressCallback, cancelFlag, 0, totalLatitudes);
}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesForOrbitSegment(
    const OrbitSegment &segment,
    double referenceDistanceAU,
    double obliquityDegrees,
    double perihelionArgumentDegrees,
    int latitudePoints,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag) const {
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
    int progressOffset,
    int totalProgress) const {
    QVector<TemperatureRangePoint> points =
        radiativeBalanceByLatitudeForSegment(latitudePoints,
                                             segmentSolarConstant,
                                             declinationDegrees,
                                             progressCallback,
                                             cancelFlag,
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
    const double pressureFromMass =
        (totalGas > 0.0)
            ? (totalGas / kEarthAtmosphereMassGt) * (surfaceGravity / 9.8) / areaScale
            : 0.0;
    // Если масса газов не задана, берём давление из UI-настройки атмосферы,
    // чтобы демпфирование и связанные формулы работали даже без явной массы.
    const double pressureAtm =
        (totalGas > 0.0) ? pressureFromMass : qMax(0.0, atmospherePressureAtm_);
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

    // Модель оптической толщины повторяет формулы из React-кода:
    // tau_CO2 = ln(1 + M_CO2 * colDensity * power * 0.001) * broadening,
    // где broadening учитывает расширение линий при росте давления.
    const double colDensity = 1.0 / qMax(1.0, areaScale);
    const double broadening = std::pow(qMax(0.5, pressureAtm * 2.0), 0.32);
    const double baseTau =
        std::log(1.0 + co2Mass * colDensity * 0.02 * 0.001) * broadening;
    double traceTau = 0.0;
    for (const auto &spec : kTraceGases) {
        const double traceMass = gasMassGigatons(atmosphere_, QLatin1String(spec.id));
        traceTau += (traceMass / qMax(1.0, areaScale)) * spec.greenhousePower * 1e-6 *
                    broadening;
    }

    const double albedo = qBound(0.0, material_.albedo, 1.0);
    double pressureClouds = pressureAtm > 0.05 ? 0.25 * (1.0 - std::exp(-pressureAtm)) : 0.0;
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
    const double cloudAlbedo = qMin(0.88, pressureClouds + waterClouds);
    const double waterTau = qMin(8.0, evaporation * 1.5);
    double totalTau = baseTau + traceTau + waterTau;
    if (!useAtmosphericModel_ && greenhouseOpacity_ > 0.0) {
        // Дополнительная непрозрачность, когда атмосферная модель выключена,
        // но задана парниковая "шторка" вручную.
        totalTau += -std::log(qMax(1e-6, 1.0 - greenhouseOpacity_));
    }
    const double ghMult = std::pow(1.0 + 0.75 * totalTau, 0.25);

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

    const bool isTidallyLocked = rotationMode_ == RotationMode::TidalLocked;

    for (int i = 0; i < latitudePoints; ++i) {
        if (cancelFlag && cancelFlag->load()) {
            return {};
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

        const double tLatRad =
            std::pow(qMax(0.1,
                          segmentSolarConstant * (1.0 - finalAlbedo) * dailyFactor) /
                         kStefanBoltzmannConstant,
                     0.25) *
            ghMult;
        const double tBase =
            tLatRad * (1.0 - meridionalTransport) + tGlobalAvg * meridionalTransport;

        // Амплитуда суточного хода переносится из исходного React-кода:
        // A = sqrt(period) * 140 / sqrt(inertia).
        const double atmDamping = pressureAtm * 20.0;
        const double waterInertia = potentialCoverage * 120.0;
        // Суточная инерция - это эффективный (не строго физический) параметр
        // модели, поэтому используем отдельное поле material_.dailyThermalInertia.
        const double totalInertia =
            qMax(1.0, material_.dailyThermalInertia + waterInertia + atmDamping);
        double amplitude = !isTidallyLocked
                               ? (std::sqrt(qMax(0.01, dayLengthDays_)) * 140.0) /
                                     std::sqrt(totalInertia)
                               : 0.0;
        if (hourAngleLimit == 0.0) {
            amplitude = 0.0;
        } else if (hourAngleLimit == kPi) {
            amplitude *= 0.2;
        } else {
            amplitude *= (tBase / 300.0);
        }

        double tMax = tBase + amplitude / 2.0;
        double tMin = tBase - amplitude / 2.0;
        const double absFloor =
            (pressureAtm > 0.5 ? 180.0 : 20.0) + 150.0 * (1.0 - std::exp(-pressureAtm * 0.2));
        if (tMin < absFloor) {
            tMin = absFloor;
        }
        if (tMax < tMin) {
            tMax = tMin;
        }

        TemperatureRangePoint point;
        point.latitudeDegrees = axisDegrees;
        point.hasInsolation = hasInsolation;
        point.minimumKelvin = tMin;
        point.maximumKelvin = tMax;
        point.meanDailyKelvin = tBase;
        point.meanDayKelvin = tMax;
        point.meanNightKelvin = tMin;
        point.minimumCelsius = tMin - kKelvinOffset;
        point.maximumCelsius = tMax - kKelvinOffset;
        point.meanDailyCelsius = tBase - kKelvinOffset;
        point.meanDayCelsius = tMax - kKelvinOffset;
        point.meanNightCelsius = tMin - kKelvinOffset;
        points.push_back(point);

        ++processedLatitudes;
        if (progressCallback) {
            progressCallback(progressOffset + processedLatitudes, totalProgress);
        }
    }

    return points;
}
