#include "surface_temperature_calculator.h"

#include "atmospheric_circulation_model.h"
#include "atmospheric_lapse_rate_model.h"
#include "atmospheric_radiation_model.h"
#include "surface_temperature_state.h"

#include <QtCore/QtMath>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kKelvinOffset = 273.15;
constexpr double kPi = 3.14159265358979323846;
constexpr double kSpaceTemperatureKelvin = 3.0;
constexpr double kSecondsPerEarthDay = 86400.0;
constexpr double kBaseStepSeconds = 900.0;
constexpr int kMinStepsPerDay = 48;
constexpr int kMaxStepsPerDay = 20000;
constexpr double kPascalPerAtm = 101325.0;
constexpr double kReferenceCpDry = 1004.0;
constexpr double kReferenceCpWet = 1850.0;

double meanDailyCosine(double latitudeRadians, double declinationRadians) {
    const double tanProduct = std::tan(latitudeRadians) * std::tan(declinationRadians);
    double hourAngleLimit = 0.0;

    if (tanProduct >= 1.0) {
        // Полярный день: Солнце не заходит, H0 = π.
        hourAngleLimit = kPi;
    } else if (tanProduct <= -1.0) {
        // Полярная ночь: Солнце не восходит, H0 = 0.
        hourAngleLimit = 0.0;
    } else {
        hourAngleLimit = std::acos(-tanProduct);
    }

    const double meanCosine =
        (hourAngleLimit * std::sin(latitudeRadians) * std::sin(declinationRadians) +
         std::cos(latitudeRadians) * std::cos(declinationRadians) * std::sin(hourAngleLimit)) / kPi;
    return qMax(0.0, meanCosine);
}

double estimateAtmosphereSpecificHeatCp(const AtmosphereComposition &atmosphere) {
    const auto fractions = atmosphere.fractions();
    if (fractions.isEmpty()) {
        return kReferenceCpDry;
    }

    double greenhouseShare = 0.0;
    double totalShare = 0.0;
    const auto gases = availableGases();
    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }
        const auto it = std::find_if(gases.begin(), gases.end(),
                                     [&fraction](const GasSpec &spec) {
                                         return spec.id == fraction.id;
                                     });
        if (it == gases.end()) {
            continue;
        }
        totalShare += fraction.share;
        if (it->isGreenhouse) {
            greenhouseShare += fraction.share;
        }
    }

    if (totalShare <= 0.0) {
        return kReferenceCpDry;
    }

    // Cp оцениваем как смесь сухого воздуха и более "влажных" компонент:
    // парниковые газы повышают теплоёмкость, поэтому смещаем Cp к верхней границе.
    const double wetShare = qBound(0.0, greenhouseShare / totalShare, 1.0);
    return kReferenceCpDry + wetShare * (kReferenceCpWet - kReferenceCpDry);
}

double atmosphereHeatCapacityPerArea(double pressureAtm,
                                     double surfaceGravity,
                                     const AtmosphereComposition &atmosphere) {
    if (pressureAtm <= 0.0 || surfaceGravity <= 0.0) {
        return 0.0;
    }

    // Используем столб атмосферы как добавочную теплоёмкость поверхности:
    // m_atm / A = P / g, поэтому C_atm = (P / g) * c_p.
    // Предполагаем, что атмосфера теплообменно связана с поверхностью
    // и эффективно участвует в суточном балансе.
    const double massPerArea = pressureAtm * kPascalPerAtm / surfaceGravity;
    const double cp = estimateAtmosphereSpecificHeatCp(atmosphere);
    return qMax(0.0, massPerArea * cp);
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

    const bool hasAtmosphere = useAtmosphericModel_;
    points.reserve(latitudePoints);

    int processedLatitudes = 0;
    const double declinationRadians = qDegreesToRadians(declinationDegrees);
    const double stepDegrees = 180.0 / static_cast<double>(latitudePoints - 1);

    for (int i = 0; i < latitudePoints; ++i) {
        if (cancelFlag && cancelFlag->load()) {
            return {};
        }
        const double axisDegrees = (rotationMode_ == RotationMode::Normal)
                                       ? (-90.0 + stepDegrees * static_cast<double>(i))
                                       : (stepDegrees * static_cast<double>(i));
        // Ось X: широта (-90..90) для обычного вращения или
        // угловое расстояние от подсолнечной точки (0..180) для приливной синхронизации.
        const double latitudeRadians = qDegreesToRadians(axisDegrees);
        // Для tidal-locked режима геометрия задается через угол от подсолнечной точки:
        // это и есть зенитный угол падения лучей, поэтому не применяем деклинацию.
        const double subsolarAngleRadians = latitudeRadians;
        const double averageCosine =
            (rotationMode_ == RotationMode::Normal)
                ? meanDailyCosine(latitudeRadians, declinationRadians)
                // В приливной синхронизации угол от подсолнечной точки равен зенитному углу,
                // поэтому cos(угол) задает локальную среднюю инсоляцию без суточного вращения.
                : qMax(0.0, std::cos(subsolarAngleRadians));
        // Признак освещенности: используем средний косинус, который уже учитывает полярную ночь.
        const bool hasInsolation = averageCosine > 0.0;
        const double dayLengthSeconds = qMax(0.01, dayLengthDays_) * kSecondsPerEarthDay;
        const double albedo = qBound(0.0, material_.albedo, 1.0);
        const double surfaceHeatCapacity = qMax(1.0, material_.heatCapacity);
        const double atmosphereHeatCapacity =
            hasAtmosphere
                ? atmosphereHeatCapacityPerArea(atmospherePressureAtm_, surfaceGravity_,
                                                atmosphere_)
                : 0.0;
        // Полная теплоёмкость слоя: поверхность + эффективно вовлечённая атмосфера.
        // C_total = C_surface + C_atm, где C_atm = (P / g) * c_p.
        const double heatCapacity = qMax(1.0, surfaceHeatCapacity + atmosphereHeatCapacity);
        // Если атмосфера отсутствует, используем чистую поверхность без радиации/циркуляции.
        const double meanSolarFlux = segmentSolarConstant * averageCosine;
        double initialTemperature = kSpaceTemperatureKelvin;
        // Атмосферный парниковый слой моделируется через эффективную оптическую толщину:
        // входящий поток ослабляется exp(-tau_sw), исходящий — exp(-tau_lw).
        std::optional<AtmosphericRadiationModel> radiationModel;
        double outgoingTransmission = 1.0;
        double greenhouseOpacity = 0.0;
        double adjustedMeanFlux = meanSolarFlux;
        double absorbedMeanFlux = qMax(0.0, adjustedMeanFlux * (1.0 - albedo));
        auto updateRadiativeTemperature = [&](double absorbedFlux, double outgoing) {
            const double greenhouseFactor = qMax(1e-6, outgoing);
            if (absorbedFlux > 0.0) {
                return std::pow(absorbedFlux / (kStefanBoltzmannConstant * greenhouseFactor), 0.25);
            }
            return kSpaceTemperatureKelvin;
        };

        // Сначала оцениваем радиационное равновесие без температурно-зависимой оптики.
        initialTemperature = updateRadiativeTemperature(absorbedMeanFlux, outgoingTransmission);

        if (hasAtmosphere) {
            // Делаем 1–2 итерации: температура -> оптика -> обновлённая температура.
            for (int iteration = 0; iteration < 2; ++iteration) {
                radiationModel.emplace(atmosphere_, atmospherePressureAtm_, initialTemperature);
                // Преобразуем эффективную оптическую толщину в прозрачность для излучения:
                // применяем ту же формулу, что и при расчёте исходящего потока.
                outgoingTransmission = radiationModel->applyOutgoingFlux(1.0);
                greenhouseOpacity = 1.0 - outgoingTransmission;
                adjustedMeanFlux = radiationModel->applyIncomingFlux(meanSolarFlux);
                absorbedMeanFlux = qMax(0.0, adjustedMeanFlux * (1.0 - albedo));
                initialTemperature =
                    updateRadiativeTemperature(absorbedMeanFlux, outgoingTransmission);
            }
        } else {
            greenhouseOpacity = qBound(0.0, greenhouseOpacity_, 0.999);
            outgoingTransmission = 1.0 - greenhouseOpacity;
            initialTemperature =
                updateRadiativeTemperature(absorbedMeanFlux, outgoingTransmission);
        }
        const double stableTimeStep = kBaseStepSeconds;
        const int stepsPerDay = qBound(
            kMinStepsPerDay,
            static_cast<int>(std::ceil(dayLengthSeconds / stableTimeStep)),
            kMaxStepsPerDay);
        const double timeStepSeconds = dayLengthSeconds / stepsPerDay;

        SurfaceTemperatureState temperatureState(initialTemperature, albedo, heatCapacity,
                                                  greenhouseOpacity);
        double minimumTemperature = std::numeric_limits<double>::max();
        double maximumTemperature = 0.0;
        double meanDailySum = 0.0;
        double meanDaySum = 0.0;
        double meanNightSum = 0.0;
        int meanDailyCount = 0;
        int meanDayCount = 0;
        int meanNightCount = 0;

        constexpr int kSpinupCycles = 3;
        const double tidalSolarFactor = qMax(0.0, std::cos(subsolarAngleRadians));
        for (int cycle = 0; cycle < kSpinupCycles; ++cycle) {
            for (int step = 0; step < stepsPerDay; ++step) {
                if (cancelFlag && cancelFlag->load()) {
                    return {};
                }
                double solarFactor = tidalSolarFactor;
                if (rotationMode_ == RotationMode::Normal) {
                    const double phase = static_cast<double>(step) / stepsPerDay;
                    const double hourAngle = 2.0 * kPi * phase - kPi;
                    // Солнечный фактор: cos(зенитного угла) с поправкой на сезонную деклинацию,
                    // где широта задаёт наклон поверхности относительно лучей звезды.
                    solarFactor = std::sin(latitudeRadians) * std::sin(declinationRadians) +
                                  std::cos(latitudeRadians) * std::cos(declinationRadians) *
                                      std::cos(hourAngle);
                }
                double solarFlux = segmentSolarConstant * qMax(0.0, solarFactor);
                if (hasAtmosphere) {
                    solarFlux = radiationModel->applyIncomingFlux(solarFlux);
                }

                temperatureState.updateTemperature(solarFlux, timeStepSeconds);

                if (cycle == kSpinupCycles - 1) {
                    const double surfaceTemperature = temperatureState.temperatureKelvin();
                    minimumTemperature = qMin(minimumTemperature, surfaceTemperature);
                    maximumTemperature = qMax(maximumTemperature, surfaceTemperature);
                    meanDailySum += surfaceTemperature;
                    ++meanDailyCount;
                    // Усредняем по шагам, разделяя день/ночь по знаку solarFactor:
                    // так видны отдельные характеристики нагрева и остывания,
                    // особенно при длительном полярном дне или ночи.
                    if (solarFactor > 0.0) {
                        meanDaySum += surfaceTemperature;
                        ++meanDayCount;
                    } else {
                        meanNightSum += surfaceTemperature;
                        ++meanNightCount;
                    }
                }
            }
        }

        const double meanDayTemperature =
            (meanDayCount > 0) ? (meanDaySum / meanDayCount)
                               : ((meanNightCount > 0) ? (meanNightSum / meanNightCount)
                                                       : initialTemperature);
        const double meanNightTemperature =
            (meanNightCount > 0) ? (meanNightSum / meanNightCount)
                                 : ((meanDayCount > 0) ? (meanDaySum / meanDayCount)
                                                       : initialTemperature);
        const double meanDailyTemperature =
            (meanDailyCount > 0) ? (meanDailySum / meanDailyCount)
                                 : initialTemperature;

        TemperatureRangePoint point;
        point.latitudeDegrees = axisDegrees;
        point.hasInsolation = hasInsolation;
        point.minimumKelvin = minimumTemperature;
        point.maximumKelvin = maximumTemperature;
        point.meanDailyKelvin = meanDailyTemperature;
        point.meanDayKelvin = meanDayTemperature;
        point.meanNightKelvin = meanNightTemperature;
        point.minimumCelsius = minimumTemperature - kKelvinOffset;
        point.maximumCelsius = maximumTemperature - kKelvinOffset;
        point.meanDailyCelsius = meanDailyTemperature - kKelvinOffset;
        point.meanDayCelsius = meanDayTemperature - kKelvinOffset;
        point.meanNightCelsius = meanNightTemperature - kKelvinOffset;
        points.push_back(point);

        ++processedLatitudes;
        if (progressCallback) {
            progressCallback(progressOffset + processedLatitudes, totalProgress);
        }
    }

    return points;
}
