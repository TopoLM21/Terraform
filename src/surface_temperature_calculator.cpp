#include "surface_temperature_calculator.h"

#include <QtCore/QtMath>

#include <cmath>
#include <limits>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kKelvinOffset = 273.15;
constexpr double kPi = 3.14159265358979323846;
constexpr double kSpaceTemperatureKelvin = 3.0;
constexpr double kInternalHeatFlux = 0.05;
constexpr double kSecondsPerEarthDay = 86400.0;
constexpr double kBaseStepSeconds = 900.0;
constexpr int kMinStepsPerDay = 48;
constexpr int kMaxStepsPerDay = 20000;
constexpr int kLayerCount = 8;
constexpr double kSurfaceDepthMeters = 1.0;
constexpr double kMaxSurfaceTemperature = 4000.0;

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
}  // namespace

SurfaceTemperatureCalculator::SurfaceTemperatureCalculator(double solarConstant,
                                                           const SurfaceMaterial &material,
                                                           double dayLengthDays)
    : solarConstant_(solarConstant), material_(material), dayLengthDays_(dayLengthDays) {}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesByLatitude(
    int stepDegrees) const {
    return temperatureRangesByLatitude(stepDegrees, ProgressCallback{}, nullptr);
}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesByLatitude(
    int stepDegrees,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag) const {
    const int totalLatitudes = stepDegrees > 0 ? 180 / stepDegrees + 1 : 0;
    return temperatureRangesByLatitudeForSegment(stepDegrees, solarConstant_, 0.0, progressCallback,
                                                 cancelFlag, 0, totalLatitudes);
}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesForOrbitSegment(
    const OrbitSegment &segment,
    double referenceDistanceAU,
    double obliquityDegrees,
    double perihelionArgumentDegrees,
    int stepDegrees,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag) const {
    if (cancelFlag && cancelFlag->load()) {
        return {};
    }
    if (stepDegrees <= 0) {
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
    const int totalLatitudes = 180 / stepDegrees + 1;

    return temperatureRangesByLatitudeForSegment(stepDegrees,
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
    int stepDegrees,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag) const {
    QVector<QVector<TemperatureRangePoint>> results;
    if (cancelFlag && cancelFlag->load()) {
        return results;
    }
    if (segments.isEmpty() || stepDegrees <= 0) {
        return results;
    }

    const double obliquityRadians = qDegreesToRadians(obliquityDegrees);
    const double perihelionArgumentRadians = qDegreesToRadians(perihelionArgumentDegrees);
    const int latitudesCount = 180 / stepDegrees + 1;
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

        results.push_back(temperatureRangesByLatitudeForSegment(stepDegrees,
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
    int stepDegrees,
    double segmentSolarConstant,
    double declinationDegrees,
    const ProgressCallback &progressCallback,
    const std::atomic_bool *cancelFlag,
    int progressOffset,
    int totalProgress) const {
    QVector<TemperatureRangePoint> points;
    if (stepDegrees <= 0) {
        return points;
    }
    if (cancelFlag && cancelFlag->load()) {
        return {};
    }

    const int steps = 180 / stepDegrees;
    points.reserve(steps + 1);

    int processedLatitudes = 0;
    const double declinationRadians = qDegreesToRadians(declinationDegrees);

    for (int latitude = -90; latitude <= 90; latitude += stepDegrees) {
        if (cancelFlag && cancelFlag->load()) {
            return {};
        }
        const double latitudeRadians = qDegreesToRadians(static_cast<double>(latitude));
        const double dayLengthSeconds = qMax(0.01, dayLengthDays_) * kSecondsPerEarthDay;
        const double layerThickness = kSurfaceDepthMeters / kLayerCount;
        const double density = qMax(1.0, material_.density);
        const double specificHeat = qMax(1.0, material_.specificHeat);
        const double thermalConductivity = qMax(1e-6, material_.thermalConductivity);
        // Тепловая инерция I = sqrt(k * rho * c) описывает устойчивость поверхности к нагреву.
        const double thermalInertia =
            std::sqrt(thermalConductivity * density * specificHeat);
        const double heatCapacity = density * specificHeat * layerThickness;
        // Коэффициент теплопереноса между слоями: k / dz.
        const double conductionFactor = thermalConductivity / layerThickness;
        const double emissivity = qMax(0.0001, material_.emissivity);
        const double spaceTemperaturePower = std::pow(kSpaceTemperatureKelvin, 4.0);
        const double maxSolarFlux = segmentSolarConstant * (1.0 - material_.albedo);
        const double maxRadiativeTemperature =
            std::pow((maxSolarFlux + kInternalHeatFlux) /
                         (emissivity * kStefanBoltzmannConstant) +
                     spaceTemperaturePower,
                     0.25);
        // Ограничиваем шаг по времени по условиям устойчивости явной схемы:
        // Δt_cond ~ (rho * c * dz^2) / (2 * k), Δt_rad ~ (rho * c * dz) / (4 * εσT^3).
        // thermalInertia влияет на выбор Δt_cond через k, rho и c.
        (void)thermalInertia;
        const double conductionTimeStep = heatCapacity / (2.0 * conductionFactor);
        const double radiativeTimeStep =
            heatCapacity /
            (4.0 * emissivity * kStefanBoltzmannConstant *
             std::pow(qMax(1.0, maxRadiativeTemperature), 3.0));
        const double stableTimeStep = qMin(kBaseStepSeconds, qMin(conductionTimeStep, radiativeTimeStep));
        const int stepsPerDay = qBound(
            kMinStepsPerDay,
            static_cast<int>(std::ceil(dayLengthSeconds / stableTimeStep)),
            kMaxStepsPerDay);
        const double timeStepSeconds = dayLengthSeconds / stepsPerDay;
        // Усредненная за сутки инсоляция: учитываем сезонную деклинацию.
        const double averageCosine = meanDailyCosine(latitudeRadians, declinationRadians);
        const double meanSolarFlux =
            segmentSolarConstant * averageCosine * (1.0 - material_.albedo);
        const double meanFlux = meanSolarFlux + kInternalHeatFlux +
                                emissivity * kStefanBoltzmannConstant * spaceTemperaturePower;
        const double initialTemperature = std::pow(meanFlux / (emissivity * kStefanBoltzmannConstant),
                                                   0.25);

        QVector<double> layers(kLayerCount, initialTemperature);
        double minimumTemperature = std::numeric_limits<double>::max();
        double maximumTemperature = 0.0;
        double meanDailySum = 0.0;
        double meanDaySum = 0.0;
        double meanNightSum = 0.0;
        int meanDailyCount = 0;
        int meanDayCount = 0;
        int meanNightCount = 0;

        constexpr int kSpinupCycles = 3;
        for (int cycle = 0; cycle < kSpinupCycles; ++cycle) {
            for (int step = 0; step < stepsPerDay; ++step) {
                if (cancelFlag && cancelFlag->load()) {
                    return {};
                }
                const double phase = static_cast<double>(step) / stepsPerDay;
                const double hourAngle = 2.0 * kPi * phase - kPi;
                // Солнечный фактор: cos(зенитного угла) с поправкой на сезонную деклинацию.
                const double solarFactor = std::sin(latitudeRadians) * std::sin(declinationRadians) +
                                           std::cos(latitudeRadians) * std::cos(declinationRadians) *
                                               std::cos(hourAngle);
                const double solarFlux = segmentSolarConstant * qMax(0.0, solarFactor) *
                                         (1.0 - material_.albedo);

                QVector<double> nextLayers = layers;
                for (int layer = 0; layer < kLayerCount; ++layer) {
                    double conductionFlux = 0.0;
                    if (layer > 0) {
                        conductionFlux += conductionFactor * (layers[layer - 1] - layers[layer]);
                    }
                    if (layer + 1 < kLayerCount) {
                        conductionFlux += conductionFactor * (layers[layer + 1] - layers[layer]);
                    }

                    double netFlux = conductionFlux;
                    if (layer == 0) {
                        const double surfaceTemperature = qMax(1.0, layers[layer]);
                        const double surfacePower = std::pow(surfaceTemperature, 4.0);
                        const double radiationLoss =
                            emissivity * kStefanBoltzmannConstant *
                            (surfacePower - spaceTemperaturePower);
                        // Радиационный баланс: излучение в холодный космос не дает остыть до 0 K.
                        netFlux += solarFlux - radiationLoss;
                    }

                    if (layer == kLayerCount - 1) {
                        // Внутренний тепловой поток подается снизу, чтобы отрабатывать по слоям.
                        netFlux += kInternalHeatFlux;
                    }

                    const double nextTemperature =
                        layers[layer] + (netFlux * timeStepSeconds) / heatCapacity;
                    // Ограничиваем диапазон, чтобы при очень длинных сутках не было числового переполнения.
                    nextLayers[layer] = qBound(kSpaceTemperatureKelvin, nextTemperature,
                                               kMaxSurfaceTemperature);
                }

                layers = nextLayers;

                if (cycle == kSpinupCycles - 1) {
                    const double surfaceTemperature = layers[0];
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
        point.latitudeDegrees = latitude;
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
