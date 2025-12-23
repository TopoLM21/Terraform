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
constexpr double kDensity = 1500.0;
constexpr double kSpecificHeat = 800.0;
constexpr double kThermalConductivity = 0.5;
constexpr double kMaxSurfaceTemperature = 4000.0;
}  // namespace

SurfaceTemperatureCalculator::SurfaceTemperatureCalculator(double solarConstant,
                                                           const SurfaceMaterial &material,
                                                           double dayLengthDays)
    : solarConstant_(solarConstant), material_(material), dayLengthDays_(dayLengthDays) {}

QVector<TemperatureRangePoint> SurfaceTemperatureCalculator::temperatureRangesByLatitude(
    int stepDegrees) const {
    QVector<TemperatureRangePoint> points;
    if (stepDegrees <= 0) {
        return points;
    }

    const int steps = 180 / stepDegrees;
    points.reserve(steps + 1);

    for (int latitude = -90; latitude <= 90; latitude += stepDegrees) {
        const double latitudeRadians = qDegreesToRadians(static_cast<double>(latitude));
        const double dayLengthSeconds = qMax(0.01, dayLengthDays_) * kSecondsPerEarthDay;
        const double layerThickness = kSurfaceDepthMeters / kLayerCount;
        const double heatCapacity = kDensity * kSpecificHeat * layerThickness;
        const double conductionFactor = kThermalConductivity / layerThickness;
        const double emissivity = qMax(0.0001, material_.emissivity);
        const double spaceTemperaturePower = std::pow(kSpaceTemperatureKelvin, 4.0);
        const double maxSolarFlux = solarConstant_ * (1.0 - material_.albedo);
        const double maxRadiativeTemperature =
            std::pow((maxSolarFlux + kInternalHeatFlux) /
                         (emissivity * kStefanBoltzmannConstant) +
                     spaceTemperaturePower,
                     0.25);
        // Ограничиваем шаг по времени по условиям устойчивости явной схемы
        // для теплопроводности и радиационного охлаждения.
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
        const double averageCosine = qMax(0.0, std::cos(latitudeRadians)) / kPi;
        // Усредненная за сутки инсоляция нужна для старта итераций без резких скачков.
        const double meanSolarFlux = solarConstant_ * averageCosine * (1.0 - material_.albedo);
        const double meanFlux = meanSolarFlux + kInternalHeatFlux +
                                emissivity * kStefanBoltzmannConstant * spaceTemperaturePower;
        const double initialTemperature = std::pow(meanFlux / (emissivity * kStefanBoltzmannConstant),
                                                   0.25);

        QVector<double> layers(kLayerCount, initialTemperature);
        double minimumTemperature = std::numeric_limits<double>::max();
        double maximumTemperature = 0.0;

        constexpr int kSpinupCycles = 3;
        for (int cycle = 0; cycle < kSpinupCycles; ++cycle) {
            for (int step = 0; step < stepsPerDay; ++step) {
                const double phase = static_cast<double>(step) / stepsPerDay;
                const double hourAngle = 2.0 * kPi * phase - kPi;
                const double solarFactor = std::cos(latitudeRadians) * std::cos(hourAngle);
                const double solarFlux = solarConstant_ * qMax(0.0, solarFactor) *
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
                    minimumTemperature = qMin(minimumTemperature, layers[0]);
                    maximumTemperature = qMax(maximumTemperature, layers[0]);
                }
            }
        }

        TemperatureRangePoint point;
        point.latitudeDegrees = latitude;
        point.minimumKelvin = minimumTemperature;
        point.maximumKelvin = maximumTemperature;
        point.minimumCelsius = minimumTemperature - kKelvinOffset;
        point.maximumCelsius = maximumTemperature - kKelvinOffset;
        points.push_back(point);
    }

    return points;
}
