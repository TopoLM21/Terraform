#include "surface_temperature_calculator.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kKelvinOffset = 273.15;
}  // namespace

SurfaceTemperatureCalculator::SurfaceTemperatureCalculator(double solarConstant,
                                                           const SurfaceMaterial &material)
    : solarConstant_(solarConstant), material_(material) {}

QVector<TemperaturePoint> SurfaceTemperatureCalculator::temperaturesByLatitude(int stepDegrees) const {
    QVector<TemperaturePoint> points;
    if (stepDegrees <= 0) {
        return points;
    }

    const int steps = 180 / stepDegrees;
    points.reserve(steps + 1);

    for (int latitude = -90; latitude <= 90; latitude += stepDegrees) {
        const double latitudeRadians = qDegreesToRadians(static_cast<double>(latitude));
        const double cosineFactor = qMax(0.0, std::cos(latitudeRadians));
        const double absorbedFlux = solarConstant_ * cosineFactor * (1.0 - material_.albedo);
        const double effectiveFlux = absorbedFlux / qMax(0.0001, material_.emissivity);
        const double temperatureKelvin = std::pow(effectiveFlux / kStefanBoltzmannConstant, 0.25);

        TemperaturePoint point;
        point.latitudeDegrees = latitude;
        point.temperatureKelvin = temperatureKelvin;
        point.temperatureCelsius = temperatureKelvin - kKelvinOffset;
        points.push_back(point);
    }

    return points;
}
