#pragma once

#include "planet_presets.h"

#include <QtCore/QVector>

struct TemperaturePoint {
    double latitudeDegrees;
    double temperatureKelvin;
    double temperatureCelsius;
};

class SurfaceTemperatureCalculator {
public:
    SurfaceTemperatureCalculator(double solarConstant, const SurfaceMaterial &material);

    QVector<TemperaturePoint> temperaturesByLatitude(int stepDegrees = 15) const;

private:
    double solarConstant_;
    SurfaceMaterial material_;
};
