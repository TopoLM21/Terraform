#pragma once

#include "planet_presets.h"

#include <QtCore/QVector>

struct TemperatureRangePoint {
    double latitudeDegrees;
    double minimumKelvin;
    double maximumKelvin;
    double minimumCelsius;
    double maximumCelsius;
};

class SurfaceTemperatureCalculator {
public:
    SurfaceTemperatureCalculator(double solarConstant, const SurfaceMaterial &material,
                                 double dayLengthDays);

    QVector<TemperatureRangePoint> temperatureRangesByLatitude(int stepDegrees = 15) const;

private:
    double solarConstant_;
    SurfaceMaterial material_;
    double dayLengthDays_;
};
