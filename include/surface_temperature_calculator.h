#pragma once

#include "orbit_segment_calculator.h"
#include "planet_presets.h"

#include <QtCore/QVector>

#include <atomic>
#include <functional>

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

    using ProgressCallback = std::function<void(int processed, int total)>;

    QVector<TemperatureRangePoint> temperatureRangesByLatitude(int stepDegrees = 1) const;
    QVector<TemperatureRangePoint> temperatureRangesByLatitude(int stepDegrees,
                                                               const ProgressCallback &progressCallback,
                                                               const std::atomic_bool *cancelFlag) const;
    QVector<QVector<TemperatureRangePoint>> temperatureRangesByOrbitSegments(
        const QVector<OrbitSegment> &segments,
        double referenceDistanceAU,
        double obliquityDegrees,
        double perihelionArgumentDegrees,
        int stepDegrees,
        const ProgressCallback &progressCallback,
        const std::atomic_bool *cancelFlag) const;

private:
    QVector<TemperatureRangePoint> temperatureRangesByLatitudeForSegment(
        int stepDegrees,
        double segmentSolarConstant,
        double declinationDegrees,
        const ProgressCallback &progressCallback,
        const std::atomic_bool *cancelFlag,
        int progressOffset,
        int totalProgress) const;

    double solarConstant_;
    SurfaceMaterial material_;
    double dayLengthDays_;
};
