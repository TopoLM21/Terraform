#pragma once

#include "atmosphere_model.h"
#include "orbit_segment_calculator.h"
#include "planet_presets.h"
#include "rotation_mode.h"

#include <QtCore/QVector>

#include <atomic>
#include <functional>

struct TemperatureRangePoint {
    double latitudeDegrees;
    bool hasInsolation = false;
    double minimumKelvin;
    double maximumKelvin;
    double meanDailyKelvin;
    double meanDayKelvin;
    double meanNightKelvin;
    double minimumCelsius;
    double maximumCelsius;
    double meanDailyCelsius;
    double meanDayCelsius;
    double meanNightCelsius;
};

struct TemperatureSummaryPoint {
    double latitudeDegrees;
    double minimumKelvin;
    double maximumKelvin;
    double meanAnnualKelvin;
    double meanAnnualDayKelvin;
    double meanAnnualNightKelvin;
    double minimumCelsius;
    double maximumCelsius;
    double meanAnnualCelsius;
    double meanAnnualDayCelsius;
    double meanAnnualNightCelsius;
};

class SurfaceTemperatureCalculator {
public:
    SurfaceTemperatureCalculator(double solarConstant, const SurfaceMaterial &material,
                                 double dayLengthDays, RotationMode rotationMode,
                                 const AtmosphereComposition &atmosphere = AtmosphereComposition{},
                                 double atmospherePressureAtm = 0.0,
                                 double surfaceGravity = 0.0,
                                 bool useAtmosphericModel = false,
                                 int meridionalTransportSteps = 8);

    void setMeridionalTransportSteps(int steps);

    using ProgressCallback = std::function<void(int processed, int total)>;

    QVector<TemperatureRangePoint> temperatureRangesByLatitude(int latitudePoints = 181) const;
    QVector<TemperatureRangePoint> temperatureRangesByLatitude(int latitudePoints,
                                                               const ProgressCallback &progressCallback,
                                                               const std::atomic_bool *cancelFlag) const;
    QVector<TemperatureRangePoint> temperatureRangesForOrbitSegment(
        const OrbitSegment &segment,
        double referenceDistanceAU,
        double obliquityDegrees,
        double perihelionArgumentDegrees,
        int latitudePoints,
        const ProgressCallback &progressCallback,
        const std::atomic_bool *cancelFlag) const;
    QVector<QVector<TemperatureRangePoint>> temperatureRangesByOrbitSegments(
        const QVector<OrbitSegment> &segments,
        double referenceDistanceAU,
        double obliquityDegrees,
        double perihelionArgumentDegrees,
        int latitudePoints,
        const ProgressCallback &progressCallback,
        const std::atomic_bool *cancelFlag) const;

private:
    QVector<TemperatureRangePoint> temperatureRangesByLatitudeForSegment(
        int latitudePoints,
        double segmentSolarConstant,
        double declinationDegrees,
        const ProgressCallback &progressCallback,
        const std::atomic_bool *cancelFlag,
        int progressOffset,
        int totalProgress) const;

    QVector<TemperatureRangePoint> radiativeBalanceByLatitudeForSegment(
        int latitudePoints,
        double segmentSolarConstant,
        double declinationDegrees,
        const ProgressCallback &progressCallback,
        const std::atomic_bool *cancelFlag,
        int progressOffset,
        int totalProgress) const;

    double solarConstant_;
    SurfaceMaterial material_;
    double dayLengthDays_;
    RotationMode rotationMode_ = RotationMode::Normal;
    AtmosphereComposition atmosphere_;
    double atmospherePressureAtm_ = 0.0;
    double surfaceGravity_ = 0.0;
    bool useAtmosphericModel_ = false;
    int meridionalTransportSteps_ = 8;
};
