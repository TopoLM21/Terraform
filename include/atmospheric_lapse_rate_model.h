#pragma once

#include "atmosphere_model.h"

struct TemperatureRangePoint;

class AtmosphericLapseRateModel {
public:
    AtmosphericLapseRateModel(double atmospherePressureAtm,
                              const AtmosphereComposition &atmosphere,
                              double surfaceGravity);

    TemperatureRangePoint applyLapseRate(const TemperatureRangePoint &point) const;

private:
    double atmospherePressureAtm_ = 0.0;
    AtmosphereComposition atmosphere_;
    double surfaceGravity_ = 0.0;

    double meanMolarMassKgPerMol() const;
    double specificGasConstant() const;
    double specificHeatCp() const;
    double relativeHumidityEstimate() const;
    double saturationVaporPressureAtm(double temperatureKelvin) const;
    double moistAdiabaticLapseRate(double temperatureKelvin) const;
    double dryAdiabaticLapseRate() const;
    double scaleHeightMeters(double temperatureKelvin) const;
    double surfaceAdjustmentHeightMeters(double temperatureKelvin) const;
};
