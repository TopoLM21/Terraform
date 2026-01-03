#pragma once

#include "atmosphere_model.h"

class AtmosphericRadiationModel {
public:
    AtmosphericRadiationModel(const AtmosphereComposition &composition,
                              double pressureAtm,
                              double baseTemperatureKelvin,
                              bool useMultiLayerRadiation = false);

    double effectiveOpticalDepth() const;
    double incomingTransmission() const;
    double outgoingTransmission() const;
    double applyIncomingFlux(double flux) const;
    double applyOutgoingFlux(double flux) const;

private:
    void computeOpticalDepths();

    AtmosphereComposition composition_;
    double pressureAtm_ = 0.0;
    double baseTemperatureKelvin_ = 0.0;
    bool useMultiLayerRadiation_ = false;
    double effectiveOpticalDepth_ = 0.0;
    double shortwaveOpticalDepth_ = 0.0;
};
