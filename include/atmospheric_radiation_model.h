#pragma once

#include "radiation_model.h"

class AtmosphericRadiationModel : public RadiationModel {
public:
    AtmosphericRadiationModel(const AtmosphereComposition &composition,
                              double pressureAtm,
                              double baseTemperatureKelvin,
                              RadiationModelType radiationModelType = RadiationModelType::Fast);

    double effectiveOpticalDepth() const override;
    double incomingTransmission() const override;
    double outgoingTransmission() const override;
    double applyIncomingFlux(double flux) const;
    double applyOutgoingFlux(double flux) const;

private:
    void computeOpticalDepths();

    AtmosphereComposition composition_;
    double pressureAtm_ = 0.0;
    double baseTemperatureKelvin_ = 0.0;
    RadiationModelType radiationModelType_ = RadiationModelType::Fast;
    double effectiveOpticalDepth_ = 0.0;
    double shortwaveOpticalDepth_ = 0.0;
};
