#pragma once

#include "radiation_model.h"

#include <QVector>

class LayeredRadiationModel : public RadiationModel {
public:
    struct Layer {
        double pressureAtm = 0.0;
        double temperatureKelvin = 0.0;
        double opticalDepth = 0.0;
        double shortwaveOpticalDepth = 0.0;
    };

    LayeredRadiationModel(const AtmosphereComposition &composition,
                          double pressureAtm,
                          double baseTemperatureKelvin);

    double effectiveOpticalDepth() const override;
    double incomingTransmission() const override;
    double outgoingTransmission() const override;
    double bottomLayerTemperatureKelvin() const;

private:
    void buildLayers();

    AtmosphereComposition composition_;
    double pressureAtm_ = 0.0;
    double baseTemperatureKelvin_ = 0.0;
    double effectiveOpticalDepth_ = 0.0;
    double shortwaveOpticalDepth_ = 0.0;
    QVector<Layer> layers_;
};
