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

    struct LongwaveFluxProfile {
        // Потоки на границах слоёв: индекс 0 соответствует поверхности, последний — верх атмосферы.
        QVector<double> upwardFluxWPerM2;
        QVector<double> downwardFluxWPerM2;
    };

    LayeredRadiationModel(const AtmosphereComposition &composition,
                          double pressureAtm,
                          double baseTemperatureKelvin,
                          double effectiveTemperatureKelvin,
                          double surfaceGravity);

    double effectiveOpticalDepth() const override;
    double incomingTransmission() const override;
    double outgoingTransmission() const override;
    double bottomLayerTemperatureKelvin() const;
    LongwaveFluxProfile longwaveFluxProfile(double surfaceEmittedFluxWPerM2) const;

private:
    void buildLayers();
    double layerLongwaveEmissivity(double opticalDepth) const;

    AtmosphereComposition composition_;
    double pressureAtm_ = 0.0;
    double baseTemperatureKelvin_ = 0.0;
    double effectiveTemperatureKelvin_ = 0.0;
    double surfaceGravity_ = 0.0;
    double effectiveOpticalDepth_ = 0.0;
    double shortwaveOpticalDepth_ = 0.0;
    QVector<Layer> layers_;
};
