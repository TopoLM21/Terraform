#pragma once

#include "atmospheric_cell_state.h"
#include "atmospheric_layer_state.h"

class SurfacePointState;

class SurfaceAtmosphereCoupler {
public:
    SurfaceAtmosphereCoupler() = default;
    explicit SurfaceAtmosphereCoupler(double heatTransferCoefficientWPerM2K);

    void setHeatTransferCoefficientWPerM2K(double coefficient);
    double heatTransferCoefficientWPerM2K() const;

    void exchangeHeat(SurfacePointState &surface,
                      AtmosphericLayerState &lowestLayer,
                      double roughnessLengthMeters,
                      double dtSeconds) const;

    void exchangeHeat(SurfacePointState &surface,
                      AtmosphericCellState &atmosphere,
                      double windSpeedMps,
                      double longwaveEmissivity,
                      double roughnessLengthMeters,
                      double dtSeconds) const;

private:
    double sensibleHeatTransferCoefficient(double windSpeedMps,
                                           double roughnessLengthMeters,
                                           double airHeatCapacityJPerM2K,
                                           double layerThicknessMeters) const;
    double longwaveEmissivityFromOpticalDepth(double opticalDepth) const;
    void exchangeHeatInternal(SurfacePointState &surface,
                              double airTemperatureKelvin,
                              double airHeatCapacityJPerM2K,
                              double windSpeedMps,
                              double longwaveEmissivity,
                              double roughnessLengthMeters,
                              double layerThicknessMeters,
                              double dtSeconds,
                              double &updatedAirTemperature) const;

    double heatTransferCoefficientWPerM2K_ = 8.0;
};
