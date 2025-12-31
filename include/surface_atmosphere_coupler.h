#pragma once

#include "atmospheric_cell_state.h"

class SurfaceTemperatureState;

class SurfaceAtmosphereCoupler {
public:
    SurfaceAtmosphereCoupler() = default;
    explicit SurfaceAtmosphereCoupler(double heatTransferCoefficientWPerM2K);

    void setHeatTransferCoefficientWPerM2K(double coefficient);
    double heatTransferCoefficientWPerM2K() const;

    void exchangeSensibleHeat(SurfaceTemperatureState &surface,
                              AtmosphericCellState &atmosphere,
                              double dtSeconds) const;

private:
    double heatTransferCoefficientWPerM2K_ = 8.0;
};
