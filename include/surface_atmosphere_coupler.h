#pragma once

#include "atmospheric_cell_state.h"

class SurfacePointState;

class SurfaceAtmosphereCoupler {
public:
    SurfaceAtmosphereCoupler() = default;
    explicit SurfaceAtmosphereCoupler(double heatTransferCoefficientWPerM2K);

    void setHeatTransferCoefficientWPerM2K(double coefficient);
    double heatTransferCoefficientWPerM2K() const;

    void exchangeSensibleHeat(SurfacePointState &surface,
                              AtmosphericCellState &atmosphere,
                              double dtSeconds) const;

private:
    double heatTransferCoefficientWPerM2K_ = 8.0;
};
