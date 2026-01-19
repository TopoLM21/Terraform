#pragma once

#include "atmosphere_model.h"
#include "atmospheric_column.h"

class ConvectiveAdjustmentSolver {
public:
    ConvectiveAdjustmentSolver(const AtmosphereComposition &composition, double gravityMps2);

    void adjust(AtmosphericColumn &column) const;

private:
    double adiabaticGradientDlnT(double temperatureKelvin, double pressureAtm) const;

    AtmosphereComposition composition_;
    bool hasWaterVapor_ = false;
    double gravityMps2_ = 0.0;
    double specificHeatCp_ = 0.0;
    double specificGasConstant_ = 0.0;
};
