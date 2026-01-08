#pragma once

#include "subsurface_temperature_solver.h"

class SurfaceMaterial;

class SurfaceTemperatureState {
public:
    SurfaceTemperatureState(double initialTemperatureKelvin,
                            double albedo,
                            double greenhouseOpacity,
                            double minTemperatureKelvin,
                            const SurfaceMaterial &material,
                            const SubsurfaceModelSettings &subsurfaceSettings);

    double temperatureKelvin() const;
    void setTemperatureKelvin(double temperatureKelvin);

    double absorbedFlux(double solarIrradiance) const;
    double emittedFlux() const;
    double topLayerHeatCapacityJPerM2K() const;

    void updateTemperature(double absorbedFlux, double emittedFlux, double dtSeconds);
    void applySurfaceFlux(double netFlux, double dtSeconds);

private:
    void clampProfile();

    double albedo_ = 0.0;
    double greenhouseOpacity_ = 0.0;
    double minTemperatureKelvin_ = 3.0;
    SubsurfaceTemperatureSolver solver_;
};
