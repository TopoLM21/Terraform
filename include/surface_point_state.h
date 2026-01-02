#pragma once

#include "subsurface_temperature_solver.h"

class SurfaceMaterial;

class SurfacePointState {
public:
    SurfacePointState() = default;
    SurfacePointState(double initialTemperatureKelvin,
                      double albedo,
                      double greenhouseOpacity,
                      double minTemperatureKelvin,
                      const SurfaceMaterial &material,
                      const SubsurfaceModelSettings &subsurfaceSettings);

    double temperatureKelvin() const;
    void setTemperatureKelvin(double temperatureKelvin);
    void setGreenhouseOpacity(double greenhouseOpacity);
    const SubsurfaceTemperatureSolver &solver() const;

    double absorbedFlux(double solarIrradiance) const;
    double emittedFlux() const;

    void updateTemperature(double absorbedFlux, double emittedFlux, double dtSeconds);

private:
    void clampProfile();

    double albedo_ = 0.0;
    double greenhouseOpacity_ = 0.0;
    double minTemperatureKelvin_ = 3.0;
    SubsurfaceTemperatureSolver solver_;
};
