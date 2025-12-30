#include "surface_temperature_state.h"

#include "planet_presets.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
}

SurfaceTemperatureState::SurfaceTemperatureState(double initialTemperatureKelvin,
                                                 double albedo,
                                                 double greenhouseOpacity,
                                                 double minTemperatureKelvin,
                                                 const SurfaceMaterial &material,
                                                 const SubsurfaceModelSettings &subsurfaceSettings)
    : albedo_(qBound(0.0, albedo, 1.0)),
      greenhouseOpacity_(qBound(0.0, greenhouseOpacity, 0.999)),
      minTemperatureKelvin_(qMax(0.0, minTemperatureKelvin)) {
    const SubsurfaceGrid grid(subsurfaceSettings.layerCount,
                              subsurfaceSettings.topLayerThicknessMeters,
                              subsurfaceSettings.bottomDepthMeters);
    solver_.reset(grid,
                  material.thermalConductivity,
                  material.density,
                  material.specificHeat,
                  subsurfaceSettings.bottomBoundary,
                  initialTemperatureKelvin);
    solver_.setInitialTemperature(qMax(minTemperatureKelvin_, initialTemperatureKelvin));
}

double SurfaceTemperatureState::temperatureKelvin() const {
    return solver_.surfaceTemperatureKelvin();
}

void SurfaceTemperatureState::setTemperatureKelvin(double temperatureKelvin) {
    solver_.setInitialTemperature(qMax(minTemperatureKelvin_, temperatureKelvin));
}

double SurfaceTemperatureState::absorbedFlux(double solarIrradiance) const {
    return solarIrradiance * (1.0 - albedo_);
}

double SurfaceTemperatureState::emittedFlux() const {
    return kStefanBoltzmannConstant * std::pow(temperatureKelvin(), 4.0) *
           (1.0 - greenhouseOpacity_);
}

void SurfaceTemperatureState::updateTemperature(double absorbedFlux,
                                                double emittedFlux,
                                                double dtSeconds) {
    // Радиативный баланс для верхнего слоя:
    // E_in = absorbedFlux, E_out = emittedFlux,
    // поток на границе = E_in - E_out.
    const double netFlux = absorbedFlux - emittedFlux;
    solver_.stepImplicit(netFlux, dtSeconds);
    clampProfile();
}

void SurfaceTemperatureState::clampProfile() {
    auto profile = solver_.temperatures();
    for (double &value : profile) {
        value = qMax(minTemperatureKelvin_, value);
    }
    solver_.setTemperatures(profile);
}
