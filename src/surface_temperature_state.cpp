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

double SurfaceTemperatureState::topLayerHeatCapacityJPerM2K() const {
    return solver_.topLayerHeatCapacity();
}

void SurfaceTemperatureState::updateTemperature(double absorbedFlux,
                                                double emittedFlux,
                                                double dtSeconds) {
    Q_UNUSED(emittedFlux)
    const double surfaceTemperature = temperatureKelvin();
    const double emittedNow = kStefanBoltzmannConstant * std::pow(surfaceTemperature, 4.0) *
        (1.0 - greenhouseOpacity_);
    const double heatCapacity = solver_.topLayerHeatCapacity();
    // Линеаризованный радиационный поток:
    // F_rad(T) ≈ F_rad(T_n) + (dF/dT)|_{T_n} (T - T_n),
    // dF/dT = 4 * sigma * (1 - tau) * T_n^3.
    // Тогда неявная стабилизация даёт
    // F_eff = (F_in - F_rad(T_n)) / (1 + alpha), alpha = (dF/dT) * dt / C.
    // При alpha > 1 подавляется смена знака потока между шагами и исчезают
    // осцилляции температуры из-за слишком сильной радиации.
    const double radiativeDerivative =
        4.0 * kStefanBoltzmannConstant * std::pow(surfaceTemperature, 3.0) *
        (1.0 - greenhouseOpacity_);
    const double alpha = (heatCapacity > 0.0)
        ? (radiativeDerivative * dtSeconds / heatCapacity)
        : 0.0;
    const double netFlux = (absorbedFlux - emittedNow) / (1.0 + qMax(0.0, alpha));
    solver_.stepImplicit(netFlux, dtSeconds);
    clampProfile();
}

void SurfaceTemperatureState::applySurfaceFlux(double netFlux, double dtSeconds) {
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
