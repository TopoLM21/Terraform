#include "surface_temperature_state.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
}

SurfaceTemperatureState::SurfaceTemperatureState(double initialTemperatureKelvin,
                                                 double albedo,
                                                 double heatCapacity,
                                                 double greenhouseOpacity,
                                                 double minTemperatureKelvin)
    : temperatureKelvin_(qMax(minTemperatureKelvin, initialTemperatureKelvin)),
      albedo_(qBound(0.0, albedo, 1.0)),
      heatCapacity_(qMax(1.0, heatCapacity)),
      greenhouseOpacity_(qBound(0.0, greenhouseOpacity, 0.999)),
      minTemperatureKelvin_(qMax(0.0, minTemperatureKelvin)) {
    temperatureKelvin_ = qMax(minTemperatureKelvin_, temperatureKelvin_);
}

double SurfaceTemperatureState::temperatureKelvin() const {
    return temperatureKelvin_;
}

void SurfaceTemperatureState::setTemperatureKelvin(double temperatureKelvin) {
    temperatureKelvin_ = qMax(minTemperatureKelvin_, temperatureKelvin);
}

void SurfaceTemperatureState::updateTemperature(double solarIrradiance, double dtSeconds) {
    // Радиативный баланс для поверхностного слоя:
    // E_in = S * (1 - A), E_out = sigma * T^4 * (1 - G),
    // dT/dt = (E_in - E_out) / C.
    const double absorbedFlux = solarIrradiance * (1.0 - albedo_);
    const double emittedFlux = kStefanBoltzmannConstant * std::pow(temperatureKelvin_, 4.0) *
                               (1.0 - greenhouseOpacity_);
    const double dT = (absorbedFlux - emittedFlux) / heatCapacity_ * dtSeconds;

    temperatureKelvin_ = qMax(minTemperatureKelvin_, temperatureKelvin_ + dT);
}
