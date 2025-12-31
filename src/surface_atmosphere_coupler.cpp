#include "surface_atmosphere_coupler.h"

#include "surface_temperature_state.h"

#include <QtCore/QtMath>

SurfaceAtmosphereCoupler::SurfaceAtmosphereCoupler(
    double heatTransferCoefficientWPerM2K)
    : heatTransferCoefficientWPerM2K_(heatTransferCoefficientWPerM2K) {}

void SurfaceAtmosphereCoupler::setHeatTransferCoefficientWPerM2K(double coefficient) {
    heatTransferCoefficientWPerM2K_ = coefficient;
}

double SurfaceAtmosphereCoupler::heatTransferCoefficientWPerM2K() const {
    return heatTransferCoefficientWPerM2K_;
}

void SurfaceAtmosphereCoupler::exchangeSensibleHeat(SurfaceTemperatureState &surface,
                                                    AtmosphericCellState &atmosphere,
                                                    double dtSeconds) const {
    if (dtSeconds <= 0.0) {
        return;
    }

    const double surfaceTemperature = surface.temperatureKelvin();
    const double airTemperature = atmosphere.airTemperatureKelvin();
    const double surfaceHeatCapacity = surface.topLayerHeatCapacityJPerM2K();
    const double airHeatCapacity = atmosphere.heatCapacityJPerM2K();
    if (surfaceHeatCapacity <= 0.0 || airHeatCapacity <= 0.0) {
        return;
    }

    // Чувствительный теплообмен между поверхностью и воздухом:
    // Q_sensible = h_c * (T_surface - T_air), где h_c в Вт/(м^2·К).
    // h_c можно трактовать как суммарную «связь» с погранслоем,
    // зависящую от ветра, шероховатости и турбулентности.
    const double heatTransfer = heatTransferCoefficientWPerM2K_;
    const double fluxWPerM2 = heatTransfer * (surfaceTemperature - airTemperature);

    // Устойчивость явного обмена: dt < C / h_c для каждой подсистемы.
    const double maxStableDt =
        0.5 * qMin(surfaceHeatCapacity / heatTransfer, airHeatCapacity / heatTransfer);
    const double stableDt = qMin(dtSeconds, maxStableDt);

    const double airDelta = fluxWPerM2 * stableDt / airHeatCapacity;

    surface.applySurfaceFlux(-fluxWPerM2, stableDt);
    atmosphere.setAirTemperatureKelvin(airTemperature + airDelta);
}
