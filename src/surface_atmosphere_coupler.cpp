#include "surface_atmosphere_coupler.h"

#include "surface_point_state.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kVonKarmanConstant = 0.4;
constexpr double kMinWindSpeedMps = 0.1;
constexpr double kMinRoughnessLengthMeters = 1.0e-4;
constexpr double kMinLayerThicknessMeters = 0.1;
} // namespace

SurfaceAtmosphereCoupler::SurfaceAtmosphereCoupler(
    double heatTransferCoefficientWPerM2K)
    : heatTransferCoefficientWPerM2K_(heatTransferCoefficientWPerM2K) {}

void SurfaceAtmosphereCoupler::setHeatTransferCoefficientWPerM2K(double coefficient) {
    heatTransferCoefficientWPerM2K_ = coefficient;
}

double SurfaceAtmosphereCoupler::heatTransferCoefficientWPerM2K() const {
    return heatTransferCoefficientWPerM2K_;
}

void SurfaceAtmosphereCoupler::exchangeHeat(SurfacePointState &surface,
                                            AtmosphericLayerState &lowestLayer,
                                            double roughnessLengthMeters,
                                            double dtSeconds,
                                            double *surfaceAirFluxWPerM2) const {
    if (dtSeconds <= 0.0) {
        if (surfaceAirFluxWPerM2) {
            *surfaceAirFluxWPerM2 = 0.0;
        }
        return;
    }

    const double windSpeed =
        std::hypot(lowestLayer.windUMps(), lowestLayer.windVMps());
    const double longwaveEmissivity =
        longwaveEmissivityFromOpticalDepth(lowestLayer.opticalDepthLongwave());

    double updatedAirTemperature = lowestLayer.temperatureKelvin();
    exchangeHeatInternal(surface,
                         updatedAirTemperature,
                         lowestLayer.heatCapacityJPerM2K(),
                         windSpeed,
                         longwaveEmissivity,
                         roughnessLengthMeters,
                         lowestLayer.thicknessMeters(),
                         dtSeconds,
                         updatedAirTemperature,
                         surfaceAirFluxWPerM2);
    lowestLayer.setTemperatureKelvin(updatedAirTemperature);
}

void SurfaceAtmosphereCoupler::exchangeHeat(SurfacePointState &surface,
                                            AtmosphericCellState &atmosphere,
                                            double windSpeedMps,
                                            double longwaveEmissivity,
                                            double roughnessLengthMeters,
                                            double dtSeconds,
                                            double *surfaceAirFluxWPerM2) const {
    if (dtSeconds <= 0.0) {
        if (surfaceAirFluxWPerM2) {
            *surfaceAirFluxWPerM2 = 0.0;
        }
        return;
    }

    double updatedAirTemperature = atmosphere.airTemperatureKelvin();
    exchangeHeatInternal(surface,
                         updatedAirTemperature,
                         atmosphere.heatCapacityJPerM2K(),
                         windSpeedMps,
                         longwaveEmissivity,
                         roughnessLengthMeters,
                         kMinLayerThicknessMeters,
                         dtSeconds,
                         updatedAirTemperature,
                         surfaceAirFluxWPerM2);
    atmosphere.setAirTemperatureKelvin(updatedAirTemperature);
}

double SurfaceAtmosphereCoupler::sensibleHeatTransferCoefficient(
    double windSpeedMps,
    double roughnessLengthMeters,
    double airHeatCapacityJPerM2K,
    double layerThicknessMeters) const {
    if (airHeatCapacityJPerM2K <= 0.0 || layerThicknessMeters <= 0.0) {
        return heatTransferCoefficientWPerM2K_;
    }

    const double effectiveThickness = qMax(layerThicknessMeters, kMinLayerThicknessMeters);
    // Cp * rho = C / dz для слоя: используем связь с толщиной нижнего слоя.
    const double volumetricHeatCapacity = airHeatCapacityJPerM2K / effectiveThickness;
    const double roughness = qMax(roughnessLengthMeters, kMinRoughnessLengthMeters);
    const double referenceHeight = qMax(1.0, 0.5 * effectiveThickness);
    const double logTerm = qLn(referenceHeight / roughness);
    if (logTerm <= 0.0) {
        return heatTransferCoefficientWPerM2K_;
    }

    // Коэффициент сопротивления для логарифмического профиля ветра:
    // C_d = (kappa / ln(z / z0))^2.
    const double dragCoefficient = std::pow(kVonKarmanConstant / logTerm, 2.0);
    const double windSpeed = qMax(std::abs(windSpeedMps), kMinWindSpeedMps);
    const double transfer = volumetricHeatCapacity * dragCoefficient * windSpeed;
    return qMax(transfer, heatTransferCoefficientWPerM2K_);
}

double SurfaceAtmosphereCoupler::longwaveEmissivityFromOpticalDepth(double opticalDepth) const {
    if (opticalDepth <= 0.0) {
        return 0.0;
    }
    return qBound(0.0, 1.0 - std::exp(-opticalDepth), 1.0);
}

void SurfaceAtmosphereCoupler::exchangeHeatInternal(SurfacePointState &surface,
                                                    double airTemperatureKelvin,
                                                    double airHeatCapacityJPerM2K,
                                                    double windSpeedMps,
                                                    double longwaveEmissivity,
                                                    double roughnessLengthMeters,
                                                    double layerThicknessMeters,
                                                    double dtSeconds,
                                                    double &updatedAirTemperature,
                                                    double *surfaceAirFluxWPerM2) const {
    if (dtSeconds <= 0.0) {
        if (surfaceAirFluxWPerM2) {
            *surfaceAirFluxWPerM2 = 0.0;
        }
        return;
    }

    const double surfaceTemperature = surface.temperatureKelvin();
    const double surfaceHeatCapacity = surface.topLayerHeatCapacityJPerM2K();
    if (surfaceHeatCapacity <= 0.0 || airHeatCapacityJPerM2K <= 0.0) {
        if (surfaceAirFluxWPerM2) {
            *surfaceAirFluxWPerM2 = 0.0;
        }
        return;
    }

    const double sensibleTransfer = sensibleHeatTransferCoefficient(windSpeedMps,
                                                                     roughnessLengthMeters,
                                                                     airHeatCapacityJPerM2K,
                                                                     layerThicknessMeters);
    const double sensibleFluxWPerM2 = sensibleTransfer * (surfaceTemperature - airTemperatureKelvin);

    const double surfaceEmissivity = qBound(0.0, surface.emissivity(), 1.0);
    const double airEmissivity = qBound(0.0, longwaveEmissivity, 1.0);
    // Лучистый обмен между поверхностью и нижним слоем:
    // Q_LW = ε_s * ε_air * σ (T_s^4 - T_air^4).
    const double longwaveFluxWPerM2 =
        surfaceEmissivity * airEmissivity * kStefanBoltzmannConstant *
        (std::pow(surfaceTemperature, 4.0) - std::pow(airTemperatureKelvin, 4.0));

    // Эффективный коэффициент связи для устойчивого явного шага:
    // h_eff = h_sensible + dQ_LW/dT ≈ h_sensible + 4 σ ε_s ε_air T^3.
    const double referenceTemperature = qMax(surfaceTemperature, airTemperatureKelvin);
    const double longwaveDerivative =
        4.0 * kStefanBoltzmannConstant * surfaceEmissivity * airEmissivity *
        std::pow(referenceTemperature, 3.0);
    const double effectiveTransfer = sensibleTransfer + qMax(0.0, longwaveDerivative);

    const double maxStableDt = (effectiveTransfer > 0.0)
        ? 0.5 * qMin(surfaceHeatCapacity / effectiveTransfer,
                     airHeatCapacityJPerM2K / effectiveTransfer)
        : dtSeconds;
    const double stableDt = qMin(dtSeconds, maxStableDt);

    // Итоговый поток (W/м²) считаем положительным при переносе энергии
    // от поверхности к воздуху.
    const double totalFluxWPerM2 = sensibleFluxWPerM2 + longwaveFluxWPerM2;
    const double airDelta = totalFluxWPerM2 * stableDt / airHeatCapacityJPerM2K;

    surface.applySurfaceFlux(-totalFluxWPerM2, stableDt);
    updatedAirTemperature = airTemperatureKelvin + airDelta;
    if (surfaceAirFluxWPerM2) {
        *surfaceAirFluxWPerM2 = totalFluxWPerM2;
    }
}
