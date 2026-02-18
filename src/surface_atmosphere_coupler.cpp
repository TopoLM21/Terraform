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
                                            double *surfaceAirFluxWPerM2,
                                            double surfaceResistanceM2KPerW) const {
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
                         surfaceAirFluxWPerM2,
                         surfaceResistanceM2KPerW);
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
                                                    double *surfaceAirFluxWPerM2,
                                                    double surfaceResistanceM2KPerW) const {
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

    // Длинноволновый обмен между поверхностью и атмосферой полностью рассчитывается
    // в LayeredRadiationSolver (двухпоточная модель Эддингтона). Каплер отвечает
    // только за турбулентный (сенсибл) теплообмен, чтобы избежать двойного учёта.
    Q_UNUSED(longwaveEmissivity);

    // Коэффициент связи для устойчивого явного шага: только сенсибл.
    // При наличии дополнительного теплового сопротивления (например, снег)
    // последовательное сопротивление: 1/h_eff = 1/h_atm + R_surface.
    const double effectiveTransfer =
        (surfaceResistanceM2KPerW > 0.0 && sensibleTransfer > 0.0)
            ? 1.0 / (1.0 / sensibleTransfer + surfaceResistanceM2KPerW)
            : sensibleTransfer;

    const double maxStableDt = (effectiveTransfer > 0.0)
        ? 0.5 * qMin(surfaceHeatCapacity / effectiveTransfer,
                     airHeatCapacityJPerM2K / effectiveTransfer)
        : dtSeconds;
    const double stableDt = qMin(dtSeconds, maxStableDt);

    // Итоговый поток (W/м²) считаем положительным при переносе энергии
    // от поверхности к воздуху. Используем effectiveTransfer (с учётом
    // снежной изоляции) вместо чистого sensibleTransfer.
    const double totalFluxWPerM2 = effectiveTransfer * (surfaceTemperature - airTemperatureKelvin);
    const double airDelta = totalFluxWPerM2 * stableDt / airHeatCapacityJPerM2K;

    surface.applySurfaceFlux(-totalFluxWPerM2, stableDt);
    updatedAirTemperature = airTemperatureKelvin + airDelta;
    if (surfaceAirFluxWPerM2) {
        *surfaceAirFluxWPerM2 = totalFluxWPerM2;
    }
}
