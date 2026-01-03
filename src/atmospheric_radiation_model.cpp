#include "atmospheric_radiation_model.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kOpticalDepthReferencePressureAtm = 0.1;
constexpr double kOpticalDepthPressureExponent = 2.0;
constexpr double kVenusSurfacePressureAtm = 92.0;
constexpr double kVenusTargetEffectiveTemperatureKelvin = 230.0;
constexpr double kVenusTargetSurfaceTemperatureKelvin = 735.0;
constexpr double kShortwaveToLongwaveRatio = 0.12;
}  // namespace

AtmosphericRadiationModel::AtmosphericRadiationModel(const AtmosphereComposition &composition,
                                                     double pressureAtm,
                                                     double baseTemperatureKelvin)
    : composition_(composition),
      pressureAtm_(pressureAtm),
      baseTemperatureKelvin_(baseTemperatureKelvin) {
    computeOpticalDepths();
}

void AtmosphericRadiationModel::computeOpticalDepths() {
    effectiveOpticalDepth_ = 0.0;
    shortwaveOpticalDepth_ = 0.0;
    if (pressureAtm_ <= 0.0) {
        return;
    }

    // Простая параметризация серой атмосферы:
    // τ(p) = τ_s * (p / p_s)^n, где n≈2.0 и τ_s — оптическая толщина на поверхности.
    // Требуем τ(p=0.1 бар) ≈ 1 ⇒ τ_s = (p_s / 0.1)^n. Далее вводим калибровочный
    // множитель, чтобы при параметрах Венеры получить T_eff ≈ 230 K
    // (OLR_TOA = σ T_eff^4 ≈ 160 Вт/м²).
    const double pressureRatio =
        qMax(0.0, pressureAtm_ / kOpticalDepthReferencePressureAtm);
    const double tauSurfaceReference =
        std::pow(pressureRatio, kOpticalDepthPressureExponent);
    // Калибровка по Венере: берем требуемую эффективную температуру (T_eff),
    // связываем её с OLR через σ T_eff^4 ≈ 160 Вт/м², и подбираем τ_s так, чтобы
    // в серой атмосфере получалась наблюдаемая T_surf ~ 735 K:
    // T_surf = T_eff * (1 + 0.75 * τ_s)^(1/4).
    const double venusTargetTauSurface =
        (std::pow(kVenusTargetSurfaceTemperatureKelvin / kVenusTargetEffectiveTemperatureKelvin,
                  4.0) -
         1.0) /
        0.75;
    const double venusReferenceTauSurface =
        std::pow(kVenusSurfacePressureAtm / kOpticalDepthReferencePressureAtm,
                 kOpticalDepthPressureExponent);
    const double tauCalibration =
        (venusReferenceTauSurface > 0.0)
            ? (venusTargetTauSurface / venusReferenceTauSurface)
            : 0.0;

    effectiveOpticalDepth_ = tauSurfaceReference * tauCalibration;
    // В коротковолновом диапазоне поглощение значительно слабее, поэтому берем долю
    // от длинноволновой оптической толщины.
    shortwaveOpticalDepth_ = effectiveOpticalDepth_ * kShortwaveToLongwaveRatio;
}

double AtmosphericRadiationModel::effectiveOpticalDepth() const {
    return effectiveOpticalDepth_;
}

double AtmosphericRadiationModel::incomingTransmission() const {
    // Закон Бугера-Ламберта: I = I0 * exp(-tau).
    return std::exp(-shortwaveOpticalDepth_);
}

double AtmosphericRadiationModel::outgoingTransmission() const {
    // Эффективная прозрачность в длинноволновом диапазоне.
    return std::exp(-effectiveOpticalDepth_);
}

double AtmosphericRadiationModel::applyIncomingFlux(double flux) const {
    return flux * incomingTransmission();
}

double AtmosphericRadiationModel::applyOutgoingFlux(double flux) const {
    return flux * outgoingTransmission();
}
