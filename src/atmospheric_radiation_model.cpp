#include "atmospheric_radiation_model.h"

#include "radiative_convective_profile.h"

namespace {
constexpr double kTwoStreamEddingtonFactor = 0.75;
}  // namespace

AtmosphericRadiationModel::AtmosphericRadiationModel(const AtmosphereComposition &composition,
                                                     double pressureAtm,
                                                     double baseTemperatureKelvin,
                                                     double effectiveTemperatureKelvin,
                                                     double surfaceGravity,
                                                     RadiationModelType radiationModelType)
    : composition_(composition),
      pressureAtm_(pressureAtm),
      baseTemperatureKelvin_(baseTemperatureKelvin),
      effectiveTemperatureKelvin_(effectiveTemperatureKelvin),
      surfaceGravity_(surfaceGravity),
      radiationModelType_(radiationModelType) {
    computeOpticalDepths();
}

void AtmosphericRadiationModel::computeOpticalDepths() {
    effectiveOpticalDepth_ = 0.0;
    shortwaveOpticalDepth_ = 0.0;
    if (pressureAtm_ <= 0.0) {
        return;
    }

    // Оптическая толщина определяется интегралом τ = ∫ κ(T, p) ρ dz,
    // где κ — массовая непрозрачность, ρ — плотность, z — высота.
    const RadiativeConvectiveProfile profile(composition_,
                                              pressureAtm_,
                                              baseTemperatureKelvin_,
                                              effectiveTemperatureKelvin_,
                                              surfaceGravity_);
    effectiveOpticalDepth_ = profile.totalOpticalDepthLongwave();
    shortwaveOpticalDepth_ = profile.totalOpticalDepthShortwave();
}

double AtmosphericRadiationModel::effectiveOpticalDepth() const {
    return effectiveOpticalDepth_;
}

double AtmosphericRadiationModel::incomingTransmission() const {
    if (shortwaveOpticalDepth_ <= 0.0) {
        return 1.0;
    }
    // Двухпоточная аппроксимация (Eddington) для коротковолнового диапазона:
    // T_sw ≈ 1 / (1 + 3/4 * τ_sw).
    return 1.0 / (1.0 + kTwoStreamEddingtonFactor * shortwaveOpticalDepth_);
}

double AtmosphericRadiationModel::outgoingTransmission() const {
    if (effectiveOpticalDepth_ <= 0.0) {
        return 1.0;
    }
    // Двухпоточная аппроксимация (Eddington) для длинноволнового диапазона:
    // T_lw ≈ 1 / (1 + 3/4 * τ_lw).
    return 1.0 / (1.0 + kTwoStreamEddingtonFactor * effectiveOpticalDepth_);
}

double AtmosphericRadiationModel::applyIncomingFlux(double flux) const {
    return flux * incomingTransmission();
}

double AtmosphericRadiationModel::applyOutgoingFlux(double flux) const {
    return flux * outgoingTransmission();
}
