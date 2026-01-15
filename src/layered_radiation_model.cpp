#include "layered_radiation_model.h"

#include "radiative_convective_profile.h"

#include <QtGlobal>

#include <cmath>

namespace {
constexpr int kLayerCount = 24;
constexpr double kTwoStreamEddingtonFactor = 0.75;
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
}  // namespace

LayeredRadiationModel::LayeredRadiationModel(const AtmosphereComposition &composition,
                                             double pressureAtm,
                                             double baseTemperatureKelvin,
                                             double effectiveTemperatureKelvin,
                                             double surfaceGravity)
    : composition_(composition),
      pressureAtm_(pressureAtm),
      baseTemperatureKelvin_(baseTemperatureKelvin),
      effectiveTemperatureKelvin_(effectiveTemperatureKelvin),
      surfaceGravity_(surfaceGravity) {
    buildLayers();
}

void LayeredRadiationModel::buildLayers() {
    layers_.clear();
    effectiveOpticalDepth_ = 0.0;
    shortwaveOpticalDepth_ = 0.0;

    if (pressureAtm_ <= 0.0) {
        return;
    }

    const RadiativeConvectiveProfile profile(composition_,
                                              pressureAtm_,
                                              baseTemperatureKelvin_,
                                              effectiveTemperatureKelvin_,
                                              surfaceGravity_,
                                              kLayerCount);
    effectiveOpticalDepth_ = profile.totalOpticalDepthLongwave();
    shortwaveOpticalDepth_ = profile.totalOpticalDepthShortwave();

    layers_.reserve(profile.layers().size());
    for (const auto &profileLayer : profile.layers()) {
        Layer layer;
        layer.pressureAtm = profileLayer.pressureAtm;
        layer.temperatureKelvin = profileLayer.temperatureKelvin;
        layer.opticalDepth = profileLayer.opticalDepthLongwave;
        layer.shortwaveOpticalDepth = profileLayer.opticalDepthShortwave;
        layers_.push_back(layer);
    }
}

double LayeredRadiationModel::effectiveOpticalDepth() const {
    return effectiveOpticalDepth_;
}

double LayeredRadiationModel::incomingTransmission() const {
    if (shortwaveOpticalDepth_ <= 0.0) {
        return 1.0;
    }
    // Двухпоточное приближение для SW (Eddington):
    // T_sw ≈ 1 / (1 + 3/4 * τ_sw).
    return 1.0 / (1.0 + kTwoStreamEddingtonFactor * shortwaveOpticalDepth_);
}

double LayeredRadiationModel::outgoingTransmission() const {
    if (effectiveOpticalDepth_ <= 0.0) {
        return 1.0;
    }
    // LW двухпоточная аппроксимация (Eddington):
    // T_lw ≈ 1 / (1 + 3/4 * τ_lw).
    return 1.0 / (1.0 + kTwoStreamEddingtonFactor * effectiveOpticalDepth_);
}

double LayeredRadiationModel::bottomLayerTemperatureKelvin() const {
    if (layers_.isEmpty()) {
        return baseTemperatureKelvin_;
    }
    // Первый слой соответствует максимальному давлению и описывает нижнюю часть атмосферы.
    return layers_.first().temperatureKelvin;
}

double LayeredRadiationModel::layerLongwaveEmissivity(double opticalDepth) const {
    if (opticalDepth <= 0.0) {
        return 0.0;
    }
    // Серость атмосферы + двухпотоковая аппроксимация: ε = 1 - T_lw,
    // где T_lw ≈ 1 / (1 + 3/4 * τ_lw).
    const double transmission = 1.0 / (1.0 + kTwoStreamEddingtonFactor * opticalDepth);
    return qBound(0.0, 1.0 - transmission, 1.0);
}

LayeredRadiationModel::LongwaveFluxProfile
LayeredRadiationModel::longwaveFluxProfile(double surfaceEmittedFluxWPerM2) const {
    LongwaveFluxProfile profile;
    const int layerCount = layers_.size();
    profile.upwardFluxWPerM2.fill(0.0, layerCount + 1);
    profile.downwardFluxWPerM2.fill(0.0, layerCount + 1);

    if (layerCount == 0) {
        profile.upwardFluxWPerM2[0] = surfaceEmittedFluxWPerM2;
        return profile;
    }

    // Двухпотоковая аппроксимация: каждый слой излучает как серое тело
    // с эмиссивностью ε(τ) в оба направления (up/down) по закону Стефана–Больцмана.
    profile.upwardFluxWPerM2[0] = surfaceEmittedFluxWPerM2;
    for (int i = 0; i < layerCount; ++i) {
        const Layer &layer = layers_.at(i);
        const double emissivity = layerLongwaveEmissivity(layer.opticalDepth);
        const double emission =
            emissivity * kStefanBoltzmannConstant * std::pow(layer.temperatureKelvin, 4.0);
        profile.upwardFluxWPerM2[i + 1] =
            profile.upwardFluxWPerM2[i] * (1.0 - emissivity) + emission;
    }

    profile.downwardFluxWPerM2[layerCount] = 0.0;
    for (int i = layerCount - 1; i >= 0; --i) {
        const Layer &layer = layers_.at(i);
        const double emissivity = layerLongwaveEmissivity(layer.opticalDepth);
        const double emission =
            emissivity * kStefanBoltzmannConstant * std::pow(layer.temperatureKelvin, 4.0);
        profile.downwardFluxWPerM2[i] =
            profile.downwardFluxWPerM2[i + 1] * (1.0 - emissivity) + emission;
    }
    return profile;
}
