#include "layered_radiation_model.h"

#include "radiative_convective_profile.h"

namespace {
constexpr int kLayerCount = 24;
constexpr double kTwoStreamEddingtonFactor = 0.75;
}  // namespace

LayeredRadiationModel::LayeredRadiationModel(const AtmosphereComposition &composition,
                                             double pressureAtm,
                                             double baseTemperatureKelvin,
                                             double surfaceGravity)
    : composition_(composition),
      pressureAtm_(pressureAtm),
      baseTemperatureKelvin_(baseTemperatureKelvin),
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
