#include "layered_radiation_model.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr int kLayerCount = 12;
constexpr double kOpticalDepthReferencePressureAtm = 0.1;
constexpr double kOpticalDepthPressureExponent = 2.0;
constexpr double kVenusSurfacePressureAtm = 92.0;
constexpr double kVenusTargetEffectiveTemperatureKelvin = 230.0;
constexpr double kVenusTargetSurfaceTemperatureKelvin = 735.0;
constexpr double kShortwaveToLongwaveRatio = 0.12;
constexpr double kTopPressureFloorAtm = 1e-4;
constexpr double kMinPressureRatio = 1e-3;
constexpr double kAdiabaticKappa = 0.286;  // R/Cp для сухого воздуха (серый газ).
}  // namespace

LayeredRadiationModel::LayeredRadiationModel(const AtmosphereComposition &composition,
                                             double pressureAtm,
                                             double baseTemperatureKelvin)
    : composition_(composition),
      pressureAtm_(pressureAtm),
      baseTemperatureKelvin_(baseTemperatureKelvin) {
    buildLayers();
}

void LayeredRadiationModel::buildLayers() {
    layers_.clear();
    effectiveOpticalDepth_ = 0.0;
    shortwaveOpticalDepth_ = 0.0;

    if (pressureAtm_ <= 0.0) {
        return;
    }

    // Простая параметризация серой атмосферы с постоянной κ:
    // τ(p) = τ_s * (p / p_s)^n, где n≈2.0.
    // τ_s подбираем по Венере через формулу серой атмосферы
    // T_surf = T_eff * (1 + 0.75 * τ_s)^(1/4).
    const double pressureRatio =
        qMax(0.0, pressureAtm_ / kOpticalDepthReferencePressureAtm);
    const double tauSurfaceReference =
        std::pow(pressureRatio, kOpticalDepthPressureExponent);
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
    shortwaveOpticalDepth_ = effectiveOpticalDepth_ * kShortwaveToLongwaveRatio;

    const double topPressureAtm =
        qMax(kTopPressureFloorAtm, pressureAtm_ * kMinPressureRatio);
    const double logSurface = std::log(pressureAtm_);
    const double logTop = std::log(topPressureAtm);

    QVector<double> edges;
    edges.reserve(kLayerCount + 1);
    for (int i = 0; i <= kLayerCount; ++i) {
        const double frac = static_cast<double>(i) / kLayerCount;
        const double logP = logSurface + (logTop - logSurface) * frac;
        edges.push_back(std::exp(logP));
    }

    layers_.reserve(kLayerCount);
    for (int i = 0; i < kLayerCount; ++i) {
        const double pLower = edges[i];
        const double pUpper = edges[i + 1];
        const double pMid = std::sqrt(pLower * pUpper);
        const double deltaP = pLower - pUpper;

        Layer layer;
        layer.pressureAtm = pMid;
        // Радиативно-конвективный профиль: сухая адиабата
        // T(p) = T_s * (p/p_s)^(R/Cp), серый газ, κ=const.
        layer.temperatureKelvin =
            baseTemperatureKelvin_ * std::pow(pMid / pressureAtm_, kAdiabaticKappa);
        // Постоянное κ ⇒ τ ∝ p, поэтому распределяем τ_s по давлению.
        layer.opticalDepth = (pressureAtm_ > 0.0)
                                 ? effectiveOpticalDepth_ * (deltaP / pressureAtm_)
                                 : 0.0;
        layer.shortwaveOpticalDepth = layer.opticalDepth * kShortwaveToLongwaveRatio;
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
    // SW ослабляется по слоям законом Бугера-Ламберта: I = I0 * exp(-Στ_sw).
    double transmission = 1.0;
    for (const auto &layer : layers_) {
        transmission *= std::exp(-layer.shortwaveOpticalDepth);
    }
    return transmission;
}

double LayeredRadiationModel::outgoingTransmission() const {
    if (effectiveOpticalDepth_ <= 0.0) {
        return 1.0;
    }
    // LW двухпоточная аппроксимация (Eddington):
    // t_lw ≈ 1 / (1 + 0.75 * τ_layer). Итоговый коэффициент — произведение слоёв.
    double transmission = 1.0;
    for (const auto &layer : layers_) {
        const double layerTransmission = 1.0 / (1.0 + 0.75 * layer.opticalDepth);
        transmission *= layerTransmission;
    }
    return transmission;
}
