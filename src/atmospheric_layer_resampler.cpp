#include "atmospheric_layer_resampler.h"

#include <QtMath>

namespace {
constexpr double kEpsilonHeightMeters = 1e-6;

double lerp(double a, double b, double t) {
    return a + (b - a) * t;
}
} // namespace

QVector<AtmosphericLayerState> AtmosphericLayerResampler::resample(
    const QVector<AtmosphericLayerState> &oldLayers,
    const QVector<TargetLayer> &targetLayers) {
    QVector<AtmosphericLayerState> result;
    result.reserve(targetLayers.size());
    for (const auto &target : targetLayers) {
        AtmosphericLayerState newLayer =
            interpolateLayerState(oldLayers, target.heightMeters, target.thicknessMeters);
        newLayer.setHeightMeters(target.heightMeters);
        newLayer.setThicknessMeters(target.thicknessMeters);
        result.push_back(newLayer);
    }

    return result;
}

QVector<AtmosphericLayerState> AtmosphericLayerResampler::resampleUniform(
    const QVector<AtmosphericLayerState> &oldLayers,
    int newLayerCount,
    double topHeightMeters) {
    if (newLayerCount <= 0 || topHeightMeters <= 0.0) {
        return {};
    }

    const double layerThickness =
        topHeightMeters / static_cast<double>(newLayerCount);
    if (layerThickness <= 0.0) {
        return {};
    }

    QVector<TargetLayer> targets;
    targets.reserve(newLayerCount);
    double bottomEdgeMeters = 0.0;
    for (int i = 0; i < newLayerCount; ++i) {
        const double heightMeters = bottomEdgeMeters + layerThickness * 0.5;
        bottomEdgeMeters += layerThickness;
        targets.push_back(TargetLayer{heightMeters, layerThickness});
    }

    return resample(oldLayers, targets);
}

AtmosphericLayerState AtmosphericLayerResampler::interpolateLayerState(
    const QVector<AtmosphericLayerState> &layers,
    double targetHeightMeters,
    double targetThicknessMeters) {
    AtmosphericLayerState result;
    if (layers.isEmpty()) {
        return result;
    }

    int upperIndex = 0;
    while (upperIndex < layers.size() &&
           layers.at(upperIndex).heightMeters() < targetHeightMeters) {
        ++upperIndex;
    }

    int lowerIndex = upperIndex;
    if (upperIndex <= 0) {
        lowerIndex = 0;
        upperIndex = 0;
    } else if (upperIndex >= layers.size()) {
        lowerIndex = layers.size() - 1;
        upperIndex = layers.size() - 1;
    } else {
        lowerIndex = upperIndex - 1;
    }

    const AtmosphericLayerState &lower = layers.at(lowerIndex);
    const AtmosphericLayerState &upper = layers.at(upperIndex);
    const double lowerHeight = lower.heightMeters();
    const double upperHeight = upper.heightMeters();
    const double heightSpan = qMax(kEpsilonHeightMeters, upperHeight - lowerHeight);
    const double t = qBound(0.0, (targetHeightMeters - lowerHeight) / heightSpan, 1.0);

    const auto perMeter = [](double value, double thickness) -> double {
        return (thickness > 0.0) ? (value / thickness) : 0.0;
    };

    const double heatCapacityPerMeter =
        lerp(perMeter(lower.heatCapacityJPerM2K(), lower.thicknessMeters()),
             perMeter(upper.heatCapacityJPerM2K(), upper.thicknessMeters()),
             t);
    const double tauSwPerMeter =
        lerp(perMeter(lower.opticalDepthShortwave(), lower.thicknessMeters()),
             perMeter(upper.opticalDepthShortwave(), upper.thicknessMeters()),
             t);
    const double tauLwPerMeter =
        lerp(perMeter(lower.opticalDepthLongwave(), lower.thicknessMeters()),
             perMeter(upper.opticalDepthLongwave(), upper.thicknessMeters()),
             t);

    const double mixingCoefficient =
        lerp(lower.convectionMixingCoefficient(), upper.convectionMixingCoefficient(), t);

    result.setTemperatureKelvin(
        lerp(lower.temperatureKelvin(), upper.temperatureKelvin(), t));
    result.setPressureAtm(lerp(lower.pressureAtm(), upper.pressureAtm(), t));
    result.setDensityKgPerM3(lerp(lower.densityKgPerM3(), upper.densityKgPerM3(), t));
    result.setWindUMps(lerp(lower.windUMps(), upper.windUMps(), t));
    result.setWindVMps(lerp(lower.windVMps(), upper.windVMps(), t));
    // Теплоёмкость и оптическая толщина масштабируются по толщине слоя.
    result.setHeatCapacityJPerM2K(heatCapacityPerMeter * targetThicknessMeters);
    result.setOpticalDepthShortwave(tauSwPerMeter * targetThicknessMeters);
    result.setOpticalDepthLongwave(tauLwPerMeter * targetThicknessMeters);
    result.setConvectionMixingCoefficient(mixingCoefficient);
    result.setConvectionEnabled(mixingCoefficient > 0.0);

    return result;
}
