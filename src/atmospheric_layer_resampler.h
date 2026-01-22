#pragma once

#include "atmospheric_layer_state.h"

#include <QVector>

class AtmosphericLayerResampler {
public:
    struct TargetLayer {
        double heightMeters = 0.0;
        double thicknessMeters = 0.0;
    };

    static QVector<AtmosphericLayerState> resample(
        const QVector<AtmosphericLayerState> &oldLayers,
        const QVector<TargetLayer> &targetLayers);

    static QVector<AtmosphericLayerState> resampleUniform(
        const QVector<AtmosphericLayerState> &oldLayers,
        int newLayerCount,
        double topHeightMeters);

private:
    static AtmosphericLayerState interpolateLayerState(
        const QVector<AtmosphericLayerState> &layers,
        double targetHeightMeters,
        double targetThicknessMeters);
};
