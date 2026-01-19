#pragma once

#include "atmospheric_layer_state.h"

#include <QVector>

class AtmosphericColumn {
public:
    QVector<AtmosphericLayerState> &layers();
    const QVector<AtmosphericLayerState> &layers() const;

    void resize(int layerCount);

private:
    QVector<AtmosphericLayerState> layers_;
};
