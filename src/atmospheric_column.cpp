#include "atmospheric_column.h"

QVector<AtmosphericLayerState> &AtmosphericColumn::layers() {
    return layers_;
}

const QVector<AtmosphericLayerState> &AtmosphericColumn::layers() const {
    return layers_;
}

void AtmosphericColumn::resize(int layerCount) {
    layers_.resize(layerCount);
}
