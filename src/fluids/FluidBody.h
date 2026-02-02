#pragma once

#include <QtCore/QVector>

#include "physics/Units.h"
#include "FluidLayer.h"

class FluidBody {
public:
    FluidBody() = default;
    explicit FluidBody(double surfaceAreaSquareMeters);

    double surfaceAreaSquareMeters() const;
    void setSurfaceAreaSquareMeters(double surfaceAreaSquareMeters);

    void addLayer(const FluidLayer &layer);
    void setLayers(const QVector<FluidLayer> &layers);
    const QVector<FluidLayer> &layers() const;

    // Полная масса по формуле m = A * Σ(ρ * h).
    double totalMassKg() const;
    double totalMassGigatons() const;

private:
    double surfaceAreaSquareMeters_ = 0.0;
    QVector<FluidLayer> layers_;
};
