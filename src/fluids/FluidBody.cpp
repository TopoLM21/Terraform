#include "FluidBody.h"

FluidBody::FluidBody(double surfaceAreaSquareMeters)
    : surfaceAreaSquareMeters_(surfaceAreaSquareMeters) {}

double FluidBody::surfaceAreaSquareMeters() const {
    return surfaceAreaSquareMeters_;
}

void FluidBody::setSurfaceAreaSquareMeters(double surfaceAreaSquareMeters) {
    surfaceAreaSquareMeters_ = surfaceAreaSquareMeters;
}

void FluidBody::addLayer(const FluidLayer &layer) {
    layers_.push_back(layer);
}

void FluidBody::setLayers(const QVector<FluidLayer> &layers) {
    layers_ = layers;
}

const QVector<FluidLayer> &FluidBody::layers() const {
    return layers_;
}

double FluidBody::totalMassKg() const {
    double columnMassSum = 0.0;
    for (const auto &layer : layers_) {
        columnMassSum += layer.columnMassKgPerM2();
    }
    return surfaceAreaSquareMeters_ * columnMassSum;
}

double FluidBody::totalMassGigatons() const {
    return totalMassKg() / PhysicsUnits::kKgPerGigaton;
}
