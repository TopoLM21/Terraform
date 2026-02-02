#include "FluidLayer.h"

FluidLayer::FluidLayer(double densityKgPerM3, double thicknessMeters)
    : densityKgPerM3_(densityKgPerM3)
    , thicknessMeters_(thicknessMeters) {}

double FluidLayer::densityKgPerM3() const {
    return densityKgPerM3_;
}

void FluidLayer::setDensityKgPerM3(double densityKgPerM3) {
    densityKgPerM3_ = densityKgPerM3;
}

double FluidLayer::thicknessMeters() const {
    return thicknessMeters_;
}

void FluidLayer::setThicknessMeters(double thicknessMeters) {
    thicknessMeters_ = thicknessMeters;
}

double FluidLayer::columnMassKgPerM2() const {
    return densityKgPerM3_ * thicknessMeters_;
}

double FluidLayer::pressurePascal(double gravityMps2) const {
    return densityKgPerM3_ * gravityMps2 * thicknessMeters_;
}
