#include "emission_layer_model.h"

#include "co2_convective_adjustment.h"

#include <algorithm>
#include <cmath>

namespace {
constexpr double kTwoThirds = 2.0 / 3.0;
}

// Температура слоя τ≈1 (если атмосфера тоньше — используем τ=τ_surface).
double emissionLayerTemperature(double surfaceTemperatureKelvin, double opticalDepth) {
    const double tauSurface = std::max(0.0, opticalDepth);
    const double tauEmission = std::min(1.0, tauSurface);
    const double denominator = tauSurface + kTwoThirds;
    if (denominator <= 0.0) {
        return surfaceTemperatureKelvin;
    }
    const double ratio = (tauEmission + kTwoThirds) / denominator;
    const double radiativeTemperature = surfaceTemperatureKelvin * std::pow(ratio, 0.25);
    return Co2ConvectiveAdjustment::adjustedEmissionTemperature(surfaceTemperatureKelvin,
                                                                tauSurface,
                                                                radiativeTemperature);
}
