#include "emission_layer_model.h"

#include "co2_convective_adjustment.h"

#include <algorithm>
#include <cmath>

namespace {
constexpr double kTwoThirds = 2.0 / 3.0;
}

// Переводим «парниковую непрозрачность» в эквивалентную оптическую толщину.
// Используем ту же связь, что и в калькуляторе: G = 1 - 1 / (1 + 0.75 * tau).
double opticalDepthFromGreenhouseOpacity(double greenhouseOpacity) {
    const double boundedOpacity = std::clamp(greenhouseOpacity, 0.0, 0.999);
    if (boundedOpacity <= 0.0) {
        return 0.0;
    }
    return (4.0 / 3.0) * (boundedOpacity / (1.0 - boundedOpacity));
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
