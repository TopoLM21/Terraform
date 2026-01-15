#include "emission_layer_model.h"

#include "co2_convective_adjustment.h"

#include <algorithm>
#include <cmath>

namespace {
constexpr double kTwoThirds = 2.0 / 3.0;
}

// Переводим «парниковую непрозрачность» в эквивалентную оптическую толщину.
// Fast и Eddington используют разные обратные формулы: для Fast T = exp(-tau),
// а для Eddington T = 1 / (1 + 0.75 * tau), поэтому инверсия зависит от модели.
double opticalDepthFromGreenhouseOpacity(double greenhouseOpacity, RadiationModelType type) {
    const double boundedOpacity = std::clamp(greenhouseOpacity, 0.0, 0.999);
    if (boundedOpacity <= 0.0) {
        return 0.0;
    }
    const double transmission = std::max(1e-6, 1.0 - boundedOpacity);
    if (type == RadiationModelType::Fast) {
        // Закон Бугера-Ламберта: T = exp(-tau).
        return -std::log(transmission);
    }
    // Двухпотоковая аппроксимация: T = 1 / (1 + 0.75 * tau).
    return (1.0 / transmission - 1.0) / 0.75;
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
