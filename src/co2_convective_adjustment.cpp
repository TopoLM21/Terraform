#include "co2_convective_adjustment.h"

#include <algorithm>
#include <cmath>

namespace {
constexpr double kCo2GravityMPerS2 = 8.87;
constexpr double kCo2GasConstantJPerKgK = 189.0;
constexpr double kCo2Gamma = 1.28;
}  // namespace

double Co2ConvectiveAdjustment::dryAdiabaticLapseRate() {
    // Γ_d = g / c_p, где c_p = R / (1 - 1/γ).
    const double cp = kCo2GasConstantJPerKgK / (1.0 - 1.0 / kCo2Gamma);
    return kCo2GravityMPerS2 / std::max(1.0, cp);
}

double Co2ConvectiveAdjustment::adjustedEmissionTemperature(
    double surfaceTemperatureKelvin,
    double opticalDepth,
    double radiativeEmissionTemperatureKelvin) {
    const double tauSurface = std::max(0.0, opticalDepth);
    if (tauSurface <= 0.0 || surfaceTemperatureKelvin <= 0.0) {
        return radiativeEmissionTemperatureKelvin;
    }

    // Связь τ и давления: τ ~ p, поэтому p/p_s = τ/τ_s.
    // Адиабатический профиль: T(τ) = T_s * (τ/τ_s)^(R/c_p).
    const double cp = kCo2GasConstantJPerKgK / (1.0 - 1.0 / kCo2Gamma);
    const double adiabaticExponent = kCo2GasConstantJPerKgK / std::max(1.0, cp);
    const double tauEmission = std::min(1.0, tauSurface);
    const double tauRatio = std::max(1e-6, tauEmission / tauSurface);
    const double adiabaticTemperature =
        surfaceTemperatureKelvin * std::pow(tauRatio, adiabaticExponent);

    // Если радиационный профиль круче сухой адиабаты (температура ниже),
    // поднимаем её до конвективной границы.
    return std::max(radiativeEmissionTemperatureKelvin, adiabaticTemperature);
}
