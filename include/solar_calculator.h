#pragma once

#include <optional>

struct StellarParameters {
    double radiusInSolarRadii;
    double temperatureKelvin;
    double distanceInAU;
};

struct BinarySystemParameters {
    StellarParameters primary;
    std::optional<StellarParameters> secondary;
};

class SolarCalculator {
public:
    static double solarConstant(const StellarParameters &parameters);
    static double solarConstant(const BinarySystemParameters &parameters);
};
