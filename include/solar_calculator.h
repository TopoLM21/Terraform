#pragma once

struct StellarParameters {
    double radiusInSolarRadii;
    double temperatureKelvin;
    double distanceInAU;
};

class SolarCalculator {
public:
    static double solarConstant(const StellarParameters &parameters);
};

