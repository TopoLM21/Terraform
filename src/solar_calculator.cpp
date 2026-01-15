#include "solar_calculator.h"

#include <cmath>

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kSolarRadiusMeters = 6.957e8;              // meters
constexpr double kAstronomicalUnitMeters = 1.496e11;         // meters
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;  // W·m−2·K−4

double stellarFluxAtPlanet(const StellarParameters &parameters) {
    return SolarCalculator::stellarFluxAtDistance(parameters.radiusInSolarRadii,
                                                  parameters.temperatureKelvin,
                                                  parameters.distanceInAU);
}
}  // namespace

double SolarCalculator::stellarFluxAtDistance(double radiusInSolarRadii,
                                              double temperatureKelvin,
                                              double distanceInAU) {
    if (radiusInSolarRadii <= 0.0 || temperatureKelvin <= 0.0 || distanceInAU <= 0.0) {
        return 0.0;
    }
    const double radiusMeters = radiusInSolarRadii * kSolarRadiusMeters;
    const double distanceMeters = distanceInAU * kAstronomicalUnitMeters;

    // Светимость звезды: L = 4πR^2σT^4.
    const double luminosity = 4.0 * kPi * radiusMeters * radiusMeters *
                              kStefanBoltzmannConstant *
                              std::pow(temperatureKelvin, 4.0);

    // Поток на расстоянии r: F = L / (4πr^2).
    return luminosity / (4.0 * kPi * distanceMeters * distanceMeters);
}

double SolarCalculator::solarConstant(const StellarParameters &parameters) {
    return stellarFluxAtPlanet(parameters);
}

double SolarCalculator::solarConstant(const BinarySystemParameters &parameters) {
    double totalFlux = stellarFluxAtPlanet(parameters.primary);

    if (parameters.secondary) {
        totalFlux += stellarFluxAtPlanet(*parameters.secondary);
    }

    return totalFlux;
}
