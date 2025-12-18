#include "solar_calculator.h"

#include <cmath>

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kSolarRadiusMeters = 6.957e8;              // meters
constexpr double kAstronomicalUnitMeters = 1.496e11;         // meters
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;  // W·m−2·K−4

double stellarFluxAtPlanet(const StellarParameters &parameters) {
    const double radiusMeters = parameters.radiusInSolarRadii * kSolarRadiusMeters;
    const double distanceMeters = parameters.distanceInAU * kAstronomicalUnitMeters;

    const double luminosity = 4.0 * kPi * radiusMeters * radiusMeters *
                              kStefanBoltzmannConstant *
                              std::pow(parameters.temperatureKelvin, 4.0);

    return luminosity / (4.0 * kPi * distanceMeters * distanceMeters);
}
}  // namespace

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
