#include "atmospheric_pressure_model.h"

#include "atmospheric_thermodynamics.h"

#include <cmath>

namespace {
constexpr double kUniversalGasConstant = 8.314462618; // Дж/(моль·К)
} // namespace

double AtmosphericPressureModel::pressureAtHeightAtm(double seaLevelPressureAtm,
                                                     double heightMeters,
                                                     double temperatureKelvin,
                                                     const AtmosphereComposition &composition,
                                                     double gravity) {
    if (seaLevelPressureAtm <= 0.0) {
        return 0.0;
    }
    if (temperatureKelvin <= 0.0 || gravity <= 0.0) {
        return seaLevelPressureAtm;
    }

    const double molarMass = AtmosphericThermodynamics::meanMolarMassKgPerMol(composition);
    const double specificGasConstant = kUniversalGasConstant / molarMass;
    const double scaleHeight = (specificGasConstant * temperatureKelvin) / gravity;
    if (scaleHeight <= 0.0) {
        return seaLevelPressureAtm;
    }

    return seaLevelPressureAtm * std::exp(-heightMeters / scaleHeight);
}
