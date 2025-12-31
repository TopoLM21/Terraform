#include "atmospheric_pressure_model.h"

#include <cmath>

namespace {
constexpr double kUniversalGasConstant = 8.314462618; // Дж/(моль·К)
constexpr double kDefaultMolarMassKgPerMol = 0.02897; // Средняя для воздуха, кг/моль.

double molarMassKgPerMol(const QString &gasId) {
    const auto gases = availableGases();
    for (const auto &gas : gases) {
        if (gas.id == gasId) {
            return gas.molarMass * 1e-3;
        }
    }
    return 0.0;
}

double meanMolarMassKgPerMol(const AtmosphereComposition &composition) {
    const auto fractions = composition.fractions();
    double inverseMassSum = 0.0;
    double shareSum = 0.0;
    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }
        const double molarMass = molarMassKgPerMol(fraction.id);
        if (molarMass <= 0.0) {
            continue;
        }
        inverseMassSum += fraction.share / molarMass;
        shareSum += fraction.share;
    }
    if (inverseMassSum <= 0.0 || shareSum <= 0.0) {
        return kDefaultMolarMassKgPerMol;
    }
    const double normalizedInverseMass = inverseMassSum / shareSum;
    return 1.0 / normalizedInverseMass;
}
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

    const double molarMass = meanMolarMassKgPerMol(composition);
    const double specificGasConstant = kUniversalGasConstant / molarMass;
    const double scaleHeight = (specificGasConstant * temperatureKelvin) / gravity;
    if (scaleHeight <= 0.0) {
        return seaLevelPressureAtm;
    }

    return seaLevelPressureAtm * std::exp(-heightMeters / scaleHeight);
}
