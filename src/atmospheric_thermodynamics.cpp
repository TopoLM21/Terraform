#include "atmospheric_thermodynamics.h"

#include <QtCore/QtMath>

#include <algorithm>
#include <cmath>

namespace {
constexpr double kKelvinOffset = 273.15;
constexpr double kPascalPerAtm = 101325.0;
constexpr double kUniversalGasConstant = 8.314462618;
constexpr double kReferenceCpDry = 1004.0;
constexpr double kDryAirMolarMassKgPerMol = 0.02897;
constexpr double kLatentHeatVaporization = 2.5e6;
constexpr double kMinLapseRate = 0.0015;
constexpr double kMaxLapseRate = 0.012;
}  // namespace

double AtmosphericThermodynamics::meanMolarMassKgPerMol(
    const AtmosphereComposition &composition) {
    const auto fractions = composition.fractions();
    if (fractions.isEmpty()) {
        return kDryAirMolarMassKgPerMol;
    }

    double totalShare = 0.0;
    double inverseMassSum = 0.0;
    const auto gases = availableGases();
    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }
        const auto it = std::find_if(gases.begin(), gases.end(),
                                     [&fraction](const GasSpec &spec) {
                                         return spec.id == fraction.id;
                                     });
        if (it == gases.end()) {
            continue;
        }
        totalShare += fraction.share;
        inverseMassSum += fraction.share / (it->molarMass / 1000.0);
    }

    if (totalShare <= 0.0 || inverseMassSum <= 0.0) {
        return kDryAirMolarMassKgPerMol;
    }
    // Для массовых долей w_i средняя молярная масса:
    // M = 1 / Σ(w_i / M_i) = (Σ w_i) / Σ(w_i / M_i).
    return totalShare / inverseMassSum;
}

double AtmosphericThermodynamics::specificGasConstant(const AtmosphereComposition &composition) {
    // R_specific = R_universal / M.
    const double molarMass = qMax(1e-6, meanMolarMassKgPerMol(composition));
    return kUniversalGasConstant / molarMass;
}

double AtmosphericThermodynamics::specificHeatCp(const AtmosphereComposition &composition) {
    const auto fractions = composition.fractions();
    if (fractions.isEmpty()) {
        return kReferenceCpDry;
    }

    // Средневзвешенная Cp по массовым долям газов: Cp_mix = Σ(w_i × Cp_i).
    // Каждый газ имеет собственную Cp, вычисленную из γ и молярной массы.
    double cpSum = 0.0;
    double totalShare = 0.0;
    const auto gases = availableGases();
    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }
        const auto it = std::find_if(gases.begin(), gases.end(),
                                     [&fraction](const GasSpec &spec) {
                                         return spec.id == fraction.id;
                                     });
        if (it == gases.end()) {
            continue;
        }
        totalShare += fraction.share;
        cpSum += fraction.share * it->specificHeatCpJPerKgK;
    }

    if (totalShare <= 0.0) {
        return kReferenceCpDry;
    }

    return cpSum / totalShare;
}

double AtmosphericThermodynamics::relativeHumidityEstimate(
    const AtmosphereComposition &composition) {
    const auto fractions = composition.fractions();
    double h2oShare = 0.0;
    for (const auto &fraction : fractions) {
        if (fraction.id == QLatin1String("h2o")) {
            h2oShare = fraction.share;
            break;
        }
    }
    // Грубая оценка относительной влажности из доли водяного пара.
    return qBound(0.0, h2oShare * 4.0, 1.0);
}

double AtmosphericThermodynamics::saturationVaporPressureAtm(double temperatureKelvin) {
    // Формула Тетенса для насыщенного давления водяного пара:
    // e_s(T) = 6.112 * exp(17.67 * (T - 273.15) / (T - 29.65)) [гПа].
    const double temperatureCelsius = temperatureKelvin - kKelvinOffset;
    const double exponent = 17.67 * temperatureCelsius / (temperatureCelsius + 243.5);
    const double saturationHpa = 6.112 * std::exp(exponent);
    const double saturationPascal = saturationHpa * 100.0;
    return saturationPascal / kPascalPerAtm;
}

double AtmosphericThermodynamics::moistAdiabaticLapseRate(
    double temperatureKelvin,
    const AtmosphereComposition &composition,
    double pressureAtm,
    double gravity) {
    if (pressureAtm <= 0.0 || gravity <= 0.0) {
        return 0.0;
    }

    const double cp = specificHeatCp(composition);
    const double rSpecific = specificGasConstant(composition);
    const double dryLapse = dryAdiabaticLapseRate(composition, gravity);

    // Влажный градиент через латентную теплоту:
    // Γ_m = Γ_d * (1 + L_v * q_s / (R * T)) / (1 + L_v^2 * q_s / (c_p * R * T^2)).
    // q_s — массовая доля насыщенного водяного пара.
    const double saturationPressureAtm = saturationVaporPressureAtm(temperatureKelvin);
    const double saturationPressure = saturationPressureAtm * kPascalPerAtm;
    const double totalPressure = qMax(1.0, pressureAtm * kPascalPerAtm);
    const double mixingRatio =
        0.622 * saturationPressure / qMax(1.0, totalPressure - saturationPressure);
    const double relativeHumidity = relativeHumidityEstimate(composition);
    const double saturationMixingRatio = mixingRatio * relativeHumidity;

    const double numerator = 1.0 + (kLatentHeatVaporization * saturationMixingRatio) /
                                       (rSpecific * temperatureKelvin);
    const double denominator =
        1.0 + (kLatentHeatVaporization * kLatentHeatVaporization *
               saturationMixingRatio) /
                  (cp * rSpecific * temperatureKelvin * temperatureKelvin);
    const double moistLapse = dryLapse * (numerator / qMax(1e-6, denominator));
    const double maxLapse = qMin(kMaxLapseRate, dryLapse);
    return qBound(kMinLapseRate, moistLapse, maxLapse);
}

double AtmosphericThermodynamics::dryAdiabaticLapseRate(
    const AtmosphereComposition &composition,
    double gravity) {
    const double cp = specificHeatCp(composition);
    // Γ_d = g / c_p.
    return (gravity > 0.0) ? (gravity / qMax(1.0, cp)) : 0.0;
}

double AtmosphericThermodynamics::adiabaticGradientDlnT(
    const AtmosphereComposition &composition,
    double pressureAtm,
    double temperatureKelvin,
    double gravity) {
    const double lapseRate = moistAdiabaticLapseRate(temperatureKelvin,
                                                     composition,
                                                     pressureAtm,
                                                     gravity);
    // d ln T / d ln P = Γ * R / g (из гидростатики и уравнения состояния).
    const double rSpecific = specificGasConstant(composition);
    if (gravity <= 0.0) {
        return 0.0;
    }
    return (lapseRate * rSpecific) / gravity;
}
