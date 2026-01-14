#include "geothermal_flux_model.h"

#include <QtMath>

namespace {
constexpr double kEarthMassKg = 5.972e24;
constexpr double kEarthRadiusM = 6.371e6;
constexpr double kEarthAgeGyr = 4.5;

// Базовые оценки земного вклада в геотермальный поток:
// радиогенный ~0.047 Вт/м², остаточный ~0.04 Вт/м².
constexpr double kEarthRadiogenicFluxWPerM2 = 0.047;
constexpr double kEarthResidualFluxWPerM2 = 0.04;

constexpr double kGravitationalConstant = 6.67430e-11;

double earthFluxToPowerW(double fluxWPerM2) {
    const double surfaceArea = 4.0 * M_PI * kEarthRadiusM * kEarthRadiusM;
    return fluxWPerM2 * surfaceArea;
}

double earthPowerPerKg(double fluxWPerM2) {
    return earthFluxToPowerW(fluxWPerM2) / kEarthMassKg;
}

double applyDecay(double powerPerKg,
                  double ageGyr,
                  double decayGyr,
                  double referenceAgeGyr = kEarthAgeGyr) {
    if (decayGyr <= 0.0) {
        return powerPerKg;
    }

    // Нормируем базовую мощность так, чтобы при возрасте Земли
    // возвращалось земное значение, а затем применяем экспоненциальный спад.
    const double normalizedPower = powerPerKg * qExp(referenceAgeGyr / decayGyr);
    return normalizedPower * qExp(-ageGyr / decayGyr);
}
} // namespace

double GeothermalFluxModel::surfaceFluxWPerM2(const GeothermalBodyParameters &parameters) const {
    if (parameters.radiusM <= 0.0) {
        return 0.0;
    }

    const double powerW =
        radiogenicPowerW(parameters) + residualPowerW(parameters) + tidalPowerW(parameters);
    const double surfaceArea = 4.0 * M_PI * parameters.radiusM * parameters.radiusM;
    return powerW / surfaceArea;
}

double GeothermalFluxModel::radiogenicPowerW(const GeothermalBodyParameters &parameters) const {
    if (parameters.massKg <= 0.0) {
        return 0.0;
    }

    // Q_rad = k_rad * M * exp(-t/τ) с нормировкой на земной поток.
    const double earthRadiogenicPerKg = earthPowerPerKg(kEarthRadiogenicFluxWPerM2);
    const double perKg =
        applyDecay(earthRadiogenicPerKg, parameters.ageGyr, parameters.radiogenicDecayGyr);
    return perKg * parameters.massKg * parameters.radiogenicScale;
}

double GeothermalFluxModel::residualPowerW(const GeothermalBodyParameters &parameters) const {
    if (parameters.massKg <= 0.0) {
        return 0.0;
    }

    // Q_res = k_res * M * exp(-t/τ) с нормировкой на земной поток.
    const double earthResidualPerKg = earthPowerPerKg(kEarthResidualFluxWPerM2);
    const double perKg =
        applyDecay(earthResidualPerKg, parameters.ageGyr, parameters.residualDecayGyr);
    return perKg * parameters.massKg * parameters.residualScale;
}

double GeothermalFluxModel::tidalPowerW(const GeothermalBodyParameters &parameters) const {
    if (parameters.parentMassKg <= 0.0 || parameters.semiMajorAxisM <= 0.0 ||
        parameters.radiusM <= 0.0) {
        return 0.0;
    }

    // Упрощенная формула приливного разогрева:
    // Q_tide ~ (G^2 * M_*^2 * R^5 / a^6) * e^2 * (k2/Q).
    const double parentMassSquared = parameters.parentMassKg * parameters.parentMassKg;
    const double radiusPower = qPow(parameters.radiusM, 5);
    const double axisPower = qPow(parameters.semiMajorAxisM, 6);
    const double eccentricitySquared = parameters.eccentricity * parameters.eccentricity;
    return kGravitationalConstant * kGravitationalConstant * parentMassSquared * radiusPower /
        axisPower * eccentricitySquared * parameters.k2OverQ;
}
