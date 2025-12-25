#include "atmosphere_model.h"

#include <QtCore/QtMath>

namespace {
constexpr double kEarthMassKg = 5.9722e24;
constexpr double kGravitationalConstant = 6.67430e-11;
constexpr double kPascalPerAtm = 101325.0;
constexpr double kKgPerTon = 1000.0;
constexpr double kKgPerGigaton = 1.0e12;
} // namespace

void AtmosphereComposition::setMassGigatons(const QString &gasId, double massGigatons) {
    for (auto &entry : m_gases) {
        if (entry.id == gasId) {
            entry.massGigatons = massGigatons;
            return;
        }
    }
    m_gases.push_back({gasId, massGigatons, 0.0});
}

double AtmosphereComposition::massGigatons(const QString &gasId) const {
    for (const auto &entry : m_gases) {
        if (entry.id == gasId) {
            return entry.massGigatons;
        }
    }
    return 0.0;
}

QVector<GasFraction> AtmosphereComposition::fractions() const {
    QVector<GasFraction> result;
    result.reserve(m_gases.size());

    const double total = totalMassGigatons();
    for (const auto &entry : m_gases) {
        GasFraction fraction = entry;
        fraction.share = total > 0.0 ? entry.massGigatons / total : 0.0;
        result.push_back(fraction);
    }
    return result;
}

double AtmosphereComposition::totalMassGigatons() const {
    double sum = 0.0;
    for (const auto &entry : m_gases) {
        sum += entry.massGigatons;
    }
    return sum;
}

double AtmosphereComposition::totalMassTons() const {
    return totalMassGigatons() * 1.0e9;
}

double AtmosphereComposition::totalMassKg() const {
    return totalMassGigatons() * kKgPerGigaton;
}

double AtmosphereComposition::totalPressureAtm(double planetMassEarths, double radiusKm) const {
    return calculatePressureAtmFromKg(totalMassKg(), planetMassEarths, radiusKm);
}

QVector<GasSpec> availableGases() {
    return {
        {QStringLiteral("n2"), QStringLiteral("N₂"), 28.014, false},
        {QStringLiteral("o2"), QStringLiteral("O₂"), 31.998, false},
        {QStringLiteral("co2"), QStringLiteral("CO₂"), 44.009, true},
        {QStringLiteral("ch4"), QStringLiteral("CH₄"), 16.043, true},
        {QStringLiteral("h2o"), QStringLiteral("H₂O"), 18.015, true},
        {QStringLiteral("nh3"), QStringLiteral("Аммиак (NH₃)"), 17.031, true},
        {QStringLiteral("sf6"), QStringLiteral("SF₆"), 146.06, true},
        {QStringLiteral("nf3"), QStringLiteral("NF₃"), 71.003, true},
        {QStringLiteral("ar"), QStringLiteral("Ar"), 39.948, false},
        {QStringLiteral("he"), QStringLiteral("He"), 4.0026, false},
        {QStringLiteral("h2"), QStringLiteral("H₂"), 2.01588, false},
    };
}

double calculatePressureAtm(double massTons, double planetMassEarths, double radiusKm) {
    return calculatePressureAtmFromKg(massTons * kKgPerTon, planetMassEarths, radiusKm);
}

double calculatePressureAtmFromKg(double massKg, double planetMassEarths, double radiusKm) {
    const double radiusMeters = radiusKm * 1000.0;
    const double planetMassKg = planetMassEarths * kEarthMassKg;

    // g = G * M / R^2
    const double surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);

    // P = (m_atm * g) / (4 * π * R^2)
    const double surfaceArea = 4.0 * M_PI * radiusMeters * radiusMeters;
    const double pressurePascal = (massKg * surfaceGravity) / surfaceArea;

    return pressurePascal / kPascalPerAtm;
}

double calculateAtmosphereMassKgFromPressureAtm(double pressureAtm,
                                                double planetMassEarths,
                                                double radiusKm) {
    if (pressureAtm <= 0.0 || planetMassEarths <= 0.0 || radiusKm <= 0.0) {
        return 0.0;
    }

    const double radiusMeters = radiusKm * 1000.0;
    const double planetMassKg = planetMassEarths * kEarthMassKg;

    // g = G * M / R^2
    const double surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);

    // m_atm = (P * 4 * π * R^2) / g, где P задано в Па.
    // Единицы: P [Па], m_atm [кг], R [м], g [м/с^2].
    const double surfaceArea = 4.0 * M_PI * radiusMeters * radiusMeters;
    const double pressurePascal = pressureAtm * kPascalPerAtm;
    return (pressurePascal * surfaceArea) / surfaceGravity;
}
