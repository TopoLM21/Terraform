#include "atmosphere_model.h"

#include <QtCore/QtMath>

#include "physics/Units.h"

namespace {
// Масса Земли для пересчёта планетарных масс в СИ.
constexpr double kEarthMassKg = 5.9722e24;
// Гравитационная постоянная Ньютона, СИ (м^3/(кг·с^2)).
constexpr double kGravitationalConstant = 6.67430e-11;
// 1 стандартная атмосфера в паскалях.
constexpr double kPascalPerAtm = 101325.0;
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
    return totalMassGigatons() * PhysicsUnits::kTonsPerGigaton;
}

double AtmosphereComposition::totalMassKg() const {
    return totalMassGigatons() * PhysicsUnits::kKgPerGigaton;
}

double AtmosphereComposition::totalPressureAtm(double planetMassEarths, double radiusKm) const {
    return calculatePressureAtmFromKg(totalMassKg(), planetMassEarths, radiusKm);
}

QVector<GasSpec> availableGases() {
    // Cp = γ·R_universal / ((γ-1)·M), где R = 8.3145 Дж/(моль·К), M в кг/моль.
    // Значения γ при STP из NIST/CRC Handbook.
    return {
        {QStringLiteral("n2"), QStringLiteral("N₂"), 28.014, false, 1040.0},   // γ=1.400
        {QStringLiteral("o2"), QStringLiteral("O₂"), 31.998, false, 919.0},    // γ=1.395
        {QStringLiteral("co2"), QStringLiteral("CO₂"), 44.009, true, 844.0},   // γ=1.289
        {QStringLiteral("ch4"), QStringLiteral("CH₄"), 16.043, true, 2226.0},  // γ=1.320
        {QStringLiteral("h2o"), QStringLiteral("H₂O"), 18.015, true, 1864.0},  // γ=1.330
        {QStringLiteral("nh3"), QStringLiteral("Аммиак (NH₃)"), 17.031, true, 2060.0}, // γ=1.310
        {QStringLiteral("sf6"), QStringLiteral("SF₆"), 146.06, true, 665.0},   // γ=1.098
        {QStringLiteral("nf3"), QStringLiteral("NF₃"), 71.003, true, 664.0},   // γ=1.214
        {QStringLiteral("ar"), QStringLiteral("Ar"), 39.948, false, 520.0},     // γ=1.667
        {QStringLiteral("he"), QStringLiteral("He"), 4.0026, false, 5193.0},    // γ=1.667
        {QStringLiteral("h2"), QStringLiteral("H₂"), 2.01588, false, 14300.0},  // γ=1.410
    };
}

double calculatePressureAtm(double massTons, double planetMassEarths, double radiusKm) {
    return calculatePressureAtmFromKg(massTons * PhysicsUnits::kKgPerTon, planetMassEarths, radiusKm);
}

double calculatePressureAtmFromKg(double massKg, double planetMassEarths, double radiusKm) {
    const double radiusMeters = radiusKm * 1000.0;
    const double planetMassKg = planetMassEarths * kEarthMassKg;

    // g = G * M / R^2 — поверхностное ускорение свободного падения.
    const double surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);

    // Давление от массы атмосферы, распределённой по сфере:
    // P = (m_atm * g) / (4 * π * R^2).
    const double surfaceArea = 4.0 * M_PI * radiusMeters * radiusMeters;
    const double pressurePascal = (massKg * surfaceGravity) / surfaceArea;

    return pressurePascal / kPascalPerAtm;
}

double calculateCellPressureAtmFromKg(double totalMassKg,
                                      double planetMassEarths,
                                      double radiusKm,
                                      double cellAreaKm2) {
    if (totalMassKg <= 0.0 || planetMassEarths <= 0.0 || radiusKm <= 0.0 || cellAreaKm2 <= 0.0) {
        return 0.0;
    }

    const double radiusMeters = radiusKm * 1000.0;
    const double planetMassKg = planetMassEarths * kEarthMassKg;
    const double cellAreaMeters2 = cellAreaKm2 * 1e6;

    // g = G * M / R^2 — поверхностное ускорение свободного падения.
    const double surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);

    // Масса столба над ячейкой пропорциональна площади: m_cell = m_total * (area_cell / total_area).
    const double totalArea = 4.0 * M_PI * radiusMeters * radiusMeters;
    const double cellMassKg = totalMassKg * (cellAreaMeters2 / totalArea);

    // P = (m_cell * g) / area_cell — давление над ячейкой.
    const double pressurePascal = (cellMassKg * surfaceGravity) / cellAreaMeters2;
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

    // g = G * M / R^2 — поверхностное ускорение свободного падения.
    const double surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);

    // m_atm = (P * 4 * π * R^2) / g, где P задано в Па.
    // Единицы: P [Па], m_atm [кг], R [м], g [м/с^2].
    const double surfaceArea = 4.0 * M_PI * radiusMeters * radiusMeters;
    const double pressurePascal = pressureAtm * kPascalPerAtm;
    return (pressurePascal * surfaceArea) / surfaceGravity;
}
