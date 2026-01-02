#include "planet_presets.h"

#include <QtCore/QPair>

namespace {
constexpr double kKgPerGigaton = 1.0e12;

AtmosphereComposition atmosphereByPressureAtm(double pressureAtm,
                                              double planetMassEarths,
                                              double radiusKm,
                                              const QVector<QPair<QString, double>> &shares) {
    AtmosphereComposition composition;
    if (pressureAtm <= 0.0 || shares.isEmpty()) {
        return composition;
    }

    const double totalMassKg =
        calculateAtmosphereMassKgFromPressureAtm(pressureAtm, planetMassEarths, radiusKm);
    if (totalMassKg <= 0.0) {
        return composition;
    }

    for (const auto &entry : shares) {
        composition.setMassGigatons(entry.first, totalMassKg * entry.second / kKgPerGigaton);
    }
    return composition;
}
} // namespace

QVector<SurfaceMaterial> surfaceMaterials() {
    return {
        {QStringLiteral("rocky"), QStringLiteral("Каменная поверхность"), 0.15, 0.95, 2.5, 2600.0,
         800.0, 2080000.0},
        {QStringLiteral("ice"), QStringLiteral("Ледяная поверхность"), 0.6, 0.98, 2.2, 920.0, 2100.0,
         1932000.0},
        {QStringLiteral("desert"), QStringLiteral("Песчаная поверхность"), 0.35, 0.93, 0.25, 1600.0,
         800.0, 1280000.0},
        {QStringLiteral("ocean"), QStringLiteral("Океан"), 0.06, 0.98, 0.6, 1000.0, 4180.0,
         4180000.0},
        {QStringLiteral("forest"), QStringLiteral("Лес"), 0.12, 0.96, 0.3, 1200.0, 1500.0,
         1800000.0},
        {QStringLiteral("metal"), QStringLiteral("Металлическая поверхность"), 0.2, 0.85, 45.0, 7800.0,
         500.0, 3900000.0},
        {QStringLiteral("regolith_mercury"), QStringLiteral("Реголит Меркурия"), 0.12, 0.95, 0.4,
         1700.0, 750.0, 1275000.0},
        {QStringLiteral("regolith_moon"), QStringLiteral("Лунный реголит"), 0.12, 0.95, 0.02, 1500.0,
         700.0, 1050000.0},
    };
}

QVector<PlanetPreset> solarSystemPresets() {
    return {
        // dayLengthDays: солнечные сутки (длительность солнечного дня), не сидерический период.
        // Например, у Венеры солнечные сутки ~116.75, а 243 дня — сидерическое вращение.
        {QStringLiteral("Меркурий"), 0.39, 176.0, 0.2056, 0.03, 29.12, 0.0553, 2439.7,
         QStringLiteral("regolith_mercury"), AtmosphereComposition{}, 0.0, 0.0, false,
         HeightSourceType::Procedural, QString(), 0.0, 1001u, true},
        {QStringLiteral("Венера"), 0.72, 116.75, 0.0068, 177.36, 54.88, 0.815, 6051.8,
         QStringLiteral("desert"),
         atmosphereByPressureAtm(92.0, 0.815, 6051.8, {{QStringLiteral("co2"), 1.0}}),
         0.99,
         0.75,
         false,
         HeightSourceType::Procedural, QString(), 0.0, 1002u, true},
        {QStringLiteral("Земля"), 1.00, 1.0, 0.0167, 23.44, 102.94, 1.0, 6371.0,
         QStringLiteral("ocean"),
         atmosphereByPressureAtm(
             1.0,
             1.0,
             6371.0,
             {
                 {QStringLiteral("n2"), 0.7808},
                 {QStringLiteral("o2"), 0.2095},
                 {QStringLiteral("ar"), 0.0093},
                 {QStringLiteral("co2"), 0.0004},
             }),
         0.0,
         0.0,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         31003u,
         true},
        {QStringLiteral("Луна"), 1.00, 29.5, 0.0549, 6.68, 0.0, 0.0123, 1737.4,
         QStringLiteral("regolith_moon"), AtmosphereComposition{}, 0.0, 0.0, false,
         HeightSourceType::Procedural, QString(), 0.0, 1004u, true},
        {QStringLiteral("Марс"), 1.52, 1.03, 0.0934, 25.19, 286.5, 0.107, 3389.5,
         QStringLiteral("desert"),
         atmosphereByPressureAtm(0.006, 0.107, 3389.5, {{QStringLiteral("co2"), 1.0}}),
         0.0,
         0.0,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         31005u,
         true},
        {QStringLiteral("Церрера"), 2.77, 0.38, 0.0758, 4.0, 73.6, 0.00015, 473.0,
         QStringLiteral("ice"), AtmosphereComposition{}, 0.0, 0.0, false,
         HeightSourceType::Procedural, QString(), 0.0, 1006u, true},
    };
}

QVector<PlanetPreset> sweetSkyPresets() {
    return {
        {QStringLiteral("Планета 1"), 0.30, 84.9, 0.0003, -10.5, 190.51, 0.578, 5151.0,
         QStringLiteral("desert"),
         atmosphereByPressureAtm(
             1.03,
             0.578,
             5151.0,
             {
                 {QStringLiteral("n2"), 0.66},
                 {QStringLiteral("o2"), 0.234},
                 {QStringLiteral("co2"), 0.093},
             }),
         0.0,
         0.0,
         true,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         2001u,
         true},
        {QStringLiteral("Планета 2"), 0.40, 2.4, 0.003, 11.94, 21.12, 0.4317, 4710.0,
         QStringLiteral("desert"),
         atmosphereByPressureAtm(
             0.77,
             0.4317,
             4710.0,
             {
                 {QStringLiteral("n2"), 0.785},
                 {QStringLiteral("o2"), 0.202},
                 {QStringLiteral("co2"), 0.012},
             }),
         0.0,
         0.0,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         2002u,
         true},
        {QStringLiteral("Планета 3"), 0.51, 1.4, 0.0001, -8.84, 343.60, 0.5173, 4979.0,
         QStringLiteral("desert"),
         atmosphereByPressureAtm(
             0.74,
             0.5173,
             4979.0,
             {
                 {QStringLiteral("n2"), 0.788},
                 {QStringLiteral("o2"), 0.201},
                 {QStringLiteral("co2"), 0.011},
             }),
         0.0,
         0.0,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         2003u,
         true},
    };
}
