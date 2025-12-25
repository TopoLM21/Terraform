#include "planet_presets.h"

QVector<SurfaceMaterial> surfaceMaterials() {
    return {
        {QStringLiteral("rocky"), QStringLiteral("Каменная поверхность"), 0.15, 0.95, 2.5, 2600.0,
         800.0},
        {QStringLiteral("ice"), QStringLiteral("Ледяная поверхность"), 0.6, 0.98, 2.2, 920.0, 2100.0},
        {QStringLiteral("desert"), QStringLiteral("Песчаная поверхность"), 0.35, 0.93, 0.25, 1600.0,
         800.0},
        {QStringLiteral("ocean"), QStringLiteral("Океан"), 0.06, 0.98, 0.6, 1000.0, 4180.0},
        {QStringLiteral("forest"), QStringLiteral("Лес"), 0.12, 0.96, 0.3, 1200.0, 1500.0},
        {QStringLiteral("metal"), QStringLiteral("Металлическая поверхность"), 0.2, 0.85, 45.0, 7800.0,
         500.0},
        {QStringLiteral("regolith_mercury"), QStringLiteral("Реголит Меркурия"), 0.12, 0.95, 0.4,
         1700.0, 750.0},
        {QStringLiteral("regolith_moon"), QStringLiteral("Лунный реголит"), 0.12, 0.95, 0.02, 1500.0,
         700.0},
    };
}

QVector<PlanetPreset> solarSystemPresets() {
    return {
        // dayLengthDays: солнечные сутки (длительность солнечного дня), не сидерический период
        {QStringLiteral("Меркурий"), 0.39, 176.0, 0.2056, 0.03, 29.12, 0.0553, 2439.7,
         QStringLiteral("regolith_mercury"), false},
        {QStringLiteral("Венера"), 0.72, 243.0, 0.0068, 177.36, 54.88, 0.815, 6051.8,
         QStringLiteral("desert"), false},
        {QStringLiteral("Земля"), 1.00, 1.0, 0.0167, 23.44, 102.94, 1.0, 6371.0,
         QStringLiteral("ocean"), false},
        {QStringLiteral("Луна"), 1.00, 29.5, 0.0549, 6.68, 0.0, 0.0123, 1737.4,
         QStringLiteral("regolith_moon"), false},
        {QStringLiteral("Марс"), 1.52, 1.03, 0.0934, 25.19, 286.5, 0.107, 3389.5,
         QStringLiteral("desert"), false},
        {QStringLiteral("Церрера"), 2.77, 0.38, 0.0758, 4.0, 73.6, 0.00015, 473.0,
         QStringLiteral("ice"), false},
    };
}

QVector<PlanetPreset> sweetSkyPresets() {
    return {
        {QStringLiteral("Планета 1"), 0.30, 84.9, 0.0003, -10.5, 190.51, 0.7, 5200.0,
         QStringLiteral("desert"), true},
        {QStringLiteral("Планета 2"), 0.40, 2.4, 0.003, 11.94, 21.12, 1.1, 6800.0,
         QStringLiteral("desert"), false},
        {QStringLiteral("Планета 3"), 0.51, 1.4, 0.0001, -8.84, 343.60, 0.9, 6000.0,
         QStringLiteral("desert"), false},
    };
}
