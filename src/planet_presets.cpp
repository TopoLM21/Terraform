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
        {QStringLiteral("Меркурий"), 0.39, 176.0, 0.2056, 0.03, 29.12,
         QStringLiteral("regolith_mercury"), RotationMode::Normal},
        {QStringLiteral("Венера"), 0.72, 243.0, 0.0068, 177.36, 54.88,
         QStringLiteral("desert"), RotationMode::Normal},
        {QStringLiteral("Земля"), 1.00, 1.0, 0.0167, 23.44, 102.94,
         QStringLiteral("ocean"), RotationMode::Normal},
        {QStringLiteral("Луна"), 1.00, 29.5, 0.0549, 6.68, 0.0,
         QStringLiteral("regolith_moon"), RotationMode::Normal},
        {QStringLiteral("Марс"), 1.52, 1.03, 0.0934, 25.19, 286.5,
         QStringLiteral("desert"), RotationMode::Normal},
        {QStringLiteral("Церрера"), 2.77, 0.38, 0.0758, 4.0, 73.6,
         QStringLiteral("ice"), RotationMode::Normal},
    };
}

QVector<PlanetPreset> sweetSkyPresets() {
    return {
        {QStringLiteral("Планета 1"), 0.30, 0.7, 0.02, 5.0, 0.0,
         QStringLiteral("rocky"), RotationMode::Normal},
        {QStringLiteral("Планета 2"), 0.40, 1.4, 0.0, 12.0, 45.0,
         QStringLiteral("forest"), RotationMode::Normal},
        {QStringLiteral("Планета 3"), 0.51, 2.2, 0.1, 20.0, 90.0,
         QStringLiteral("ocean"), RotationMode::Normal},
    };
}
