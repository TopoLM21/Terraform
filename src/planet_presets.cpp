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
        {QStringLiteral("Меркурий"), 0.39, 58.6, QStringLiteral("regolith_mercury")},
        {QStringLiteral("Венера"), 0.72, 243.0, QStringLiteral("desert")},
        {QStringLiteral("Земля"), 1.00, 1.0, QStringLiteral("ocean")},
        {QStringLiteral("Луна"), 1.00, 27.3, QStringLiteral("regolith_moon")},
        {QStringLiteral("Марс"), 1.52, 1.03, QStringLiteral("desert")},
        {QStringLiteral("Церрера"), 2.77, 0.38, QStringLiteral("ice")},
    };
}

QVector<PlanetPreset> sweetSkyPresets() {
    return {
        {QStringLiteral("Планета 1"), 0.30, 0.7, QStringLiteral("rocky")},
        {QStringLiteral("Планета 2"), 0.40, 1.4, QStringLiteral("forest")},
        {QStringLiteral("Планета 3"), 0.51, 2.2, QStringLiteral("ocean")},
    };
}
