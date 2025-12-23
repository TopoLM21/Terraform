#include "planet_presets.h"

QVector<SurfaceMaterial> surfaceMaterials() {
    return {
        {QStringLiteral("rocky"), QStringLiteral("Каменная поверхность"), 0.15, 0.95},
        {QStringLiteral("ice"), QStringLiteral("Ледяная поверхность"), 0.6, 0.98},
        {QStringLiteral("desert"), QStringLiteral("Песчаная поверхность"), 0.35, 0.93},
        {QStringLiteral("ocean"), QStringLiteral("Океан"), 0.06, 0.98},
        {QStringLiteral("forest"), QStringLiteral("Лес"), 0.12, 0.96},
        {QStringLiteral("metal"), QStringLiteral("Металлическая поверхность"), 0.2, 0.85},
    };
}

QVector<PlanetPreset> solarSystemPresets() {
    return {
        {QStringLiteral("Меркурий"), 0.39, 58.6, QStringLiteral("rocky")},
        {QStringLiteral("Венера"), 0.72, 243.0, QStringLiteral("desert")},
        {QStringLiteral("Земля"), 1.00, 1.0, QStringLiteral("ocean")},
        {QStringLiteral("Луна"), 1.00, 27.3, QStringLiteral("rocky")},
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
