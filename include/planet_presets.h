#pragma once

#include <QtCore/QString>
#include <QtCore/QVector>

struct SurfaceMaterial {
    QString id;
    QString name;
    double albedo;
    double emissivity;
    double thermalConductivity;
    double density;
    double specificHeat;
};

struct PlanetPreset {
    QString name;
    double semiMajorAxis;
    double dayLengthDays;
    QString surfaceMaterialId;
};

QVector<SurfaceMaterial> surfaceMaterials();
QVector<PlanetPreset> solarSystemPresets();
QVector<PlanetPreset> sweetSkyPresets();
