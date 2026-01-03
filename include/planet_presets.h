#pragma once

#include "atmosphere_model.h"
#include "rotation_mode.h"

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
    double heatCapacity;
};

enum class HeightSourceType {
    Procedural,
    HeightmapEquirectangular
};

struct PlanetPreset {
    QString name;
    double semiMajorAxis;
    double dayLengthDays;
    double eccentricity;
    double obliquityDegrees;
    double perihelionArgumentDegrees;
    // Масса в массах Земли (M⊕).
    double massEarths;
    // Радиус в километрах.
    double radiusKm;
    QString surfaceMaterialId;
    AtmosphereComposition atmosphere;
    // Прозрачность парникового слоя (0..1), где 0.99 соответствует плотной атмосфере Венеры.
    double greenhouseOpacity = 0.0;
    // Включает ручную непрозрачность поверх атмосферной модели (для проверки гипотез).
    bool manualGreenhouseOnTopOfAtmosphere = false;
    // Отражательная способность облаков (0..1), учитывает даже неводные облака.
    double cloudAlbedo = 0.0;
    // Приливная синхронизация задается отдельно: длина суток может совпадать
    // с резонансом вращения (например, 3:2) и не означает жесткую блокировку.
    bool tidallyLocked = false;
    // Источник высот поверхности.
    HeightSourceType heightSourceType = HeightSourceType::Procedural;
    // Путь к equirectangular heightmap (строго 2:1), используется только для HeightmapEquirectangular.
    QString heightmapPath;
    // Масштаб высот в километрах: значение 1.0 соответствует диапазону (-1..+1) км.
    double heightmapScaleKm = 0.0;
    // Seed для процедурного рельефа.
    quint32 heightSeed = 0;
    // Управляет разделением суши и океана (морской уровень), не зависит от материала.
    bool hasSeaLevel = false;
    // Процедурный режим "continents" (маска суши + горы/равнины).
    bool useContinentsHeight = false;
};

QVector<SurfaceMaterial> surfaceMaterials();
QVector<PlanetPreset> solarSystemPresets();
QVector<PlanetPreset> sweetSkyPresets();
