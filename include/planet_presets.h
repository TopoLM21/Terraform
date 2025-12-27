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
    // Объемная теплоемкость (J/(м³·K)), используется для теплового баланса в
    // динамическом шаге (нагрев/охлаждение поверхности).
    double heatCapacity;
    // Эффективная суточная тепловая инерция модели для формулы амплитуды
    // дневного хода температур (условные единицы из упрощенной модели).
    // Это не физическая объемная теплоемкость, а подгоночный параметр,
    // задающий сопротивление поверхности суточным колебаниям температуры.
    double dailyThermalInertia;
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
    // Приливная синхронизация задается отдельно: длина суток может совпадать
    // с резонансом вращения (например, 3:2) и не означает жесткую блокировку.
    bool tidallyLocked = false;
};

QVector<SurfaceMaterial> surfaceMaterials();
QVector<PlanetPreset> solarSystemPresets();
QVector<PlanetPreset> sweetSkyPresets();
