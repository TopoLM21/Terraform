#pragma once

#include "atmosphere_model.h"
#include "rotation_mode.h"
#include "subsurface_temperature_solver.h"

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
    // Использовать температурную зависимость теплопроводности реголита.
    bool temperatureDependentConductivity = false;
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
    // Приливная синхронизация задается отдельно: длина суток может совпадать
    // с резонансом вращения (например, 3:2) и не означает жесткую блокировку.
    bool tidallyLocked = false;
    // Источник высот поверхности.
    HeightSourceType heightSourceType = HeightSourceType::Procedural;
    // Путь к equirectangular heightmap (строго 2:1).
    QString heightmapPath;
    // Масштаб высот в километрах: значение 1.0 соответствует диапазону (-1..+1) км.
    double heightmapScaleKm = 0.0;
    // Seed для процедурного рельефа (для heightmap игнорируется).
    quint32 heightSeed = 0;
    // Процедурный режим "continents" (маска суши + горы/равнины).
    bool useContinentsHeight = false;
    // Профиль тепловых свойств подповерхности (если нужен отличия от материала).
    SubsurfaceModelSettings subsurfaceSettings;
};

// Теплопроводность реголита часто аппроксимируют формулой:
// λ(T) = λ_contact + b * T^3,
// где λ_contact — контактная проводимость зерен, а член b * T^3 описывает
// радиационный перенос через поры.
// Для Луны типично λ ≈ 0.01..0.02 W/(m·K) у поверхности, с ростом до ~0.05..0.1 W/(m·K)
// на глубине ~1-2 м из-за уплотнения; плотность растет от ~1100 до ~1800 kg/m³,
// удельная теплоемкость ~680..800 J/(kg·K).
inline SubsurfaceModelSettings moonRegolithSubsurfaceSettings() {
    SubsurfaceModelSettings settings;
    settings.thermalConductivityByLayer = {0.012, 0.08};
    settings.densityByLayer = {1100.0, 1800.0};
    settings.specificHeatByLayer = {680.0, 800.0};
    settings.temperatureDependentConductivity = true;
    settings.regolithLambdaContact = 0.01;
    settings.regolithRadiativeCoefficientB = 2.5e-10;
    return settings;
}

inline SubsurfaceModelSettings mercuryRegolithSubsurfaceSettings() {
    SubsurfaceModelSettings settings;
    settings.thermalConductivityByLayer = {0.02, 0.12};
    settings.densityByLayer = {1400.0, 2000.0};
    settings.specificHeatByLayer = {700.0, 850.0};
    settings.temperatureDependentConductivity = true;
    settings.regolithLambdaContact = 0.018;
    settings.regolithRadiativeCoefficientB = 3.0e-10;
    return settings;
}

QVector<SurfaceMaterial> surfaceMaterials();
QVector<PlanetPreset> solarSystemPresets();
QVector<PlanetPreset> sweetSkyPresets();
