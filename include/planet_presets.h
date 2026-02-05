#pragma once

#include "atmosphere_model.h"
#include "rotation_mode.h"
#include "solar_calculator.h"
#include "thermal_conductivity_model.h"

#include <cmath>
#include <optional>
#include <QtCore/QString>
#include <QtCore/QVector>
#include <QtGui/QColor>

struct SurfaceMaterial {
    QString id;
    QString name;
    QColor baseColor;
    double albedo;
    double emissivity;
    // Характерная шероховатость поверхности (м), влияет на турбулентный обмен с атмосферой.
    double roughnessLengthMeters;
    // Базовая теплопроводность (Вт/м·К) при referenceTemperatureKelvin.
    double thermalConductivity;
    double density;
    double specificHeat;
    double heatCapacity;
    ThermalConductivityModel thermalConductivityModel;
};

enum class HeightSourceType {
    Procedural,
    HeightmapEquirectangular
};

struct PlanetPreset {
    QString name;
    double semiMajorAxis;
    double dayLengthDays;
    // Длина года в земных сутках.
    double yearLengthDays;
    double eccentricity;
    double obliquityDegrees;
    double perihelionArgumentDegrees;
    // Масса в массах Земли (M⊕).
    double massEarths;
    // Радиус в километрах.
    double radiusKm;
    QString surfaceMaterialId;
    AtmosphereComposition atmosphere;
    // Идентификатор звездного пресета, вокруг которого обращается планета.
    QString starPresetId;
    // Параметры звезды, соответствующие пресету.
    StellarParameters primaryStar{1.0, 5772.0, 1.0};
    std::optional<StellarParameters> secondaryStar = std::nullopt;
    // Ручная непрозрачность парникового слоя (0..1): не является физической моделью.
    // При включённой атмосферной модели ручное значение игнорируется, чтобы избежать
    // двойного учёта парникового эффекта.
    double greenhouseOpacity = 0.0;
    // Включает ручную непрозрачность поверх атмосферной модели (для проверки гипотез).
    bool manualGreenhouseOnTopOfAtmosphere = false;
    // Отражательная способность облаков (0..1), учитывает даже неводные облака.
    double cloudAlbedo = 0.0;
    // Множитель визуальной непрозрачности облаков в отрисовке (0..2).
    double cloudOpacityBoost = 1.0;
    // Временный режим "без воды и облаков" для последующей детальной модели.
    bool disableWaterAndClouds = false;
    // Масса поверхностной воды (гигатонны, 1 гт = 10^9 кг): используется для оценки океанского
    // покрытия, облачности и альбедо в упрощенной гидросфере.
    double surfaceWaterGigatons = 0.0;
    // Геотермальный поток снизу (Вт/м²) для подповерхностной модели.
    // По умолчанию задан земной порядок (~0.02 Вт/м²) как физически реалистичный фон.
    double geothermalFluxWPerM2 = 0.02;
    // Резонанс спин-орбит p:q — безразмерное отношение спиновой частоты к орбитальной.
    // Значение 1:1 соответствует синхронному вращению, 3:2 — меркурианскому резонансу.
    int spinOrbitP = 1;
    int spinOrbitQ = 1;
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
    // Явно выбранный пользователем режим вращения переопределяет авто-оценку.
    bool rotationModeOverride = false;
    // Параметры бинарной орбиты для визуализации (вокруг барицентра).
    bool hasBinaryOrbit = false;
    double binarySemiMajorAxisAu = 0.0;
    double binaryPeriodDays = 0.0;
    double binaryEccentricity = 0.0;
    double binaryInclinationDegrees = 0.0;
    double binaryArgumentPericenterDegreesA = 0.0;
    double binaryArgumentPericenterDegreesB = 0.0;
    // Минимальная толщина нижнего атмосферного слоя (м) для уплотнения профиля.
    double minBottomLayerThicknessMeters = 250.0;
    // Целевое давление на верхней границе атмосферы (атм): физически это
    // порог, ниже которого атмосфера считается «разреженной» для профиля.
    // 0.001 атм ≈ 1 гПа соответствует нижней стратосфере/верхней тропосфере.
    double minTopPressureAtm = 0.001;
    // Нижняя граница температуры для сверхплотных атмосфер (K), 0 — авто.
    double minDenseAtmosphereTemperatureK = 0.0;
    // Kz задаёт интенсивность вертикального турбулентного обмена (м²/с).
    double verticalWindMixingCoefficient = 1.0;
    // Настройки фоновой звезды для визуализации.
    double starBackgroundAngularDiameterDeg = 0.5;
    double starBackgroundIntensity = 1.0;
    QColor starBackgroundCoreColor = QColor(255, 244, 234);
    QColor starBackgroundEdgeColor = QColor(255, 244, 234);
    bool starBackgroundAutoDiameter = true;
    bool starBackgroundAutoColors = true;
    // Флаг "обитаемости" для пресета (служебный маркер для интерфейса/логики).
    bool isHabitable = false;
    // Эвристический сдвиг для ночного охлаждения поверхности (K).
    double diurnalCoolingBiasK = 0.0;
};

// Допуск 2% от орбитального периода: сглаживает округления и различия между
// солнечными и сидерическими сутками, не маскируя резонансы вращения.
constexpr double kTidalLockToleranceFraction = 0.02;

inline bool isTidallyLockedByPeriods(double spinPeriodDays, double orbitalPeriodDays) {
    if (spinPeriodDays <= 0.0 || orbitalPeriodDays <= 0.0) {
        return false;
    }
    const double toleranceDays = orbitalPeriodDays * kTidalLockToleranceFraction;
    return std::abs(spinPeriodDays - orbitalPeriodDays) <= toleranceDays;
}

QVector<SurfaceMaterial> surfaceMaterials();
QVector<PlanetPreset> solarSystemPresets();
QVector<PlanetPreset> sweetSkyPresets();
