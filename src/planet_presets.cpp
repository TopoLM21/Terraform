#include "planet_presets.h"

#include "physics/Units.h"

#include "rotation_period_utils.h"

#include "geothermal_flux_model.h"

#include <QtCore/QPair>
#include <QtCore/QtGlobal>

namespace {
constexpr double kEarthMassKg = 5.972e24;
// Для пресетов используем эталонные параметры звезд: расстояние держим равным 1 а.е.,
// так как реальная дистанция до планеты задаётся отдельно через большую полуось орбиты.
constexpr StellarParameters kSolarPrimary{1.0, 5772.0, 1.0};
const QString kSolarStarPresetId = QStringLiteral("solar_primary");
// Значения для "Сладкого Неба" совпадают с исторически выбранными в UI и сохраняют
// характер освещения системы (двойной красный карлик).
constexpr StellarParameters kSweetSkyPrimary{0.3761, 2576.0, 1.0};
constexpr StellarParameters kSweetSkySecondary{0.3741, 2349.0, 1.0};
const QString kSweetSkyStarPresetId = QStringLiteral("sweet_sky_binary");
constexpr double kDefaultPlanetAgeGyr = 4.5;

double estimateGeothermalFluxWPerM2(const GeothermalFluxModel &model,
                                    double massEarths,
                                    double radiusKm,
                                    double ageGyr,
                                    double radiogenicScale,
                                    double residualScale,
                                    double radiogenicDecayGyr = 8.0,
                                    double residualDecayGyr = 4.0,
                                    double parentMassKg = 0.0,
                                    double semiMajorAxisM = 0.0,
                                    double eccentricity = 0.0,
                                    double k2OverQ = 0.0) {
    GeothermalBodyParameters parameters;
    parameters.massKg = massEarths * kEarthMassKg;
    parameters.radiusM = radiusKm * 1000.0;
    parameters.ageGyr = ageGyr;
    parameters.radiogenicScale = radiogenicScale;
    parameters.residualScale = residualScale;
    parameters.radiogenicDecayGyr = radiogenicDecayGyr;
    parameters.residualDecayGyr = residualDecayGyr;
    parameters.parentMassKg = parentMassKg;
    parameters.semiMajorAxisM = semiMajorAxisM;
    parameters.eccentricity = eccentricity;
    parameters.k2OverQ = k2OverQ;
    return model.surfaceFluxWPerM2(parameters);
}

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
        composition.setMassGigatons(entry.first, totalMassKg * entry.second / PhysicsUnits::kKgPerGigaton);
    }
    return composition;
}

QVector<QPair<QString, double>> normalizeAtmosphereShares(
    const QVector<QPair<QString, double>> &shares) {
    double totalShare = 0.0;
    for (const auto &entry : shares) {
        totalShare += entry.second;
    }
    if (totalShare <= 0.0) {
        return shares;
    }

    // Нормализуем доли, чтобы проценты/доли сходились к 1.0 без ручного пересчета.
    QVector<QPair<QString, double>> normalized;
    normalized.reserve(shares.size());
    for (const auto &entry : shares) {
        normalized.append({entry.first, entry.second / totalShare});
    }
    return normalized;
}
} // namespace

QVector<SurfaceMaterial> surfaceMaterials() {
    // Базовые цвета взяты как характерные видимые оттенки материала, чтобы визуально
    // различать типы поверхности без строгой фотометрии.
    QVector<SurfaceMaterial> materials = {
        {QStringLiteral("rocky"), QStringLiteral("Каменная поверхность"), QColor(150, 120, 90), 0.15,
         0.95, 0.03, 2.5, 2600.0, 800.0, 2080000.0},
        {QStringLiteral("ice"), QStringLiteral("Ледяная поверхность"), QColor(220, 235, 245), 0.6, 0.98,
         0.001, 2.2, 920.0, 2100.0, 1932000.0},
        {QStringLiteral("desert"), QStringLiteral("Песчаная поверхность"), QColor(210, 182, 120), 0.35,
         0.93, 0.02, 0.25, 1600.0, 800.0, 1280000.0},
        {QStringLiteral("ocean"), QStringLiteral("Океан"), QColor(40, 140, 200), 0.06, 0.98, 0.0002, 0.6,
         1000.0, 4180.0, 4180000.0},
        {QStringLiteral("forest"), QStringLiteral("Лес"), QColor(52, 112, 66), 0.12, 0.96, 0.5, 0.3,
         1200.0, 1500.0, 1800000.0},
        {QStringLiteral("metal"), QStringLiteral("Металлическая поверхность"), QColor(140, 140, 150), 0.2,
         0.85, 0.001, 45.0, 7800.0, 500.0, 3900000.0},
        {QStringLiteral("regolith_mercury"), QStringLiteral("Реголит Меркурия"), QColor(135, 125, 115),
         0.12, 0.95, 0.01, 0.4, 1700.0, 750.0, 1275000.0,
         {ThermalConductivityModelType::Linear, 250.0, 2.0e-4}},
        {QStringLiteral("regolith_moon"), QStringLiteral("Лунный реголит"), QColor(135, 125, 115), 0.12,
         0.95, 0.01, 0.05, 1700.0, 800.0, 1360000.0,
         // Рост теплопроводности с температурой описывает радиационный теплообмен
         // между зернами реголита и уплотнение контактов при нагреве.
         {ThermalConductivityModelType::Linear, 250.0, 2.5e-4}},
    };
    return materials;
}

QVector<PlanetPreset> solarSystemPresets() {
    // Пресеты Солнечной системы используют параметры Солнца (R☉ = 1, T = 5772 K),
    // поскольку планеты вращаются вокруг одной звезды.
    const GeothermalFluxModel geothermalModel;
    // Эмпирические коэффициенты подобраны так, чтобы простая модель совпала
    // с характерными значениями для известных тел.
    const double mercuryGeothermalFlux =
        estimateGeothermalFluxWPerM2(geothermalModel, 0.0553, 2439.7, kDefaultPlanetAgeGyr, 0.9, 0.9);
    const double venusGeothermalFlux =
        estimateGeothermalFluxWPerM2(geothermalModel, 0.815, 6051.8, kDefaultPlanetAgeGyr, 1.0, 1.0);
    const double earthGeothermalFlux =
        estimateGeothermalFluxWPerM2(geothermalModel, 1.0, 6371.0, kDefaultPlanetAgeGyr, 1.0, 1.0);
    const double moonGeothermalFlux =
        estimateGeothermalFluxWPerM2(geothermalModel,
                                     0.0123,
                                     1737.4,
                                     kDefaultPlanetAgeGyr,
                                     1.4,
                                     1.4,
                                     8.0,
                                     4.0,
                                     kEarthMassKg,
                                     384400000.0,
                                     0.0549,
                                     1.0e-3);
    const double marsGeothermalFlux =
        estimateGeothermalFluxWPerM2(geothermalModel, 0.107, 3389.5, kDefaultPlanetAgeGyr, 0.76, 0.76);
    const double ceresGeothermalFlux =
        estimateGeothermalFluxWPerM2(geothermalModel, 0.00015, 473.0, kDefaultPlanetAgeGyr, 1.0, 1.0);
    const double mercurySolarDayDays = 176.0;
    const double mercuryOrbitalPeriodDays = 87.97;
    const double mercurySiderealPeriodDays =
        solarToSiderealPeriodDays(mercurySolarDayDays,
                                  mercuryOrbitalPeriodDays,
                                  isRetrogradeRotation(0.03));
    Q_ASSERT(mercurySiderealPeriodDays > 0.0);

    QVector<PlanetPreset> presets = {
        // dayLengthDays: солнечные сутки (длительность солнечного дня), не сидерический период.
        // Например, у Венеры солнечные сутки ~116.75, а 243 дня — сидерическое вращение.
        {QStringLiteral("Меркурий"), 0.39, mercurySolarDayDays, mercuryOrbitalPeriodDays, 0.2056, 0.03,
         29.12, 0.0553, 2439.7,
         QStringLiteral("regolith_mercury"), AtmosphereComposition{}, kSolarStarPresetId,
         kSolarPrimary, std::nullopt, 0.0, false, 0.0, 1.0, false, 0.0, mercuryGeothermalFlux,
         // Резонанс 3:2 задан по сидерическому периоду (~58.65 суток), а не по солнечным суткам.
         3, 2, false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         1001u,
         false,
         true,
         false,
         false,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
        {QStringLiteral("Венера"), 0.72, 116.75, 224.70, 0.0068, 177.36, 54.88, 0.815, 6051.8,
         QStringLiteral("desert"),
         atmosphereByPressureAtm(92.0, 0.815, 6051.8, {{QStringLiteral("co2"), 1.0}}),
         kSolarStarPresetId,
         kSolarPrimary,
         std::nullopt,
         0.0,
         false,
         0.75,
         1.0,
         false,
         0.0,
         venusGeothermalFlux,
         // Свободное вращение: резонанс 1:1 задаётся в сидерическом смысле.
         1,
         1,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         1002u,
         true,
         true,
         false,
         false,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
        {QStringLiteral("Земля"), 1.00, 1.0, 365.25, 0.0167, 23.44, 102.94, 1.0, 6371.0,
         QStringLiteral("rocky"),
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
         kSolarStarPresetId,
         kSolarPrimary,
         std::nullopt,
         0.0,
         false,
         0.0,
         1.0,
         false,
         1.4e9,
         earthGeothermalFlux,
         // Свободное вращение: резонанс 1:1 задаётся в сидерическом смысле.
         1,
         1,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         31003u,
         true,
         true,
         false,
         false,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
        {QStringLiteral("Луна"),
         1.00,
         29.5,
         365.25,
         0.0167, // Эксцентриситет солнечной орбиты, не лунной.
         6.68, 0.0, 0.0123, 1737.4,
         QStringLiteral("regolith_moon"), AtmosphereComposition{}, kSolarStarPresetId,
         kSolarPrimary, std::nullopt, 0.0, false, 0.0, 1.0, false, 0.0, moonGeothermalFlux,
         1,
         1,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         1004u,
         false,
         true,
         false,
         false,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
        {QStringLiteral("Марс"), 1.52, 1.03, 686.98, 0.0934, 25.19, 286.5, 0.107, 3389.5,
         QStringLiteral("desert"),
         atmosphereByPressureAtm(0.006, 0.107, 3389.5, {{QStringLiteral("co2"), 1.0}}),
         kSolarStarPresetId,
         kSolarPrimary,
         std::nullopt,
         0.0,
         false,
         0.0,
         1.0,
         false,
         0.0,
         marsGeothermalFlux,
         // Свободное вращение: резонанс 1:1 задаётся в сидерическом смысле.
         1,
         1,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         31005u,
         true,
         true,
         false,
         false,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
        {QStringLiteral("Церрера"), 2.77, 0.38, 1681.63, 0.0758, 4.0, 73.6, 0.00015, 473.0,
         QStringLiteral("ice"), AtmosphereComposition{}, kSolarStarPresetId, kSolarPrimary,
         std::nullopt, 0.0, false, 0.0, 1.0, false, 0.0, ceresGeothermalFlux, 1, 1, false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         1006u,
         false,
         true,
         false,
         false,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
    };
    for (auto &preset : presets) {
        if (preset.name == QStringLiteral("Венера")) {
            preset.minDenseAtmosphereTemperatureK = 700.0;
        }
    }
    return presets;
}

QVector<PlanetPreset> sweetSkyPresets() {
    // "Сладкое Небо" — бинарная система красных карликов: используем одинаковые
    // звездные параметры для всех трёх планет, чтобы акцент оставался на орбитах.
    return {
        {QStringLiteral("Планета 1"), 0.30, 84.9, 84.9, 0.0003, -10.5, 190.51, 0.578, 5151.0,
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
         kSweetSkyStarPresetId,
         kSweetSkyPrimary,
         kSweetSkySecondary,
         0.0,
         false,
         0.0,
         1.0,
         false,
         0.0,
         0.02,
         1,
         1,
         true,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         2001u,
         true,
         true,
         false,
         true,
         0.05,
         15.5,
         0.0226,
         -10.20,
         275.06,
         95.06,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
        {QStringLiteral("Планета 2"), 0.40, 2.4, 92.40, 0.003, 11.94, 21.12, 0.4317, 4710.0,
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
         kSweetSkyStarPresetId,
         kSweetSkyPrimary,
         kSweetSkySecondary,
         0.0,
         false,
         0.0,
         1.0,
         false,
         0.0,
         0.02,
         1,
         1,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         2002u,
         true,
         true,
         false,
         true,
         0.05,
         15.5,
         0.0226,
         -10.20,
         275.06,
         95.06,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
        {QStringLiteral("Планета 3"), 0.51, 1.4, 133.03, 0.0001, -8.84, 343.60, 0.5173, 4979.0,
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
         kSweetSkyStarPresetId,
         kSweetSkyPrimary,
         kSweetSkySecondary,
         0.0,
         false,
         0.0,
         1.0,
         false,
         0.0,
         0.02,
         1,
         1,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         2003u,
         true,
         true,
         false,
         true,
         0.05,
         15.5,
         0.0226,
         -10.20,
         275.06,
         95.06,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
        {QStringLiteral("Планета 4"), 0.94, 1.9, 478.9, 0.0104, 169.82, 50.99, 2.3259, 9649.0,
         QStringLiteral("ice"),
         atmosphereByPressureAtm(
             656.39,
             2.3259,
             9649.0,
             normalizeAtmosphereShares({{QStringLiteral("n2"), 99.9}, {QStringLiteral("ar"), 1.0}})),
         kSweetSkyStarPresetId,
         kSweetSkyPrimary,
         kSweetSkySecondary,
         0.0,
         false,
         0.0,
         1.0,
         false,
         0.0,
         0.02,
         1,
         1,
         false,
         HeightSourceType::Procedural,
         QString(),
         0.0,
         2004u,
         true,
         true,
         false,
         true,
         0.05,
         15.5,
         0.0226,
         -10.20,
         275.06,
         95.06,
         250.0,
         0.001,
         0.0,
         1.0,
         0.5,
         1.0,
         QColor(255, 244, 234),
         QColor(255, 244, 234),
         true,
         true},
    };
}
