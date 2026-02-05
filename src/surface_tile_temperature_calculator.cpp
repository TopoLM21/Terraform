#include "surface_tile_temperature_calculator.h"

#include "radiation_model_utils.h"
#include <QtCore/QLoggingCategory>
#include <QtCore/QtMath>

#include <cmath>

namespace {
Q_LOGGING_CATEGORY(solarRadiationLog, "solar.radiation")
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kMeanInsolationFactor = 0.25;
constexpr double kOceanBottomDepthMeters = 60.0;
constexpr int kOceanLayerCount = 64;
constexpr double kOceanMinTopLayerThicknessMeters = 0.5;

double estimateMeanEquilibriumTemperature(double solarConstant,
                                          double albedo,
                                          double minTemperatureKelvin) {
    if (solarConstant <= 0.0) {
        return minTemperatureKelvin;
    }
    // Средняя инсоляция по сфере: <S> = S0 / 4. Температура равновесия:
    // T_eq = ((S0 * (1 - A) / 4) / sigma)^(1/4).
    const double meanAbsorbedFlux =
        solarConstant * kMeanInsolationFactor * qMax(0.0, 1.0 - albedo);
    if (meanAbsorbedFlux <= 0.0) {
        return minTemperatureKelvin;
    }
    return qMax(minTemperatureKelvin,
                std::pow(meanAbsorbedFlux / kStefanBoltzmannConstant, 0.25));
}

SubsurfaceModelSettings resolveSubsurfaceSettings(const SubsurfaceModelSettings &settings) {
    SubsurfaceModelSettings resolved = settings;
    const SubsurfaceGrid resolvedGrid(resolved.layerCount,
                                      resolved.topLayerThicknessMeters,
                                      resolved.bottomDepthMeters);
    resolved.layerCount = resolvedGrid.layerCount();
    resolved.bottomDepthMeters = resolvedGrid.bottomDepthMeters();
    const auto &resolvedThicknesses = resolvedGrid.layerThicknessesMeters();
    if (!resolvedThicknesses.isEmpty()) {
        resolved.topLayerThicknessMeters = resolvedThicknesses.front();
    }
    return resolved;
}

SubsurfaceModelSettings oceanSubsurfaceSettings(const SubsurfaceModelSettings &settings) {
    // Подбираем океаническую сетку по тем же константам, что и в SolarDisplay.
    SubsurfaceModelSettings ocean = settings;
    ocean.bottomDepthMeters = kOceanBottomDepthMeters;
    ocean.layerCount = kOceanLayerCount;
    ocean.topLayerThicknessMeters =
        qMax(ocean.topLayerThicknessMeters, kOceanMinTopLayerThicknessMeters);
    return ocean;
}
}

SurfaceTileTemperatureResult SurfaceTileTemperatureCalculator::initializeSurface(
    PlanetSurfaceGrid &grid,
    const SurfaceTileTemperatureSettings &settings,
    const SurfaceTileTemperatureDefaults &defaults,
    const QHash<QString, SurfaceMaterial> &materialsById,
    const SurfaceMaterial &fallbackMaterial) const {
    qCInfo(solarRadiationLog) << "Surface tile defaults"
                              << "segmentSolarConstant=" << settings.segmentSolarConstant
                              << "cloudAlbedo=" << settings.cloudAlbedo
                              << "greenhouseOpacity=" << defaults.greenhouseOpacity
                              << "minTemperatureKelvin=" << defaults.minTemperatureKelvin
                              << "diurnalCoolingBiasK=" << defaults.diurnalCoolingBiasK;
    SurfaceTileTemperatureResult result;
    const int pointCount = grid.points().size();
    if (pointCount == 0) {
        return result;
    }

    const SubsurfaceModelSettings resolvedSubsurfaceSettings =
        resolveSubsurfaceSettings(defaults.subsurfaceSettings);
    result.resolvedSubsurfaceSettings = resolvedSubsurfaceSettings;

    result.blendedInsolations.reserve(pointCount);
    result.baselineAirTemperatures.reserve(pointCount);

    const double clampedCloudAlbedo = qBound(0.0, settings.cloudAlbedo, 1.0);
    double materialAlbedoSum = 0.0;
    for (const auto &point : grid.points()) {
        const auto materialIt = materialsById.constFind(point.materialId);
        const SurfaceMaterial material =
            (materialIt != materialsById.cend()) ? *materialIt : fallbackMaterial;
        materialAlbedoSum += qBound(0.0, material.albedo, 1.0);
    }
    const double meanMaterialAlbedo = materialAlbedoSum /
        qMax(1, grid.points().size());
    // Облака перекрывают поверхность, поэтому берём максимум отражения.
    const double globalAlbedo = qMax(meanMaterialAlbedo, clampedCloudAlbedo);

    const double baseEquilibriumTemperature =
        settings.hasSolarConstant
            ? estimateMeanEquilibriumTemperature(settings.segmentSolarConstant,
                                                 globalAlbedo,
                                                 defaults.minTemperatureKelvin)
            : defaults.minTemperatureKelvin;
    // Парниковый множитель применяется ровно один раз в радиационной модели,
    // поэтому стартовую температуру берём по равновесию без повторного
    // «серого» пересчёта (иначе получим двойной разогрев).
    const double globalTemperature =
        qMax(defaults.minTemperatureKelvin,
             baseEquilibriumTemperature);
    const double diurnalCoolingBiasK = qMax(0.0, defaults.diurnalCoolingBiasK);
    // Эти параметры (segmentSolarConstant, cloudAlbedo, greenhouseOpacity,
    // minTemperatureKelvin, diurnalCoolingBiasK) задают единый стартовый
    // adjustedGlobalTemperature: без их изменения стартовая температура не меняется.
    // Грубая поправка на суточное охлаждение: компенсируем отсутствие суточного цикла
    // и инерции грунта в стартовой оценке, чтобы не переоценивать среднюю T.
    const double adjustedGlobalTemperature =
        (diurnalCoolingBiasK > 0.0)
            ? qMax(defaults.minTemperatureKelvin, globalTemperature - diurnalCoolingBiasK)
            : globalTemperature;
    const double meanInsolation =
        settings.hasSolarConstant ? settings.segmentSolarConstant * kMeanInsolationFactor : 0.0;

    double minTemperature = adjustedGlobalTemperature;
    double maxTemperature = adjustedGlobalTemperature;

    for (auto &point : grid.points()) {
        const auto materialIt = materialsById.constFind(point.materialId);
        const SurfaceMaterial material =
            (materialIt != materialsById.cend()) ? *materialIt : fallbackMaterial;
        const double materialAlbedo = qBound(0.0, material.albedo, 1.0);
        // Облака перекрывают поверхность, поэтому берём максимум отражения.
        const double albedo = qMax(materialAlbedo, clampedCloudAlbedo);
        const SubsurfaceModelSettings pointSubsurfaceSettings =
            (material.id == QStringLiteral("ocean"))
                ? resolveSubsurfaceSettings(oceanSubsurfaceSettings(defaults.subsurfaceSettings))
                : resolvedSubsurfaceSettings;
        SurfacePointState state(adjustedGlobalTemperature,
                                albedo,
                                defaults.greenhouseOpacity,
                                defaults.radiationModelType,
                                defaults.minTemperatureKelvin,
                                material,
                                pointSubsurfaceSettings);
        point.state = state;
        point.state.setProfileTemperatureKelvin(adjustedGlobalTemperature);
        point.temperatureK = adjustedGlobalTemperature;
        point.airTemperatureK = adjustedGlobalTemperature;
        result.blendedInsolations.push_back(meanInsolation);
        result.baselineAirTemperatures.push_back(adjustedGlobalTemperature);
    }

    result.hasTemperatureRange = settings.hasSolarConstant && minTemperature <= maxTemperature;
    if (result.hasTemperatureRange) {
        result.minSurfaceTemperatureK = minTemperature;
        result.maxSurfaceTemperatureK = maxTemperature;
    }

    return result;
}
