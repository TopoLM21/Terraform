#include "surface_tile_temperature_calculator.h"

#include "radiation_model_utils.h"
#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kMeanInsolationFactor = 0.25;
constexpr double kTwoThirds = 2.0 / 3.0;

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
}

SurfaceTileTemperatureResult SurfaceTileTemperatureCalculator::initializeSurface(
    PlanetSurfaceGrid &grid,
    const SurfaceTileTemperatureSettings &settings,
    const SurfaceTileTemperatureDefaults &defaults,
    const QHash<QString, SurfaceMaterial> &materialsById,
    const SurfaceMaterial &fallbackMaterial) const {
    SurfaceTileTemperatureResult result;
    const int pointCount = grid.points().size();
    if (pointCount == 0) {
        return result;
    }

    SubsurfaceModelSettings resolvedSubsurfaceSettings = defaults.subsurfaceSettings;
    const SubsurfaceGrid resolvedGrid(resolvedSubsurfaceSettings.layerCount,
                                      resolvedSubsurfaceSettings.topLayerThicknessMeters,
                                      resolvedSubsurfaceSettings.bottomDepthMeters);
    resolvedSubsurfaceSettings.layerCount = resolvedGrid.layerCount();
    resolvedSubsurfaceSettings.bottomDepthMeters = resolvedGrid.bottomDepthMeters();
    const auto &resolvedThicknesses = resolvedGrid.layerThicknessesMeters();
    if (!resolvedThicknesses.isEmpty()) {
        resolvedSubsurfaceSettings.topLayerThicknessMeters = resolvedThicknesses.front();
    }
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
    const double tauSurface =
        opticalDepthFromGreenhouseOpacity(defaults.greenhouseOpacity,
                                          defaults.radiationModelType);
    const double tauEmission = qMin(1.0, qMax(0.0, tauSurface));
    const double greenhouseFactor = (tauEmission + kTwoThirds > 0.0)
        // Инвертируем серую формулу эмиссионного слоя:
        // T_em^4 = T_surface^4 * (τ_em + 2/3) / (τ_surface + 2/3),
        // поэтому T_surface = T_eq * ((τ_surface + 2/3)/(τ_em + 2/3))^(1/4).
        ? std::pow((tauSurface + kTwoThirds) / (tauEmission + kTwoThirds), 0.25)
        : 1.0;
    const double globalTemperature =
        qMax(defaults.minTemperatureKelvin,
             baseEquilibriumTemperature * greenhouseFactor);
    const double meanInsolation =
        settings.hasSolarConstant ? settings.segmentSolarConstant * kMeanInsolationFactor : 0.0;

    double minTemperature = globalTemperature;
    double maxTemperature = globalTemperature;

    for (auto &point : grid.points()) {
        const auto materialIt = materialsById.constFind(point.materialId);
        const SurfaceMaterial material =
            (materialIt != materialsById.cend()) ? *materialIt : fallbackMaterial;
        const double materialAlbedo = qBound(0.0, material.albedo, 1.0);
        // Облака перекрывают поверхность, поэтому берём максимум отражения.
        const double albedo = qMax(materialAlbedo, clampedCloudAlbedo);
        SurfacePointState state(globalTemperature,
                                albedo,
                                defaults.greenhouseOpacity,
                                defaults.radiationModelType,
                                defaults.minTemperatureKelvin,
                                material,
                                defaults.subsurfaceSettings);
        point.state = state;
        point.state.setTemperatureKelvin(globalTemperature);
        point.temperatureK = globalTemperature;
        point.airTemperatureK = globalTemperature;
        result.blendedInsolations.push_back(meanInsolation);
        result.baselineAirTemperatures.push_back(globalTemperature);
    }

    result.hasTemperatureRange = settings.hasSolarConstant && minTemperature <= maxTemperature;
    if (result.hasTemperatureRange) {
        result.minSurfaceTemperatureK = minTemperature;
        result.maxSurfaceTemperatureK = maxTemperature;
    }

    return result;
}
