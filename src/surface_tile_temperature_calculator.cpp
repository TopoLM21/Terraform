#include "surface_tile_temperature_calculator.h"

#include <QtCore/QtMath>

#include <cmath>
#include <limits>

namespace {
constexpr double kPi = 3.14159265358979323846;
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
    const int stepsPerDay = qMax(1, qRound(settings.dayLengthDays * 24.0));
    const double dayLengthSeconds = qMax(0.01, settings.dayLengthDays) * 86400.0;
    const double timeStepSeconds =
        (stepsPerDay > 0) ? (dayLengthSeconds / static_cast<double>(stepsPerDay)) : 0.0;
    const int spinUpDays = qMax(0, settings.spinUpDays);
    const bool allowInsolation =
        !settings.initializeWithMinTemperatureOnly &&
        settings.hasSolarConstant && settings.segmentSolarConstant > 0.0 && timeStepSeconds > 0.0;
    // Если инсоляции нет, не поднимаем стартовый уровень искусственно: температура должна
    // остывать естественно, оставаясь ограниченной только физическим минимумом модели.
    const double visualizationStartTemperature = defaults.minTemperatureKelvin;

    const double declinationRadians = qDegreesToRadians(settings.declinationDegrees);
    const double sinDeclination = std::sin(declinationRadians);
    const double cosDeclination = std::cos(declinationRadians);

    const bool isTidallyLocked = settings.rotationMode == RotationMode::TidalLocked;
    const auto resolveSubstellarLongitude = [isTidallyLocked, stepsPerDay](int hourIndex) {
        if (isTidallyLocked) {
            return 0.0;
        }
        const double phase = 2.0 * kPi * (static_cast<double>(hourIndex) + 0.5) /
                             static_cast<double>(stepsPerDay);
        const double hourAngle = phase - kPi;
        // Субзвёздная долгота — меридиан с нулевым часовым углом.
        return -hourAngle;
    };

    double minTemperature = std::numeric_limits<double>::max();
    double maxTemperature = std::numeric_limits<double>::lowest();

    for (auto &point : grid.points()) {
        const auto materialIt = materialsById.constFind(point.materialId);
        const SurfaceMaterial material =
            (materialIt != materialsById.cend()) ? *materialIt : fallbackMaterial;
        const double materialAlbedo = qBound(0.0, material.albedo, 1.0);
        // Облака перекрывают поверхность, поэтому берём максимум отражения.
        const double albedo = qMax(materialAlbedo, clampedCloudAlbedo);

        SurfacePointState state(visualizationStartTemperature,
                                albedo,
                                defaults.greenhouseOpacity,
                                visualizationStartTemperature,
                                material,
                                defaults.subsurfaceSettings);

        auto applyRadiativeStep = [&](int hourIndex) {
            const double substellarLongitude = resolveSubstellarLongitude(hourIndex);
            const double localHourAngle = point.longitudeRadians - substellarLongitude;
            const double cosZenith =
                point.sinLatitude * sinDeclination +
                point.cosLatitude * cosDeclination * std::cos(localHourAngle);
            // S_inst = S0 * cos(zenith) при освещении, иначе 0.
            const double localInsolation =
                settings.segmentSolarConstant * qMax(0.0, cosZenith);
            // Перенос тепла через глобальную инсоляцию отключён как нефизичный костыль;
            // допускается только атмосферный механизм переноса.
            const double blendedInsolation = localInsolation;
            const double absorbedFlux = state.absorbedFlux(blendedInsolation);
            const double emittedFlux = state.emittedFlux();
            state.updateTemperature(absorbedFlux, emittedFlux, timeStepSeconds);
            return blendedInsolation;
        };

        if (allowInsolation) {
            // Спин-ап прогревает точку до устойчивого цикла перед показом карты.
            for (int day = 0; day < spinUpDays; ++day) {
                for (int step = 0; step < stepsPerDay; ++step) {
                    applyRadiativeStep(step);
                }
            }
        }

        double blendedInsolation = 0.0;
        if (allowInsolation) {
            const int hourIndex = qBound(0, settings.currentHourIndex, stepsPerDay - 1);
            blendedInsolation = applyRadiativeStep(hourIndex);
        }

        point.state = state;
        point.temperatureK = state.temperatureKelvin();
        point.airTemperatureK = point.temperatureK;
        result.blendedInsolations.push_back(blendedInsolation);
        result.baselineAirTemperatures.push_back(point.temperatureK);

        minTemperature = qMin(minTemperature, point.temperatureK);
        maxTemperature = qMax(maxTemperature, point.temperatureK);
    }

    result.hasTemperatureRange = settings.hasSolarConstant && minTemperature <= maxTemperature;
    if (result.hasTemperatureRange) {
        result.minSurfaceTemperatureK = minTemperature;
        result.maxSurfaceTemperatureK = maxTemperature;
    }

    return result;
}
