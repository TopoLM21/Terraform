#pragma once

#include "planet_presets.h"
#include "planet_surface_grid.h"
#include "rotation_mode.h"
#include "subsurface_temperature_solver.h"

#include <QtCore/QHash>
#include <QtCore/QVector>

struct SurfaceTileTemperatureDefaults {
    double minTemperatureKelvin = 3.0;
    double greenhouseOpacity = 0.0;
    SubsurfaceModelSettings subsurfaceSettings;
};

struct SurfaceTileTemperatureSettings {
    double segmentSolarConstant = 0.0;
    double dayLengthDays = 0.0;
    RotationMode rotationMode = RotationMode::Normal;
    double declinationDegrees = 0.0;
    double orbitalPhaseRadians = 0.0;
    double atmospherePressureAtm = 0.0;
    double cloudAlbedo = 0.0;
    int spinUpDays = 6;
    int currentHourIndex = 0;
    bool hasSolarConstant = false;
    // Используется при смене планеты: показываем «холодный старт» без прогрева.
    bool initializeWithMinTemperatureOnly = false;
};

struct SurfaceTileTemperatureResult {
    bool hasTemperatureRange = false;
    double minSurfaceTemperatureK = 0.0;
    double maxSurfaceTemperatureK = 0.0;
    QVector<double> blendedInsolations;
    QVector<double> baselineAirTemperatures;
    SubsurfaceModelSettings resolvedSubsurfaceSettings;
};

class SurfaceTileTemperatureCalculator {
public:
    SurfaceTileTemperatureResult initializeSurface(PlanetSurfaceGrid &grid,
                                                   const SurfaceTileTemperatureSettings &settings,
                                                   const SurfaceTileTemperatureDefaults &defaults,
                                                   const QHash<QString, SurfaceMaterial> &materialsById,
                                                   const SurfaceMaterial &fallbackMaterial) const;
};
