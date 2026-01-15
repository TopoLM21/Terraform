#pragma once

#include "planet_presets.h"
#include "radiation_model.h"
#include "planet_surface_grid.h"
#include "rotation_mode.h"
#include "subsurface_temperature_solver.h"

#include <QtCore/QHash>
#include <QtCore/QVector>

struct SurfaceTileTemperatureDefaults {
    double minTemperatureKelvin = 3.0;
    double greenhouseOpacity = 0.0;
    RadiationModelType radiationModelType = RadiationModelType::Fast;
    SubsurfaceModelSettings subsurfaceSettings;
};

struct SurfaceTileTemperatureSettings {
    double segmentSolarConstant = 0.0;
    double dayLengthDays = 0.0;
    RotationMode rotationMode = RotationMode::Normal;
    double declinationDegrees = 0.0;
    double orbitalPhaseRadians = 0.0;
    // Резонанс вращения p:q, безразмерный множитель спиновой фазы.
    int spinOrbitP = 1;
    int spinOrbitQ = 1;
    AtmosphereComposition atmosphere;
    double atmospherePressureAtm = 0.0;
    double surfaceGravity = 0.0;
    double cloudAlbedo = 0.0;
    int spinUpDays = 6;
    int currentHourIndex = 0;
    // Абсолютный индекс времени в часах для непрерывной подсолнечной долготы.
    int currentAbsoluteHourIndex = 0;
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
