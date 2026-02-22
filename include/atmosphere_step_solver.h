#pragma once

#include "atmosphere_model.h"
#include "atmosphere/EvaporationModel.h"
#include "atmospheric_advection_solver.h"
#include "atmospheric_dynamics_solver.h"
#include "convective_adjustment_solver.h"
#include "layered_radiation_solver.h"
#include "planet_presets.h"
#include "planet_surface_grid.h"
#include "ocean_current_solver.h"
#include "surface/VegetationModel.h"
#include "vertical_wind_mixing_solver.h"
#include "vertical_moisture_mixing_solver.h"

#include <QtCore/QHash>
#include <QtCore/QVector>

class AtmosphereStepSolver {
public:
    struct LayeredStepInput {
        PlanetSurfaceGrid &surfaceGrid;
        AtmosphericGrid3D &atmosphereGrid;
        const QVector<double> &localInsolations;
        const QVector<double> &localCosZeniths;
        const QHash<QString, SurfaceMaterial> &materialsById;
        const SurfaceMaterial &defaultMaterial;
        double cloudShortwaveTransmission = 1.0;
        double heatTransferCoefficientWPerM2K = 0.0;
        double verticalWindMixingCoefficientKz = 1.0;
        double verticalMoistureMixingCoefficientKz = 1.0;
        int logPointIndex = 0;
    };

    AtmosphereStepSolver(const AtmosphereComposition &composition,
                         double gravityMps2,
                         double timeStepSeconds,
                         double dayLengthSeconds,
                         bool isRetrograde,
                         double planetMassEarths = 1.0,
                         double planetRadiusKm = 6371.0);

    void runLayeredStep(const LayeredStepInput &input);

    // Разложение газов и атмосферное убегание за один временной шаг.
    // Изменяет composition и возвращает true, если состав изменился.
    static bool applyGasChemistryAndEscape(AtmosphereComposition &composition,
                                           double timeStepSeconds,
                                           double planetMassEarths,
                                           double planetRadiusKm,
                                           double exosphericTemperatureK);

private:
    const SurfaceMaterial &materialForPoint(const QHash<QString, SurfaceMaterial> &materialsById,
                                            const SurfaceMaterial &defaultMaterial,
                                            const QString &materialId) const;
    void updateLayerPressures(double surfacePressureAtm,
                              QVector<AtmosphericLayerState> &layers) const;

    LayeredRadiationSolver radiationSolver_;
    ConvectiveAdjustmentSolver convectiveSolver_;
    AtmosphericDynamicsSolver dynamicsSolver_;
    AtmosphericAdvectionSolver advectionSolver_;
    VerticalWindMixingSolver verticalWindMixingSolver_;
    VerticalMoistureMixingSolver verticalMoistureMixingSolver_;
    EvaporationModel evaporationModel_;
    VegetationModel vegetationModel_;
    OceanCurrentSolver oceanCurrentSolver_;
    double co2Share_ = 0.0;
    double sf6Share_ = 0.0;
    double nf3Share_ = 0.0;
    double ch4Share_ = 0.0;
    double nh3Share_ = 0.0;
    double h2Share_ = 0.0;
    double gravityMps2_ = 0.0;
    double rSpecific_ = 0.0;
    double specificHeatCp_ = 0.0;
    double timeStepSeconds_ = 0.0;
    double dayLengthSeconds_ = 0.0;
    bool isRetrograde_ = false;
};
