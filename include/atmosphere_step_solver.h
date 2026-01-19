#pragma once

#include "atmosphere_model.h"
#include "atmospheric_advection_solver.h"
#include "atmospheric_dynamics_solver.h"
#include "convective_adjustment_solver.h"
#include "layered_radiation_solver.h"
#include "planet_presets.h"
#include "planet_surface_grid.h"

#include <QtCore/QHash>
#include <QtCore/QVector>

class AtmosphereStepSolver {
public:
    struct LayeredStepInput {
        PlanetSurfaceGrid &surfaceGrid;
        AtmosphericGrid3D &atmosphereGrid;
        const QVector<double> &localInsolations;
        const QHash<QString, SurfaceMaterial> &materialsById;
        const SurfaceMaterial &defaultMaterial;
        double cloudShortwaveTransmission = 1.0;
        double heatTransferCoefficientWPerM2K = 0.0;
    };

    AtmosphereStepSolver(const AtmosphereComposition &composition,
                         double gravityMps2,
                         double timeStepSeconds,
                         double dayLengthSeconds,
                         bool isRetrograde);

    void runLayeredStep(const LayeredStepInput &input);

private:
    const SurfaceMaterial &materialForPoint(const QHash<QString, SurfaceMaterial> &materialsById,
                                            const SurfaceMaterial &defaultMaterial,
                                            const QString &materialId) const;

    LayeredRadiationSolver radiationSolver_;
    ConvectiveAdjustmentSolver convectiveSolver_;
    AtmosphericDynamicsSolver dynamicsSolver_;
    AtmosphericAdvectionSolver advectionSolver_;
    double timeStepSeconds_ = 0.0;
    double dayLengthSeconds_ = 0.0;
    bool isRetrograde_ = false;
};
