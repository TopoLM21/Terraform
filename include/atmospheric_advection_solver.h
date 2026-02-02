#pragma once

#include "atmospheric_grid_3d.h"
#include "planet_surface_grid.h"

#include <QVector>

class AtmosphericAdvectionSolver {
public:
    void advectLayerWinds(const PlanetSurfaceGrid &grid,
                          AtmosphericGrid3D &atmosphereGrid,
                          double dayLengthSeconds,
                          double dtSeconds,
                          int smoothingIterations = 1) const;
    void advectLayerMoisture(const PlanetSurfaceGrid &grid,
                             AtmosphericGrid3D &atmosphereGrid,
                             double dtSeconds) const;

private:
    void ensureNeighbors(const PlanetSurfaceGrid &grid) const;

    mutable const PlanetSurfaceGrid *cachedGrid_ = nullptr;
    mutable QVector<QVector<int>> neighborIndices_;
};
