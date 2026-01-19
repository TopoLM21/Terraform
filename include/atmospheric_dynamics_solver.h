#pragma once

#include "atmospheric_grid_3d.h"
#include "planet_surface_grid.h"

#include <QVector>

class AtmosphericDynamicsSolver {
public:
    void updateLayerWinds(const PlanetSurfaceGrid &grid,
                          AtmosphericGrid3D &atmosphereGrid,
                          double dayLengthSeconds,
                          bool isRetrograde,
                          double dtSeconds,
                          int smoothingIterations = 1) const;

private:
    void ensureNeighbors(const PlanetSurfaceGrid &grid) const;

    mutable const PlanetSurfaceGrid *cachedGrid_ = nullptr;
    mutable double cachedRadiusKm_ = -1.0;
    mutable QVector<QVector<int>> neighborIndices_;
};
