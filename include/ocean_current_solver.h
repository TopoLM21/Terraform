#pragma once

#include "planet_surface_grid.h"

#include <QVector>

class OceanCurrentSolver {
public:
    // Вычисляет океанские течения из ветра и переносит тепло горизонтально
    // между соседними океаническими тайлами по всем подповерхностным слоям.
    void step(PlanetSurfaceGrid &grid,
              double dayLengthSeconds,
              bool isRetrograde,
              double dtSeconds);

private:
    void ensureNeighbors(const PlanetSurfaceGrid &grid) const;

    mutable const PlanetSurfaceGrid *cachedGrid_ = nullptr;
    mutable double cachedRadiusKm_ = -1.0;
    mutable QVector<QVector<int>> neighborIndices_;
};
