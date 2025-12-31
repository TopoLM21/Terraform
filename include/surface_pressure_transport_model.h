#pragma once

#include "planet_surface_grid.h"

#include <QVector>

class SurfacePressureTransportModel {
public:
    QVector<double> advectPressure(const PlanetSurfaceGrid &grid,
                                   const QVector<double> &pressureAtm,
                                   const QVector<double> &windEastMps,
                                   const QVector<double> &windNorthMps,
                                   double dtSeconds,
                                   int smoothingIterations,
                                   double minPressureAtm) const;

private:
    void ensureNeighbors(const PlanetSurfaceGrid &grid) const;

    mutable const PlanetSurfaceGrid *cachedGrid_ = nullptr;
    mutable QVector<QVector<int>> neighborIndices_;
};
