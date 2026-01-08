#pragma once

#include "planet_surface_grid.h"

#include <QVector>

class SurfaceAdvectionModel {
public:
    QVector<double> advectTemperature(const PlanetSurfaceGrid &grid,
                                      const QVector<double> &temperatureK,
                                      const QVector<double> &windEastMps,
                                      const QVector<double> &windNorthMps,
                                      double dtSeconds,
                                      int smoothingIterations = 1,
                                      double minTemperatureK = 1.0) const;

private:
    void ensureNeighbors(const PlanetSurfaceGrid &grid) const;

    mutable const PlanetSurfaceGrid *cachedGrid_ = nullptr;
    mutable QVector<QVector<int>> neighborIndices_;
};
