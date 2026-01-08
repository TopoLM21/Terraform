#pragma once

#include "planet_surface_grid.h"

#include <QVector>

struct WindVector {
    double eastMps = 0.0;
    double northMps = 0.0;
};

class WindFieldModel {
public:
    QVector<WindVector> buildField(const PlanetSurfaceGrid &grid,
                                   const QVector<double> &pressureAtm,
                                   const QVector<double> &temperatureK,
                                   double dayLengthSeconds,
                                   int smoothingIterations = 1) const;
};
