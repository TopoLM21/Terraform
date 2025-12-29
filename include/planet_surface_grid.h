#pragma once

#include "planet_surface_point.h"
#include "surface_cell.h"

#include <QVector>

class PlanetSurfaceGrid {
public:
    void setRadiusKm(double radiusKm);
    double radiusKm() const;

    void generateFibonacciPoints(int pointCount);
    void generateIcosahedronGrid(int subdivisionLevel);

    int pointCount() const;
    double pointAreaKm2() const;

    const QVector<SurfacePoint> &points() const;
    QVector<SurfacePoint> &points();
    const SurfacePoint *pointAt(int index) const;

    const QVector<SurfaceCell> &cells() const;
    const SurfaceCell *cellAt(int index) const;

private:
    void rebuildIcosahedronCells(int subdivisionLevel);
    void applyHeightModel();

    double radiusKm_ = 0.0;
    double pointAreaKm2_ = 0.0;
    QVector<SurfacePoint> points_;
    QVector<SurfaceCell> cells_;
};
