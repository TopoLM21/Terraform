#pragma once

#include "planet_surface_point.h"

#include <QVector>

class PlanetSurfaceGrid {
public:
    void setRadiusKm(double radiusKm);
    double radiusKm() const;

    void generateFibonacciPoints(int pointCount);

    int pointCount() const;
    double pointAreaKm2() const;

    const QVector<SurfacePoint> &points() const;
    QVector<SurfacePoint> &points();
    const SurfacePoint *pointAt(int index) const;

private:
    double radiusKm_ = 0.0;
    double pointAreaKm2_ = 0.0;
    QVector<SurfacePoint> points_;
};
