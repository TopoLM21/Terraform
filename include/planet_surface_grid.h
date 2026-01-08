#pragma once

#include "planet_surface_point.h"
#include "surface_cell.h"
#include "planet_presets.h"

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

    void setHeightSource(HeightSourceType sourceType,
                         const QString &heightmapPath,
                         double heightmapScaleKm,
                         quint32 heightSeed,
                         bool useContinentsHeight,
                         bool hasSeaLevel);

private:
    void rebuildIcosahedronCells(int subdivisionLevel);
    void applyHeightModel();

    double radiusKm_ = 0.0;
    double pointAreaKm2_ = 0.0;
    HeightSourceType heightSourceType_ = HeightSourceType::Procedural;
    QString heightmapPath_;
    double heightmapScaleKm_ = 0.0;
    quint32 heightSeed_ = 0;
    bool useContinentsHeight_ = false;
    bool hasSeaLevel_ = true;
    QVector<SurfacePoint> points_;
    QVector<SurfaceCell> cells_;
};
