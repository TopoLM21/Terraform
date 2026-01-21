#pragma once

#include "atmospheric_grid_3d.h"
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
    double tileAreaKm2(int pointIndex) const;
    double tileEdgeLengthKm(int pointIndex) const;

    const QVector<SurfacePoint> &points() const;
    QVector<SurfacePoint> &points();
    const SurfacePoint *pointAt(int index) const;

    const QVector<SurfaceCell> &cells() const;
    const SurfaceCell *cellAt(int index) const;

    AtmosphericGrid3D &atmosphericGrid();
    const AtmosphericGrid3D &atmosphericGrid() const;
    void initializeAtmosphericGrid(const AtmosphereComposition &composition,
                                   double planetMassEarths,
                                   double baseTemperatureKelvin,
                                   int layerCount = 0,
                                   double minBottomLayerThicknessMeters = 0.0,
                                   double minTopPressureAtm = 0.0);

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
    AtmosphericGrid3D atmosphericGrid_;
};
