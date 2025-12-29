#pragma once

#include "mollweide_projection.h"
#include "planet_surface_grid.h"

#include <QPoint>
#include <QSize>
#include <QRect>

class SurfaceMapInterpolator {
public:
    SurfaceMapInterpolator(const PlanetSurfaceGrid *grid, const MollweideProjection *projection);

    double interpolateTemperatureForPixel(const QPoint &pixel,
                                          const QSize &imageSize,
                                          int neighborCount,
                                          double power,
                                          bool *insideProjection = nullptr) const;

private:
    bool pixelInsideProjection(const QPoint &pixel,
                               const QSize &imageSize,
                               QPointF *projected) const;

    const PlanetSurfaceGrid *grid_ = nullptr;
    const MollweideProjection *projection_ = nullptr;
};
