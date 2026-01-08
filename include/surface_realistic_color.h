#pragma once

#include "planet_surface_point.h"

#include <QColor>

QColor realisticSurfaceColor(const SurfacePoint &point,
                             double minHeightKm,
                             double maxHeightKm);
