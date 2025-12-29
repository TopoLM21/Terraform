#pragma once

#include <QVector>

struct SurfaceVertex {
    double latitudeDeg = 0.0;
    double longitudeDeg = 0.0;
};

struct SurfaceCell {
    QVector<SurfaceVertex> polygon;
    int pointIndex = -1;
};
