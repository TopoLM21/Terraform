#pragma once

#include <QColor>
#include <QVector>

struct HeightColorStop {
    double position = 0.0;
    QColor color;
};

const QVector<HeightColorStop> &heightColorStops();
QColor heightColorForRatio(double ratio);
