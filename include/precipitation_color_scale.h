#pragma once

#include <QColor>
#include <QVector>

struct PrecipitationColorStop {
    double position = 0.0;
    QColor color;
};

const QVector<PrecipitationColorStop> &precipitationColorStops();
QColor precipitationColorForRatio(double ratio);
