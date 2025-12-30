#pragma once

#include <QColor>
#include <QVector>

struct TemperatureColorStop {
    double position = 0.0;
    QColor color;
};

const QVector<TemperatureColorStop> &temperatureColorStops();
QColor temperatureColorForRatio(double ratio);
