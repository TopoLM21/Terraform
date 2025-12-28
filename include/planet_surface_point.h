#pragma once

#include <QVector>
#include <QString>

struct SurfacePoint {
    double latitudeDeg = 0.0;
    double longitudeDeg = 0.0;
    double radiusKm = 0.0;
    double temperatureK = 0.0;
    double heightKm = 0.0;
    QVector<double> layerHeightsKm;
    QString materialId;
    QVector<int> neighborIndices;
};
