#pragma once

#include "surface_point_state.h"

#include <QVector>
#include <QString>

struct SurfacePoint {
    double latitudeDeg = 0.0;
    double longitudeDeg = 0.0;
    double latitudeRadians = 0.0;
    double longitudeRadians = 0.0;
    double sinLatitude = 0.0;
    double cosLatitude = 1.0;
    double radiusKm = 0.0;
    double temperatureK = 0.0;
    double airTemperatureK = 0.0;
    // Состояние хранится отдельно для каждой точки, чтобы учитывать локальную
    // тепловую инерцию и парниковую поправку без глобального пересчёта.
    SurfacePointState state;
    double heightKm = 0.0;
    double pressureAtm = 0.0;
    double windEastMps = 0.0;
    double windNorthMps = 0.0;
    double windSpeedMps = 0.0;
    QVector<double> layerHeightsKm;
    QString materialId;
    QVector<int> neighborIndices;
};
