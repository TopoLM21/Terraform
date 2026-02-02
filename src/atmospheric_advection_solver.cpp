#include "atmospheric_advection_solver.h"

#include "wind_field_model.h"
#include "atmosphere/EvaporationModel.h"

#include <QVector3D>
#include <QtMath>

#include <cmath>

namespace {
constexpr int kNeighborCount = 6;
constexpr double kHalfPi = 0.5 * M_PI;
constexpr double kMinCosLatitude = 1.0e-5;
constexpr double kSmoothingFactor = 0.3;
constexpr double kMaxWindSpeedMps = 150.0;

QVector3D latLonToCartesian(double latitudeRad, double longitudeRad) {
    const double cosLat = qCos(latitudeRad);
    return QVector3D(static_cast<float>(cosLat * qSin(longitudeRad)),
                     static_cast<float>(qSin(latitudeRad)),
                     static_cast<float>(cosLat * qCos(longitudeRad)));
}

QVector<QVector<int>> buildNeighborIndices(const QVector<SurfacePoint> &points) {
    QVector<QVector<int>> neighbors;
    neighbors.resize(points.size());
    QVector<QVector3D> positions;
    positions.reserve(points.size());
    for (const auto &point : points) {
        positions.push_back(latLonToCartesian(point.latitudeRadians, point.longitudeRadians));
    }

    for (int i = 0; i < points.size(); ++i) {
        QVector<int> indices;
        indices.reserve(kNeighborCount);
        QVector<double> bestDistances;
        bestDistances.reserve(kNeighborCount);
        const QVector3D origin = positions.at(i);
        for (int j = 0; j < points.size(); ++j) {
            if (i == j) {
                continue;
            }
            const QVector3D delta = positions.at(j) - origin;
            const double distanceSquared = QVector3D::dotProduct(delta, delta);
            if (indices.size() < kNeighborCount) {
                indices.push_back(j);
                bestDistances.push_back(distanceSquared);
                continue;
            }

            int worstIndex = 0;
            double worstDistance = bestDistances.at(0);
            for (int k = 1; k < bestDistances.size(); ++k) {
                if (bestDistances.at(k) > worstDistance) {
                    worstDistance = bestDistances.at(k);
                    worstIndex = k;
                }
            }
            if (distanceSquared < worstDistance) {
                indices[worstIndex] = j;
                bestDistances[worstIndex] = distanceSquared;
            }
        }
        neighbors[i] = indices;
    }

    return neighbors;
}

double normalizeLongitude(double longitudeRad) {
    double value = std::fmod(longitudeRad + M_PI, 2.0 * M_PI);
    if (value < 0.0) {
        value += 2.0 * M_PI;
    }
    return value - M_PI;
}

WindVector clampWind(const WindVector &wind) {
    const double speed = std::hypot(wind.eastMps, wind.northMps);
    if (speed <= kMaxWindSpeedMps || speed <= 0.0) {
        return wind;
    }
    const double scale = kMaxWindSpeedMps / speed;
    return {wind.eastMps * scale, wind.northMps * scale};
}

int findBestSourceIndex(const PlanetSurfaceGrid &grid,
                        const QVector<QVector<int>> &neighborIndices,
                        int pointIndex,
                        double sourceLat,
                        double sourceLon) {
    const int pointCount = grid.pointCount();
    const double sinSource = qSin(sourceLat);
    const double cosSource = qCos(sourceLat);
    const QVector<int> &neighbors =
        (pointIndex < neighborIndices.size()) ? neighborIndices.at(pointIndex) : QVector<int>();
    int bestIndex = pointIndex;
    double bestCos = -1.0;
    for (const int neighborIndex : neighbors) {
        if (neighborIndex < 0 || neighborIndex >= pointCount) {
            continue;
        }
        const SurfacePoint &candidate = grid.points().at(neighborIndex);
        const double cosDistance = sinSource * candidate.sinLatitude +
            cosSource * candidate.cosLatitude *
            qCos(sourceLon - candidate.longitudeRadians);
        if (cosDistance > bestCos) {
            bestCos = cosDistance;
            bestIndex = neighborIndex;
        }
    }
    if (pointIndex >= 0 && pointIndex < pointCount) {
        const SurfacePoint &candidate = grid.points().at(pointIndex);
        const double cosDistance = sinSource * candidate.sinLatitude +
            cosSource * candidate.cosLatitude *
            qCos(sourceLon - candidate.longitudeRadians);
        if (cosDistance > bestCos) {
            bestCos = cosDistance;
            bestIndex = pointIndex;
        }
    }
    return bestIndex;
}
} // namespace

void AtmosphericAdvectionSolver::advectLayerWinds(const PlanetSurfaceGrid &grid,
                                                  AtmosphericGrid3D &atmosphereGrid,
                                                  double dayLengthSeconds,
                                                  double dtSeconds,
                                                  int smoothingIterations) const {
    const int pointCount = grid.pointCount();
    if (pointCount <= 0 || dtSeconds <= 0.0) {
        return;
    }

    if (atmosphereGrid.columnCount() != pointCount || atmosphereGrid.layerCount() <= 0) {
        return;
    }

    const double radiusMeters = grid.radiusKm() * 1000.0;
    if (radiusMeters <= 0.0) {
        return;
    }

    ensureNeighbors(grid);

    const int layerCount = atmosphereGrid.layerCount();
    const QVector<AtmosphericColumn> &columns = atmosphereGrid.columns();
    WindFieldModel windModel;

    for (int layerIndex = 0; layerIndex < layerCount; ++layerIndex) {
        QVector<double> pressuresAtm;
        QVector<double> temperaturesK;
        pressuresAtm.reserve(pointCount);
        temperaturesK.reserve(pointCount);
        for (const auto &column : columns) {
            if (layerIndex >= column.layers().size()) {
                pressuresAtm.push_back(0.0);
                temperaturesK.push_back(0.0);
                continue;
            }
            const auto &layer = column.layers().at(layerIndex);
            pressuresAtm.push_back(qMax(0.0, layer.pressureAtm()));
            temperaturesK.push_back(qMax(0.0, layer.temperatureKelvin()));
        }

        QVector<WindVector> wind =
            windModel.buildField(grid,
                                 pressuresAtm,
                                 temperaturesK,
                                 qMax(0.0, dayLengthSeconds),
                                 1);
        if (wind.size() != pointCount) {
            continue;
        }

        QVector<WindVector> advected;
        advected.resize(pointCount);
        for (int i = 0; i < pointCount; ++i) {
            const SurfacePoint &point = grid.points().at(i);
            // Полулагранжева схема для поля скорости:
            // назад по ветру на шаг dt, затем выбираем ближайший тайл в сетке.
            const double cosLat = qMax(kMinCosLatitude, std::abs(point.cosLatitude));
            const double dLat = (wind.at(i).northMps * dtSeconds) / radiusMeters;
            const double dLon = (wind.at(i).eastMps * dtSeconds) / (radiusMeters * cosLat);
            const double sourceLat = qBound(-kHalfPi, point.latitudeRadians - dLat, kHalfPi);
            const double sourceLon = normalizeLongitude(point.longitudeRadians - dLon);

            const int bestIndex = findBestSourceIndex(grid,
                                                      neighborIndices_,
                                                      i,
                                                      sourceLat,
                                                      sourceLon);
            advected[i] = wind.at(bestIndex);
        }

        const int iterations = qBound(0, smoothingIterations, 3);
        if (iterations > 0 && !neighborIndices_.isEmpty()) {
            QVector<WindVector> smoothed = advected;
            for (int iter = 0; iter < iterations; ++iter) {
                for (int i = 0; i < pointCount; ++i) {
                    const QVector<int> &neighbors = neighborIndices_.at(i);
                    if (neighbors.isEmpty()) {
                        smoothed[i] = advected.at(i);
                        continue;
                    }
                    double sumEast = 0.0;
                    double sumNorth = 0.0;
                    int count = 0;
                    for (const int neighborIndex : neighbors) {
                        if (neighborIndex < 0 || neighborIndex >= pointCount) {
                            continue;
                        }
                        sumEast += advected.at(neighborIndex).eastMps;
                        sumNorth += advected.at(neighborIndex).northMps;
                        ++count;
                    }
                    if (count <= 0) {
                        smoothed[i] = advected.at(i);
                        continue;
                    }
                    const double avgEast = sumEast / static_cast<double>(count);
                    const double avgNorth = sumNorth / static_cast<double>(count);
                    // Ограничитель/сглаживание для устойчивости: тянем поле
                    // к среднему по соседям, сохраняя динамику переноса.
                    smoothed[i].eastMps =
                        advected.at(i).eastMps + kSmoothingFactor * (avgEast - advected.at(i).eastMps);
                    smoothed[i].northMps =
                        advected.at(i).northMps + kSmoothingFactor * (avgNorth - advected.at(i).northMps);
                    smoothed[i] = clampWind(smoothed[i]);
                }
                advected.swap(smoothed);
            }
        }

        for (int i = 0; i < pointCount; ++i) {
            auto &layers = atmosphereGrid.columns()[i].layers();
            if (layerIndex >= layers.size()) {
                continue;
            }
            const WindVector clamped = clampWind(advected.at(i));
            layers[layerIndex].setWindUMps(clamped.eastMps);
            layers[layerIndex].setWindVMps(clamped.northMps);
        }
    }
}

void AtmosphericAdvectionSolver::advectLayerMoisture(const PlanetSurfaceGrid &grid,
                                                     AtmosphericGrid3D &atmosphereGrid,
                                                     double dtSeconds) const {
    const int pointCount = grid.pointCount();
    if (pointCount <= 0 || dtSeconds <= 0.0) {
        return;
    }

    if (atmosphereGrid.columnCount() != pointCount || atmosphereGrid.layerCount() <= 0) {
        return;
    }

    const double radiusMeters = grid.radiusKm() * 1000.0;
    if (radiusMeters <= 0.0) {
        return;
    }

    ensureNeighbors(grid);

    const int layerCount = atmosphereGrid.layerCount();
    for (int layerIndex = 0; layerIndex < layerCount; ++layerIndex) {
        QVector<double> vaporAdvected;
        QVector<double> liquidAdvected;
        vaporAdvected.resize(pointCount);
        liquidAdvected.resize(pointCount);

        for (int i = 0; i < pointCount; ++i) {
            const auto &layers = atmosphereGrid.columns().at(i).layers();
            if (layerIndex >= layers.size()) {
                vaporAdvected[i] = 0.0;
                liquidAdvected[i] = 0.0;
                continue;
            }

            const SurfacePoint &point = grid.points().at(i);
            const AtmosphericLayerState &layer = layers.at(layerIndex);
            // Полулагранжева адвекция влаги: ищем точку источника, откуда приходит
            // воздух за шаг dt, и переносим состав (пар/капли) без явного диффузионного
            // сглаживания — это стабильнее для больших dt на сферической сетке.
            const double cosLat = qMax(kMinCosLatitude, std::abs(point.cosLatitude));
            const double dLat = (layer.windVMps() * dtSeconds) / radiusMeters;
            const double dLon = (layer.windUMps() * dtSeconds) / (radiusMeters * cosLat);
            const double sourceLat = qBound(-kHalfPi, point.latitudeRadians - dLat, kHalfPi);
            const double sourceLon = normalizeLongitude(point.longitudeRadians - dLon);

            const int bestIndex = findBestSourceIndex(grid,
                                                      neighborIndices_,
                                                      i,
                                                      sourceLat,
                                                      sourceLon);
            const auto &sourceLayers = atmosphereGrid.columns().at(bestIndex).layers();
            if (layerIndex >= sourceLayers.size()) {
                vaporAdvected[i] = 0.0;
                liquidAdvected[i] = 0.0;
                continue;
            }

            const AtmosphericLayerState &sourceLayer = sourceLayers.at(layerIndex);
            vaporAdvected[i] = sourceLayer.waterVaporKgPerM2();
            liquidAdvected[i] = sourceLayer.liquidWaterKgPerM2();
        }

        for (int i = 0; i < pointCount; ++i) {
            auto &layers = atmosphereGrid.columns()[i].layers();
            if (layerIndex >= layers.size()) {
                continue;
            }
            const double vapor = qMax(0.0, vaporAdvected.at(i));
            const double liquid = qMax(0.0, liquidAdvected.at(i));
            layers[layerIndex].setWaterVaporKgPerM2(vapor);
            layers[layerIndex].setLiquidWaterKgPerM2(liquid);
            const double maxVapor =
                EvaporationModel::saturationWaterVaporKgPerM2(
                    layers[layerIndex].temperatureKelvin(),
                    layers[layerIndex].thicknessMeters());
            const double relativeHumidity =
                (maxVapor > 0.0) ? qBound(0.0, vapor / maxVapor, 1.0) : 0.0;
            layers[layerIndex].setRelativeHumidity(relativeHumidity);
        }
    }
}

void AtmosphericAdvectionSolver::ensureNeighbors(const PlanetSurfaceGrid &grid) const {
    if (cachedGrid_ == &grid && neighborIndices_.size() == grid.pointCount()) {
        return;
    }

    neighborIndices_ = buildNeighborIndices(grid.points());
    cachedGrid_ = &grid;
}
