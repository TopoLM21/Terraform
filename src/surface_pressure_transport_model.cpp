#include "surface_pressure_transport_model.h"

#include <QVector3D>
#include <QtMath>

#include <cmath>

namespace {
constexpr int kNeighborCount = 6;
constexpr double kHalfPi = 0.5 * M_PI;
constexpr double kSmoothingFactor = 0.3;
constexpr double kMinCosLatitude = 1.0e-5;

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
} // namespace

QVector<double> SurfacePressureTransportModel::advectPressure(const PlanetSurfaceGrid &grid,
                                                              const QVector<double> &pressureAtm,
                                                              const QVector<double> &windEastMps,
                                                              const QVector<double> &windNorthMps,
                                                              double dtSeconds,
                                                              int smoothingIterations,
                                                              double minPressureAtm) const {
    const int pointCount = grid.pointCount();
    if (pointCount <= 0 || pressureAtm.size() != pointCount ||
        windEastMps.size() != pointCount || windNorthMps.size() != pointCount) {
        return pressureAtm;
    }

    if (dtSeconds <= 0.0) {
        return pressureAtm;
    }

    const double radiusMeters = grid.radiusKm() * 1000.0;
    if (radiusMeters <= 0.0) {
        return pressureAtm;
    }

    ensureNeighbors(grid);

    QVector<double> advected;
    advected.resize(pointCount);

    for (int i = 0; i < pointCount; ++i) {
        const SurfacePoint &point = grid.points().at(i);
        // Переносим поверхностное давление на уровне рельефа.
        // Используется полулагранжева схема: назад по ветру на шаг dt.
        // dφ = v_n / R * dt, dλ = v_e / (R cos φ) * dt.
        const double cosLat = qMax(kMinCosLatitude, std::abs(point.cosLatitude));
        const double dLat = (windNorthMps.at(i) * dtSeconds) / radiusMeters;
        const double dLon = (windEastMps.at(i) * dtSeconds) / (radiusMeters * cosLat);
        const double sourceLat = qBound(-kHalfPi, point.latitudeRadians - dLat, kHalfPi);
        const double sourceLon = normalizeLongitude(point.longitudeRadians - dLon);

        const double sinSource = qSin(sourceLat);
        const double cosSource = qCos(sourceLat);
        const QVector<int> &neighbors = (i < neighborIndices_.size())
                                            ? neighborIndices_.at(i)
                                            : QVector<int>();
        int bestIndex = i;
        double bestCos = -1.0;
        for (const int neighborIndex : neighbors) {
            if (neighborIndex < 0 || neighborIndex >= pointCount) {
                continue;
            }
            const SurfacePoint &candidate = grid.points().at(neighborIndex);
            const double cosDistance = sinSource * candidate.sinLatitude +
                cosSource * candidate.cosLatitude * qCos(sourceLon - candidate.longitudeRadians);
            if (cosDistance > bestCos) {
                bestCos = cosDistance;
                bestIndex = neighborIndex;
            }
        }
        if (i >= 0 && i < pointCount) {
            const SurfacePoint &candidate = grid.points().at(i);
            const double cosDistance = sinSource * candidate.sinLatitude +
                cosSource * candidate.cosLatitude * qCos(sourceLon - candidate.longitudeRadians);
            if (cosDistance > bestCos) {
                bestCos = cosDistance;
                bestIndex = i;
            }
        }

        advected[i] = pressureAtm.at(bestIndex);
    }

    const int iterations = qBound(0, smoothingIterations, 3);
    if (iterations > 0 && !neighborIndices_.isEmpty()) {
        QVector<double> smoothed = advected;
        for (int iter = 0; iter < iterations; ++iter) {
            for (int i = 0; i < pointCount; ++i) {
                const QVector<int> &neighbors = neighborIndices_.at(i);
                if (neighbors.isEmpty()) {
                    smoothed[i] = advected.at(i);
                    continue;
                }
                double sum = 0.0;
                int count = 0;
                for (const int neighborIndex : neighbors) {
                    if (neighborIndex < 0 || neighborIndex >= pointCount) {
                        continue;
                    }
                    sum += advected.at(neighborIndex);
                    ++count;
                }
                if (count <= 0) {
                    smoothed[i] = advected.at(i);
                    continue;
                }
                const double average = sum / static_cast<double>(count);
                // Стабилизирующее сглаживание для избежания неустойчивых колебаний.
                smoothed[i] = advected.at(i) + kSmoothingFactor * (average - advected.at(i));
            }
            advected.swap(smoothed);
        }
    }

    const double floorPressure = qMax(0.0, minPressureAtm);
    for (double &value : advected) {
        value = qMax(value, floorPressure);
    }

    return advected;
}

void SurfacePressureTransportModel::ensureNeighbors(const PlanetSurfaceGrid &grid) const {
    if (cachedGrid_ == &grid && neighborIndices_.size() == grid.pointCount()) {
        return;
    }

    neighborIndices_ = buildNeighborIndices(grid.points());
    cachedGrid_ = &grid;
}
