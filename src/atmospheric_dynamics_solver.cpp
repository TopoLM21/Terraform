#include "atmospheric_dynamics_solver.h"
#include "wind_field_model.h"

#include <QVector3D>
#include <QtConcurrent/QtConcurrentMap>
#include <QtMath>

#include <cmath>
#include <numeric>

namespace {
constexpr int kNeighborCount = 6;
constexpr double kStandardPressurePa = 101325.0;
constexpr double kSpecificGasConstant = 287.05; // Дж/(кг·К), сухой воздух.
constexpr double kMinCoriolis = 1.0e-6; // 1/с, мягкий порог для экватора.
constexpr double kViscosityFactor = 0.25;
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
} // namespace

void AtmosphericDynamicsSolver::updateLayerWinds(const PlanetSurfaceGrid &grid,
                                                 AtmosphericGrid3D &atmosphereGrid,
                                                 double dayLengthSeconds,
                                                 bool isRetrograde,
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

    const double omegaMagnitude = (dayLengthSeconds > 0.0)
                                      ? (2.0 * M_PI / dayLengthSeconds)
                                      : 0.0;
    const double rotationSign = isRetrograde ? -1.0 : 1.0;
    const double omega = rotationSign * omegaMagnitude;

    ensureNeighbors(grid);

    const int layerCount = atmosphereGrid.layerCount();
    const QVector<AtmosphericColumn> &columns = atmosphereGrid.columns();
    QVector<WindVector> nextWinds;
    QVector<WindVector> smoothed;

    QVector<int> pointIndices(pointCount);
    std::iota(pointIndices.begin(), pointIndices.end(), 0);

    for (int layerIndex = 0; layerIndex < layerCount; ++layerIndex) {
        nextWinds.resize(pointCount);

        // Вычисление градиентов давления и обновление ветра — параллельно по точкам.
        // Каждая точка читает только свои и соседние данные (read-only),
        // записывает в свой элемент nextWinds[i].
        QtConcurrent::blockingMap(pointIndices, [&](int i) {
            if (layerIndex >= columns.at(i).layers().size()) {
                nextWinds[i] = WindVector{};
                return;
            }

            const SurfacePoint &point = grid.points().at(i);
            const AtmosphericLayerState &layer = columns.at(i).layers().at(layerIndex);
            const double pressurePa = qMax(0.0, layer.pressureAtm()) * kStandardPressurePa;
            const double temperature = qMax(1.0, layer.temperatureKelvin());
            const double density = (pressurePa > 0.0)
                                       ? (pressurePa / (kSpecificGasConstant * temperature))
                                       : 0.0;
            if (density <= 0.0) {
                nextWinds[i] = WindVector{};
                return;
            }

            const QVector<int> &neighbors = (i < neighborIndices_.size())
                                                ? neighborIndices_.at(i)
                                                : QVector<int>();
            double gradLonPaPerM = 0.0;
            double gradLatPaPerM = 0.0;
            if (!neighbors.isEmpty()) {
                double sumLon = 0.0;
                double sumLat = 0.0;
                int count = 0;
                const double cosLat = qMax(1.0e-5, std::abs(point.cosLatitude));
                for (const int neighborIndex : neighbors) {
                    if (neighborIndex < 0 || neighborIndex >= pointCount) {
                        continue;
                    }
                    if (layerIndex >= columns.at(neighborIndex).layers().size()) {
                        continue;
                    }
                    const SurfacePoint &neighborPoint = grid.points().at(neighborIndex);
                    const AtmosphericLayerState &neighborLayer =
                        columns.at(neighborIndex).layers().at(layerIndex);
                    const double neighborPressurePa =
                        qMax(0.0, neighborLayer.pressureAtm()) * kStandardPressurePa;
                    const double dLon =
                        normalizeLongitude(neighborPoint.longitudeRadians - point.longitudeRadians);
                    const double dLat = neighborPoint.latitudeRadians - point.latitudeRadians;
                    const double distLon = radiusMeters * cosLat * dLon;
                    const double distLat = radiusMeters * dLat;
                    if (std::abs(distLon) > 1.0e-6) {
                        sumLon += (neighborPressurePa - pressurePa) / distLon;
                    }
                    if (std::abs(distLat) > 1.0e-6) {
                        sumLat += (neighborPressurePa - pressurePa) / distLat;
                    }
                    ++count;
                }
                if (count > 0) {
                    gradLonPaPerM = sumLon / static_cast<double>(count);
                    gradLatPaPerM = sumLat / static_cast<double>(count);
                }
            }

            // Уравнения движения (горизонтальная плоскость):
            // du/dt = -(1/ρ) ∂p/∂x + f v + ν ∇²u
            // dv/dt = -(1/ρ) ∂p/∂y - f u + ν ∇²v
            const double coriolis = 2.0 * omega * qSin(point.latitudeRadians);
            const double fAbs = std::abs(coriolis);
            const double fSafe = (coriolis >= 0.0 ? 1.0 : -1.0) * qMax(fAbs, kMinCoriolis);

            WindVector updated{layer.windUMps(), layer.windVMps()};
            const double accelU = -(gradLonPaPerM / density) + fSafe * updated.northMps;
            const double accelV = -(gradLatPaPerM / density) - fSafe * updated.eastMps;
            updated.eastMps += accelU * dtSeconds;
            updated.northMps += accelV * dtSeconds;
            nextWinds[i] = clampWind(updated);
        });

        // Сглаживание — параллельно по точкам.
        // Каждая точка читает nextWinds (read-only), записывает в smoothed[i].
        const int iterations = qBound(0, smoothingIterations, 3);
        if (iterations > 0 && !neighborIndices_.isEmpty()) {
            smoothed = nextWinds;
            // QVector implicit sharing: после копирования оба вектора делят данные.
            // Параллельная запись в smoothed[i] вызовет detach() из разных потоков —
            // гонка данных. Принудительно отделяем ДО параллельной секции.
            smoothed.detach();
            for (int iter = 0; iter < iterations; ++iter) {
                QtConcurrent::blockingMap(pointIndices, [&](int i) {
                    const QVector<int> &neighbors = neighborIndices_.at(i);
                    if (neighbors.isEmpty()) {
                        smoothed[i] = nextWinds.at(i);
                        return;
                    }
                    double sumEast = 0.0;
                    double sumNorth = 0.0;
                    int count = 0;
                    for (const int neighborIndex : neighbors) {
                        if (neighborIndex < 0 || neighborIndex >= pointCount) {
                            continue;
                        }
                        sumEast += nextWinds.at(neighborIndex).eastMps;
                        sumNorth += nextWinds.at(neighborIndex).northMps;
                        ++count;
                    }
                    if (count <= 0) {
                        smoothed[i] = nextWinds.at(i);
                        return;
                    }
                    const double avgEast = sumEast / static_cast<double>(count);
                    const double avgNorth = sumNorth / static_cast<double>(count);
                    // Вязкость/диффузия: подтягиваем ветер к среднему по соседям.
                    smoothed[i].eastMps =
                        nextWinds.at(i).eastMps + kViscosityFactor * (avgEast - nextWinds.at(i).eastMps);
                    smoothed[i].northMps =
                        nextWinds.at(i).northMps + kViscosityFactor * (avgNorth - nextWinds.at(i).northMps);
                    smoothed[i] = clampWind(smoothed[i]);
                });
                nextWinds.swap(smoothed);
            }
        }

        for (int i = 0; i < pointCount; ++i) {
            auto &layers = atmosphereGrid.columns()[i].layers();
            if (layerIndex >= layers.size()) {
                continue;
            }
            const WindVector clamped = clampWind(nextWinds.at(i));
            layers[layerIndex].setWindUMps(clamped.eastMps);
            layers[layerIndex].setWindVMps(clamped.northMps);
        }
    }
}

void AtmosphericDynamicsSolver::ensureNeighbors(const PlanetSurfaceGrid &grid) const {
    const double radiusKm = grid.radiusKm();
    if (cachedGrid_ == &grid && neighborIndices_.size() == grid.pointCount() &&
        qFuzzyCompare(radiusKm, cachedRadiusKm_)) {
        return;
    }

    neighborIndices_ = buildNeighborIndices(grid.points());
    cachedGrid_ = &grid;
    cachedRadiusKm_ = radiusKm;
}
