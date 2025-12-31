#include "wind_field_model.h"

#include <QVector3D>
#include <QtMath>

#include <cmath>

namespace {
constexpr double kStandardPressurePa = 101325.0;
constexpr double kSpecificGasConstant = 287.05; // Дж/(кг·К), сухой воздух.
constexpr double kMinCoriolis = 1.0e-6; // 1/с, мягкий порог для экватора.
constexpr double kMaxWindSpeedMps = 150.0;
constexpr int kNeighborCount = 6;
constexpr double kViscosityFactor = 0.45;

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

QVector3D eastVector(double latitudeRad, double longitudeRad) {
    Q_UNUSED(latitudeRad)
    return QVector3D(static_cast<float>(qCos(longitudeRad)),
                     0.0f,
                     static_cast<float>(-qSin(longitudeRad)));
}

QVector3D northVector(double latitudeRad, double longitudeRad) {
    return QVector3D(static_cast<float>(-qSin(latitudeRad) * qSin(longitudeRad)),
                     static_cast<float>(qCos(latitudeRad)),
                     static_cast<float>(-qSin(latitudeRad) * qCos(longitudeRad)));
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

QVector<WindVector> WindFieldModel::buildField(const PlanetSurfaceGrid &grid,
                                               const QVector<double> &pressureAtm,
                                               const QVector<double> &temperatureK,
                                               double dayLengthSeconds,
                                               int smoothingIterations) const {
    const int pointCount = grid.pointCount();
    if (pointCount <= 0 || pressureAtm.size() != pointCount ||
        temperatureK.size() != pointCount) {
        return {};
    }

    const double radiusMeters = grid.radiusKm() * 1000.0;
    const double omega = (dayLengthSeconds > 0.0)
                             ? (2.0 * M_PI / dayLengthSeconds)
                             : 0.0;
    const QVector<QVector<int>> neighbors = buildNeighborIndices(grid.points());

    QVector<WindVector> wind;
    wind.resize(pointCount);
    QVector<QVector3D> positions;
    positions.reserve(pointCount);
    for (const auto &point : grid.points()) {
        positions.push_back(latLonToCartesian(point.latitudeRadians, point.longitudeRadians));
    }

    for (int i = 0; i < pointCount; ++i) {
        const SurfacePoint &point = grid.points().at(i);
        const double pressurePa = pressureAtm.at(i) * kStandardPressurePa;
        const double temperature = qMax(1.0, temperatureK.at(i));
        const double density = (pressurePa > 0.0)
                                   ? (pressurePa / (kSpecificGasConstant * temperature))
                                   : 0.0;
        if (density <= 0.0 || radiusMeters <= 0.0) {
            wind[i] = WindVector{};
            continue;
        }

        const QVector3D east = eastVector(point.latitudeRadians, point.longitudeRadians);
        const QVector3D north = northVector(point.latitudeRadians, point.longitudeRadians);
        const QVector3D origin = positions.at(i);
        const double basePressure = pressurePa;

        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;
        double spx = 0.0;
        double spy = 0.0;
        for (const int neighborIndex : neighbors.at(i)) {
            if (neighborIndex < 0 || neighborIndex >= pointCount) {
                continue;
            }
            const QVector3D delta = (positions.at(neighborIndex) - origin) * radiusMeters;
            const double dx = QVector3D::dotProduct(delta, east);
            const double dy = QVector3D::dotProduct(delta, north);
            const double dp = pressureAtm.at(neighborIndex) * kStandardPressurePa - basePressure;
            sxx += dx * dx;
            syy += dy * dy;
            sxy += dx * dy;
            spx += dp * dx;
            spy += dp * dy;
        }

        double gradX = 0.0;
        double gradY = 0.0;
        const double det = sxx * syy - sxy * sxy;
        if (std::abs(det) > 1.0e-8) {
            gradX = (spx * syy - spy * sxy) / det;
            gradY = (spy * sxx - spx * sxy) / det;
        }

        const double coriolis = 2.0 * omega * qSin(point.latitudeRadians);
        const double fAbs = std::abs(coriolis);
        const double fSafe = (coriolis >= 0.0 ? 1.0 : -1.0) * qMax(fAbs, kMinCoriolis);
        // Геострофический баланс: v ≈ (1 / (ρ f)) * (k × ∇p).
        // Здесь k × ∇p даёт поворот градиента давления на 90° в плоскости горизонта.
        WindVector localWind{-gradY / (density * fSafe), gradX / (density * fSafe)};
        // Мягкое затухание у экватора при малом f, чтобы избежать бесконечных скоростей.
        const double equatorDamping = fAbs / (fAbs + kMinCoriolis);
        localWind.eastMps *= equatorDamping;
        localWind.northMps *= equatorDamping;
        wind[i] = clampWind(localWind);
    }

    const int iterations = qBound(0, smoothingIterations, 3);
    if (iterations == 0) {
        return wind;
    }

    QVector<WindVector> smoothed = wind;
    for (int iter = 0; iter < iterations; ++iter) {
        for (int i = 0; i < pointCount; ++i) {
            const auto &neighborList = neighbors.at(i);
            if (neighborList.isEmpty()) {
                smoothed[i] = wind.at(i);
                continue;
            }
            double sumEast = 0.0;
            double sumNorth = 0.0;
            int count = 0;
            for (const int neighborIndex : neighborList) {
                if (neighborIndex < 0 || neighborIndex >= pointCount) {
                    continue;
                }
                sumEast += wind.at(neighborIndex).eastMps;
                sumNorth += wind.at(neighborIndex).northMps;
                ++count;
            }
            if (count <= 0) {
                smoothed[i] = wind.at(i);
                continue;
            }
            const double avgEast = sumEast / static_cast<double>(count);
            const double avgNorth = sumNorth / static_cast<double>(count);
            // Лапласовское сглаживание (вихревая вязкость):
            // v_{n+1} = v_n + ν (v̄ - v_n), где v̄ — среднее по соседям.
            smoothed[i].eastMps = wind.at(i).eastMps + kViscosityFactor * (avgEast - wind.at(i).eastMps);
            smoothed[i].northMps = wind.at(i).northMps +
                                   kViscosityFactor * (avgNorth - wind.at(i).northMps);
            smoothed[i] = clampWind(smoothed[i]);
        }
        wind.swap(smoothed);
    }

    return wind;
}
