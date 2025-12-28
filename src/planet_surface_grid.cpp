#include "planet_surface_grid.h"

#include <QtMath>

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kGoldenAngle = kPi * (3.0 - qSqrt(5.0));

double radiansToDegrees(double radians) {
    return radians * 180.0 / kPi;
}
} // namespace

void PlanetSurfaceGrid::setRadiusKm(double radiusKm) {
    radiusKm_ = radiusKm;
    if (points_.isEmpty()) {
        return;
    }
    const double surfaceAreaKm2 = 4.0 * kPi * radiusKm_ * radiusKm_;
    pointAreaKm2_ = surfaceAreaKm2 / static_cast<double>(points_.size());
    for (auto &point : points_) {
        point.radiusKm = radiusKm_;
    }
}

double PlanetSurfaceGrid::radiusKm() const {
    return radiusKm_;
}

void PlanetSurfaceGrid::generateFibonacciPoints(int pointCount) {
    points_.clear();
    if (pointCount <= 0 || radiusKm_ <= 0.0) {
        pointAreaKm2_ = 0.0;
        return;
    }

    points_.reserve(pointCount);
    for (int i = 0; i < pointCount; ++i) {
        // Фибоначчиева сфера: равномерное распределение по площади.
        // t задаёт равномерный шаг по sin(lat), а golden angle избегает слипания точек.
        const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(pointCount);
        const double latitude = qAsin(1.0 - 2.0 * t);
        const double longitude = qAtan2(qSin(kGoldenAngle * i), qCos(kGoldenAngle * i));

        SurfacePoint point;
        point.latitudeDeg = radiansToDegrees(latitude);
        point.longitudeDeg = radiansToDegrees(longitude);
        point.radiusKm = radiusKm_;
        points_.push_back(point);
    }

    const double surfaceAreaKm2 = 4.0 * kPi * radiusKm_ * radiusKm_;
    pointAreaKm2_ = surfaceAreaKm2 / static_cast<double>(points_.size());
}

int PlanetSurfaceGrid::pointCount() const {
    return points_.size();
}

double PlanetSurfaceGrid::pointAreaKm2() const {
    return pointAreaKm2_;
}

const QVector<SurfacePoint> &PlanetSurfaceGrid::points() const {
    return points_;
}

QVector<SurfacePoint> &PlanetSurfaceGrid::points() {
    return points_;
}

const SurfacePoint *PlanetSurfaceGrid::pointAt(int index) const {
    if (index < 0 || index >= points_.size()) {
        return nullptr;
    }
    return &points_[index];
}
