#include "surface_map_interpolator.h"

#include <QtMath>
namespace {
const double kMaxX = 2.0 * qSqrt(2.0);
const double kMaxY = qSqrt(2.0);
} // namespace

SurfaceMapInterpolator::SurfaceMapInterpolator(const PlanetSurfaceGrid *grid,
                                               const MollweideProjection *projection)
    : grid_(grid),
      projection_(projection) {
}

double SurfaceMapInterpolator::interpolateTemperatureForPixel(const QPoint &pixel,
                                                               const QSize &imageSize,
                                                               int neighborCount,
                                                               double power,
                                                               bool *insideProjection) const {
    if (!grid_ || !projection_ || grid_->points().isEmpty()) {
        if (insideProjection) {
            *insideProjection = false;
        }
        return 0.0;
    }

    QPointF projected;
    const bool inside = pixelInsideProjection(pixel, imageSize, &projected);
    if (insideProjection) {
        *insideProjection = inside;
    }
    if (!inside) {
        return 0.0;
    }

    struct Neighbor {
        double distanceSquared = 0.0;
        double temperatureK = 0.0;
    };

    const int maxNeighbors = qMax(1, neighborCount);
    QVector<Neighbor> neighbors;
    neighbors.reserve(maxNeighbors);

    for (const auto &point : grid_->points()) {
        const QPointF pointProjected = projection_->project(point.latitudeDeg, point.longitudeDeg);
        const double dx = pointProjected.x() - projected.x();
        const double dy = pointProjected.y() - projected.y();
        const double distanceSquared = dx * dx + dy * dy;
        if (distanceSquared < 1e-12) {
            return point.temperatureK;
        }

        if (neighbors.size() < maxNeighbors) {
            neighbors.push_back({distanceSquared, point.temperatureK});
            continue;
        }

        int worstIndex = 0;
        double worstDistance = neighbors[0].distanceSquared;
        for (int i = 1; i < neighbors.size(); ++i) {
            if (neighbors[i].distanceSquared > worstDistance) {
                worstDistance = neighbors[i].distanceSquared;
                worstIndex = i;
            }
        }

        if (distanceSquared < worstDistance) {
            neighbors[worstIndex] = {distanceSquared, point.temperatureK};
        }
    }

    if (neighbors.isEmpty()) {
        return 0.0;
    }

    double weightedSum = 0.0;
    double weightSum = 0.0;
    // IDW: T(x) = sum(w_i * T_i) / sum(w_i), где w_i = 1 / d_i^p.
    // Степень p управляет спадом влияния: p = 2 даёт гладкую картину без
    // чрезмерного размазывания локальных аномалий.
    for (const auto &neighbor : neighbors) {
        const double distance = qSqrt(neighbor.distanceSquared);
        if (distance <= 0.0) {
            return neighbor.temperatureK;
        }
        const double weight = 1.0 / qPow(distance, power);
        weightedSum += weight * neighbor.temperatureK;
        weightSum += weight;
    }

    if (qFuzzyIsNull(weightSum)) {
        return 0.0;
    }

    return weightedSum / weightSum;
}

bool SurfaceMapInterpolator::pixelInsideProjection(const QPoint &pixel,
                                                   const QSize &imageSize,
                                                   QPointF *projected) const {
    if (imageSize.isEmpty()) {
        return false;
    }

    const QRect imageRect(QPoint(0, 0), imageSize);
    if (!imageRect.contains(pixel)) {
        return false;
    }

    const int width = imageSize.width();
    const int height = imageSize.height();
    const double denomX = qMax(1, width - 1);
    const double denomY = qMax(1, height - 1);
    const double normalizedX = static_cast<double>(pixel.x()) / denomX;
    const double normalizedY = static_cast<double>(pixel.y()) / denomY;

    const double x = normalizedX * (2.0 * kMaxX) - kMaxX;
    const double y = kMaxY - normalizedY * (2.0 * kMaxY);

    const double ellipseValue = (x * x) / (kMaxX * kMaxX) + (y * y) / (kMaxY * kMaxY);
    if (ellipseValue > 1.0) {
        return false;
    }

    if (projected) {
        *projected = QPointF(x, y);
    }
    return true;
}
