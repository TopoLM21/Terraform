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

QVector<PixelWeights> SurfaceMapInterpolator::buildPixelWeights(const QSize &imageSize,
                                                                int neighborCount,
                                                                double power,
                                                                QVector<quint8> *insideMask) const {
    QVector<PixelWeights> result;
    if (!grid_ || !projection_ || grid_->points().isEmpty() || imageSize.isEmpty()) {
        if (insideMask) {
            insideMask->clear();
        }
        return result;
    }

    const int width = imageSize.width();
    const int height = imageSize.height();
    const int pixelCount = width * height;
    result.resize(pixelCount);
    if (insideMask) {
        insideMask->fill(0, pixelCount);
    }

    const int maxNeighbors = qMax(1, neighborCount);
    // Предрасчёт IDW-весов ускоряет последующие кадры: геометрия (проекция и
    // положение точек сетки) не меняется, поэтому можно повторно использовать
    // индексы соседей и веса при обновлении температур.
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const int pixelIndex = y * width + x;
            QPointF projected;
            const bool inside = pixelInsideProjection(QPoint(x, y), imageSize, &projected);
            if (insideMask) {
                (*insideMask)[pixelIndex] = inside ? 1 : 0;
            }
            if (!inside) {
                continue;
            }

            struct Neighbor {
                int index = -1;
                double distanceSquared = 0.0;
            };

            QVector<Neighbor> neighbors;
            neighbors.reserve(maxNeighbors);

            bool exactMatch = false;
            for (int i = 0; i < grid_->points().size(); ++i) {
                const auto &point = grid_->points()[i];
                const QPointF pointProjected = projection_->project(point.latitudeDeg,
                                                                    point.longitudeDeg);
                const double dx = pointProjected.x() - projected.x();
                const double dy = pointProjected.y() - projected.y();
                const double distanceSquared = dx * dx + dy * dy;
                if (distanceSquared < 1e-12) {
                    PixelWeights weights;
                    weights.indices.push_back(i);
                    weights.weights.push_back(1.0f);
                    result[pixelIndex] = std::move(weights);
                    exactMatch = true;
                    break;
                }

                if (neighbors.size() < maxNeighbors) {
                    neighbors.push_back({i, distanceSquared});
                    continue;
                }

                int worstIndex = 0;
                double worstDistance = neighbors[0].distanceSquared;
                for (int n = 1; n < neighbors.size(); ++n) {
                    if (neighbors[n].distanceSquared > worstDistance) {
                        worstDistance = neighbors[n].distanceSquared;
                        worstIndex = n;
                    }
                }

                if (distanceSquared < worstDistance) {
                    neighbors[worstIndex] = {i, distanceSquared};
                }
            }

            if (exactMatch || neighbors.isEmpty()) {
                continue;
            }

            PixelWeights weights;
            weights.indices.reserve(neighbors.size());
            weights.weights.reserve(neighbors.size());

            double weightSum = 0.0;
            // IDW: T(x) = sum(w_i * T_i) / sum(w_i), где w_i = 1 / d_i^p.
            // Степень p регулирует спад влияния соседей и позволяет балансировать
            // между локальными деталями и сглаживанием.
            for (const auto &neighbor : neighbors) {
                const double distance = qSqrt(neighbor.distanceSquared);
                if (distance <= 0.0) {
                    continue;
                }
                const double weight = 1.0 / qPow(distance, power);
                weights.indices.push_back(neighbor.index);
                weights.weights.push_back(static_cast<float>(weight));
                weightSum += weight;
            }

            if (qFuzzyIsNull(weightSum)) {
                continue;
            }

            for (auto &weight : weights.weights) {
                weight = static_cast<float>(weight / weightSum);
            }

            result[pixelIndex] = std::move(weights);
        }
    }

    return result;
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
