#include "surface_map_interpolator.h"

#include <QtMath>
namespace {
const double kMaxX = 2.0 * qSqrt(2.0);
const double kMaxY = qSqrt(2.0);

struct ProjectedNeighbor {
    int index = -1;
    double distanceSquared = 0.0;
};
} // namespace

struct SurfaceMapInterpolator::SpatialIndex {
    explicit SpatialIndex(const QVector<QPointF> &points) {
        if (points.isEmpty()) {
            return;
        }

        double minX = points.front().x();
        double maxX = minX;
        double minY = points.front().y();
        double maxY = minY;
        for (const auto &point : points) {
            minX = qMin(minX, point.x());
            maxX = qMax(maxX, point.x());
            minY = qMin(minY, point.y());
            maxY = qMax(maxY, point.y());
        }

        const double spanX = qMax(1e-6, maxX - minX);
        const double spanY = qMax(1e-6, maxY - minY);
        const int gridSide = qMax(1, static_cast<int>(qSqrt(points.size())));

        minX_ = minX;
        minY_ = minY;
        cellsX_ = qMax(1, gridSide);
        cellsY_ = qMax(1, gridSide);
        cellWidth_ = spanX / static_cast<double>(cellsX_);
        cellHeight_ = spanY / static_cast<double>(cellsY_);
        cells_.resize(cellsX_ * cellsY_);

        for (int i = 0; i < points.size(); ++i) {
            const int cellX = cellIndexX(points[i].x());
            const int cellY = cellIndexY(points[i].y());
            cells_[flatIndex(cellX, cellY)].push_back(i);
        }
    }

    void findNearest(const QPointF &target,
                     int maxNeighbors,
                     const QVector<QPointF> &points,
                     QVector<ProjectedNeighbor> *neighbors,
                     int *exactIndex) const {
        neighbors->clear();
        neighbors->reserve(maxNeighbors);
        if (cells_.isEmpty()) {
            return;
        }

        if (exactIndex) {
            *exactIndex = -1;
        }

        const int targetCellX = cellIndexX(target.x());
        const int targetCellY = cellIndexY(target.y());
        int minCellX = targetCellX;
        int maxCellX = targetCellX;
        int minCellY = targetCellY;
        int maxCellY = targetCellY;

        const int maxRadius = qMax(cellsX_, cellsY_);
        for (int radius = 0; radius <= maxRadius; ++radius) {
            const int startX = qMax(0, targetCellX - radius);
            const int endX = qMin(cellsX_ - 1, targetCellX + radius);
            const int startY = qMax(0, targetCellY - radius);
            const int endY = qMin(cellsY_ - 1, targetCellY + radius);

            auto scanCell = [&](int cellX, int cellY) -> bool {
                const auto &cell = cells_[flatIndex(cellX, cellY)];
                for (int index : cell) {
                    const QPointF point = points[index];
                    const double dx = point.x() - target.x();
                    const double dy = point.y() - target.y();
                    const double distanceSquared = dx * dx + dy * dy;
                    if (distanceSquared < 1e-12) {
                        if (exactIndex) {
                            *exactIndex = index;
                        }
                        return true;
                    }

                    if (neighbors->size() < maxNeighbors) {
                        neighbors->push_back({index, distanceSquared});
                        continue;
                    }

                    int worstIndex = 0;
                    double worstDistance = (*neighbors)[0].distanceSquared;
                    for (int i = 1; i < neighbors->size(); ++i) {
                        if ((*neighbors)[i].distanceSquared > worstDistance) {
                            worstDistance = (*neighbors)[i].distanceSquared;
                            worstIndex = i;
                        }
                    }

                    if (distanceSquared < worstDistance) {
                        (*neighbors)[worstIndex] = {index, distanceSquared};
                    }
                }
                return false;
            };

            for (int x = startX; x <= endX; ++x) {
                if (scanCell(x, startY)) {
                    return;
                }
                if (startY != endY && scanCell(x, endY)) {
                    return;
                }
            }

            for (int y = startY + 1; y <= endY - 1; ++y) {
                if (scanCell(startX, y)) {
                    return;
                }
                if (startX != endX && scanCell(endX, y)) {
                    return;
                }
            }

            minCellX = startX;
            maxCellX = endX;
            minCellY = startY;
            maxCellY = endY;

            if (neighbors->size() >= maxNeighbors) {
                double worstDistance = (*neighbors)[0].distanceSquared;
                for (int i = 1; i < neighbors->size(); ++i) {
                    worstDistance = qMax(worstDistance, (*neighbors)[i].distanceSquared);
                }

                const double minDistanceOutside = minDistanceToOutside(target,
                                                                       minCellX,
                                                                       maxCellX,
                                                                       minCellY,
                                                                       maxCellY);
                if (minDistanceOutside * minDistanceOutside >= worstDistance) {
                    return;
                }
            }

            if (startX == 0 && endX == cellsX_ - 1 && startY == 0 && endY == cellsY_ - 1) {
                return;
            }
        }
    }

private:
    int cellIndexX(double x) const {
        const double normalized = (x - minX_) / cellWidth_;
        return qBound(0, static_cast<int>(normalized), cellsX_ - 1);
    }

    int cellIndexY(double y) const {
        const double normalized = (y - minY_) / cellHeight_;
        return qBound(0, static_cast<int>(normalized), cellsY_ - 1);
    }

    int flatIndex(int cellX, int cellY) const {
        return cellY * cellsX_ + cellX;
    }

    double minDistanceToOutside(const QPointF &target,
                                int minCellX,
                                int maxCellX,
                                int minCellY,
                                int maxCellY) const {
        const double rectMinX = minX_ + static_cast<double>(minCellX) * cellWidth_;
        const double rectMaxX = minX_ + static_cast<double>(maxCellX + 1) * cellWidth_;
        const double rectMinY = minY_ + static_cast<double>(minCellY) * cellHeight_;
        const double rectMaxY = minY_ + static_cast<double>(maxCellY + 1) * cellHeight_;

        const double distanceLeft = target.x() - rectMinX;
        const double distanceRight = rectMaxX - target.x();
        const double distanceBottom = target.y() - rectMinY;
        const double distanceTop = rectMaxY - target.y();
        const double distanceToSide = qMin(qMin(distanceLeft, distanceRight),
                                           qMin(distanceBottom, distanceTop));
        return qMax(0.0, distanceToSide);
    }

    double minX_ = 0.0;
    double minY_ = 0.0;
    double cellWidth_ = 1.0;
    double cellHeight_ = 1.0;
    int cellsX_ = 0;
    int cellsY_ = 0;
    QVector<QVector<int>> cells_;
};

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
    QVector<QPointF> projectedPoints(grid_->points().size());
    // Проекция зависит только от геометрии сетки и одинакова для всех пикселей,
    // поэтому вычисляем её один раз, аналогично пояснению формулы IDW.
    for (int i = 0; i < grid_->points().size(); ++i) {
        const auto &point = grid_->points()[i];
        projectedPoints[i] = projection_->project(point.latitudeDeg, point.longitudeDeg);
    }

    SpatialIndex spatialIndex(projectedPoints);
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

            QVector<ProjectedNeighbor> neighbors;
            int exactIndex = -1;
            // Пространственный индекс ускоряет поиск кандидатов ближайших точек в
            // проекции, не меняя результат IDW по сравнению с полным перебором.
            spatialIndex.findNearest(projected,
                                     maxNeighbors,
                                     projectedPoints,
                                     &neighbors,
                                     &exactIndex);

            if (exactIndex >= 0) {
                PixelWeights weights;
                weights.indices.push_back(exactIndex);
                weights.weights.push_back(1.0f);
                result[pixelIndex] = std::move(weights);
                continue;
            }

            if (neighbors.isEmpty()) {
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
