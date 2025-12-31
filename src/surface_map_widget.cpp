#include "surface_map_widget.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QToolTip>
#include <QTransform>
#include <QtMath>

#include "height_color_scale.h"
#include "temperature_color_scale.h"

namespace {
const double kMaxX = 2.0 * qSqrt(2.0);
const double kMaxY = qSqrt(2.0);

QRgb encodeId(int id) {
    const int value = id + 1;
    const int r = (value >> 16) & 0xff;
    const int g = (value >> 8) & 0xff;
    const int b = value & 0xff;
    return qRgba(r, g, b, 0xff);
}

int decodeId(QRgb pixel) {
    const int r = qRed(pixel);
    const int g = qGreen(pixel);
    const int b = qBlue(pixel);
    const int value = (r << 16) | (g << 8) | b;
    return value - 1;
}

QPointF projectToPixel(const QPointF &projected, const QSize &imageSize) {
    const double normalizedX = (projected.x() + kMaxX) / (2.0 * kMaxX);
    const double normalizedY = (kMaxY - projected.y()) / (2.0 * kMaxY);
    return QPointF(normalizedX * (imageSize.width() - 1),
                   normalizedY * (imageSize.height() - 1));
}

QVector<QPainterPath> buildCellPaths(const SurfaceCell &cell,
                                     const MollweideProjection &projection,
                                     const QSize &imageSize) {
    QVector<QPainterPath> paths;
    if (cell.polygon.size() < 3 || imageSize.isEmpty()) {
        return paths;
    }

    QVector<QPointF> pixelPoints;
    pixelPoints.reserve(cell.polygon.size());

    double previousLon = cell.polygon.first().longitudeDeg;
    for (const auto &vertex : cell.polygon) {
        double lon = vertex.longitudeDeg;
        // Разворачиваем долготы, чтобы избежать скачка на ±180° и сохранить
        // непрерывный контур для проекции Mollweide.
        while (lon - previousLon > 180.0) {
            lon -= 360.0;
        }
        while (lon - previousLon < -180.0) {
            lon += 360.0;
        }
        previousLon = lon;
        const QPointF projected = projection.project(vertex.latitudeDeg, lon);
        pixelPoints.push_back(projectToPixel(projected, imageSize));
    }

    QPainterPath basePath;
    basePath.moveTo(pixelPoints.first());
    for (int i = 1; i < pixelPoints.size(); ++i) {
        basePath.lineTo(pixelPoints[i]);
    }
    basePath.closeSubpath();
    paths.push_back(basePath);

    const double imageWidth = imageSize.width();
    const QRectF bounds = basePath.boundingRect();
    // Если контур ушёл за пределы карты, рисуем копию с переносом на ширину карты.
    if (bounds.left() < 0.0) {
        QTransform shift;
        shift.translate(imageWidth, 0.0);
        paths.push_back(shift.map(basePath));
    }
    if (bounds.right() > imageWidth) {
        QTransform shift;
        shift.translate(-imageWidth, 0.0);
        paths.push_back(shift.map(basePath));
    }

    return paths;
}
} // namespace

SurfaceMapWidget::SurfaceMapWidget(QWidget *parent)
    : QWidget(parent) {
    setMouseTracking(true);
}

void SurfaceMapWidget::setGrid(const PlanetSurfaceGrid *grid) {
    grid_ = grid;
    geometryCache_.clear();
    rebuildImages();
}

void SurfaceMapWidget::setMapMode(SurfaceMapMode mode) {
    if (mapMode_ == mode) {
        return;
    }
    mapMode_ = mode;
    rebuildImages();
}

void SurfaceMapWidget::setTemperatureRange(double minK, double maxK) {
    minTemperatureK_ = minK;
    maxTemperatureK_ = maxK;
    rebuildImages();
}

void SurfaceMapWidget::setWindRange(double minMps, double maxMps) {
    minWindSpeedMps_ = minMps;
    maxWindSpeedMps_ = maxMps;
    rebuildImages();
}

void SurfaceMapWidget::setInterpolationEnabled(bool enabled) {
    if (interpolationEnabled_ == enabled) {
        return;
    }
    interpolationEnabled_ = enabled;
    rebuildImages();
}

void SurfaceMapWidget::setRenderScale(double scale) {
    const double clampedScale = qMax(0.1, scale);
    if (qFuzzyCompare(renderScale_ + 1.0, clampedScale + 1.0)) {
        return;
    }
    renderScale_ = clampedScale;
    geometryCache_.clear();
    rebuildImages();
}

void SurfaceMapWidget::setInterpolationNeighborCount(int neighborCount) {
    const int clampedNeighbors = qMax(1, neighborCount);
    if (neighborCount_ == clampedNeighbors) {
        return;
    }
    neighborCount_ = clampedNeighbors;
    geometryCache_.clear();
    rebuildImages();
}

void SurfaceMapWidget::setInterpolationPower(double power) {
    const double clampedPower = qMax(0.1, power);
    if (qFuzzyCompare(interpolationPower_ + 1.0, clampedPower + 1.0)) {
        return;
    }
    interpolationPower_ = clampedPower;
    geometryCache_.clear();
    rebuildImages();
}

void SurfaceMapWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.fillRect(rect(), palette().window());
    if (!colorImage_.isNull()) {
        // Масштабирование изображения позволяет уменьшить renderScale ради скорости,
        // сохранив прозрачный для UI компромисс между качеством и производительностью.
        painter.drawImage(rect(), colorImage_);
    }
}

void SurfaceMapWidget::mousePressEvent(QMouseEvent *event) {
    if (!grid_ || idImage_.isNull()) {
        return;
    }
    const int id = pointIdAt(event->pos());
    if (id < 0 || id >= grid_->pointCount()) {
        return;
    }
    const SurfacePoint *point = grid_->pointAt(id);
    if (!point) {
        return;
    }
    QToolTip::showText(event->globalPos(), formatPointTooltip(*point), this);
}

void SurfaceMapWidget::resizeEvent(QResizeEvent *event) {
    QWidget::resizeEvent(event);
    geometryCache_.clear();
    rebuildImages();
}

void SurfaceMapWidget::rebuildImages() {
    if (!grid_ || grid_->points().isEmpty() || width() <= 0 || height() <= 0) {
        colorImage_ = QImage();
        idImage_ = QImage();
        geometryCache_.clear();
        update();
        return;
    }

    const QSize scaledSize = scaledRenderSize();
    colorImage_ = QImage(scaledSize, QImage::Format_ARGB32);
    colorImage_.fill(Qt::transparent);
    idImage_ = QImage(size(), QImage::Format_ARGB32);
    idImage_.fill(Qt::transparent);

    if (grid_->pointCount() > 0) {
        double minTemp = grid_->points().first().temperatureK;
        double maxTemp = minTemp;
        double minHeight = grid_->points().first().heightKm;
        double maxHeight = minHeight;
        double minWind = grid_->points().first().windSpeedMps;
        double maxWind = minWind;
        for (const auto &point : grid_->points()) {
            minTemp = qMin(minTemp, point.temperatureK);
            maxTemp = qMax(maxTemp, point.temperatureK);
            minHeight = qMin(minHeight, point.heightKm);
            maxHeight = qMax(maxHeight, point.heightKm);
            minWind = qMin(minWind, point.windSpeedMps);
            maxWind = qMax(maxWind, point.windSpeedMps);
        }
        if (minTemp != maxTemp) {
            minTemperatureK_ = minTemp;
            maxTemperatureK_ = maxTemp;
        }
        minHeightKm_ = minHeight;
        maxHeightKm_ = maxHeight;
        minWindSpeedMps_ = minWind;
        maxWindSpeedMps_ = maxWind;
    }

    const double baseRadius = pointRadiusPx(grid_->pointCount(), size());
    const double baseScaledRadius = pointRadiusPx(grid_->pointCount(), scaledSize);
    // В режиме точек увеличиваем радиус чуть сильнее, чтобы уменьшить заметные швы.
    const double dotScale = interpolationEnabled_ ? 1.0 : 1.25;
    const double pointRadius = qMin(baseRadius * dotScale, 6.0);
    const double scaledPointRadius = qMin(baseScaledRadius * dotScale, 6.0);

    QPainter idPainter(&idImage_);
    idPainter.setRenderHint(QPainter::Antialiasing, true);
    idPainter.setPen(Qt::NoPen);

    if (interpolationEnabled_) {
        const int neighborCount = qMin(qMax(1, neighborCount_), grid_->pointCount());
        const SurfaceMapCacheKey cacheKey = currentCacheKey(neighborCount);
        if (!geometryCache_.isValidFor(cacheKey, grid_)) {
            geometryCache_.rebuild(grid_, &projection_, cacheKey);
        }

        const auto &weights = geometryCache_.pixelWeights();
        const auto &mask = geometryCache_.insideMask();
        const auto &points = grid_->points();

        const int imageWidth = colorImage_.width();
        for (int y = 0; y < colorImage_.height(); ++y) {
            auto *scanLine = reinterpret_cast<QRgb *>(colorImage_.scanLine(y));
            int pixelIndex = y * imageWidth;
            for (int x = 0; x < imageWidth; ++x, ++pixelIndex) {
                if (pixelIndex >= mask.size() || !mask[pixelIndex]) {
                    scanLine[x] = qRgba(0, 0, 0, 0);
                    continue;
                }
                const PixelWeights &pixelWeights = weights.value(pixelIndex);
                if (pixelWeights.indices.isEmpty()) {
                    scanLine[x] = qRgba(0, 0, 0, 0);
                    continue;
                }

                double value = 0.0;
                for (int i = 0; i < pixelWeights.indices.size(); ++i) {
                    const int pointIndex = pixelWeights.indices[i];
                    if (pointIndex < 0 || pointIndex >= points.size()) {
                        continue;
                    }
                    double sample = 0.0;
                    if (mapMode_ == SurfaceMapMode::Temperature) {
                        sample = points[pointIndex].temperatureK;
                    } else if (mapMode_ == SurfaceMapMode::Height) {
                        sample = points[pointIndex].heightKm;
                    } else {
                        sample = points[pointIndex].windSpeedMps;
                    }
                    value += pixelWeights.weights[i] * sample;
                }

                if (mapMode_ == SurfaceMapMode::Temperature) {
                    scanLine[x] = temperatureToColor(value);
                } else if (mapMode_ == SurfaceMapMode::Height) {
                    scanLine[x] = heightToColor(value);
                } else {
                    scanLine[x] = windToColor(value);
                }
            }
        }
    }

    QPainter colorPainter;
    if (!interpolationEnabled_) {
        colorPainter.begin(&colorImage_);
        colorPainter.setRenderHint(QPainter::Antialiasing, true);
        colorPainter.setPen(Qt::NoPen);
    }

    if (!grid_->cells().isEmpty()) {
        for (const SurfaceCell &cell : grid_->cells()) {
            const int pointIndex = cell.pointIndex;
            if (pointIndex < 0 || pointIndex >= grid_->points().size()) {
                continue;
            }
            const SurfacePoint &point = grid_->points().at(pointIndex);
            const auto idPaths = buildCellPaths(cell, projection_, size());
            if (!interpolationEnabled_) {
                const auto colorPaths = buildCellPaths(cell, projection_, scaledSize);
                QRgb color = 0;
                if (mapMode_ == SurfaceMapMode::Temperature) {
                    color = temperatureToColor(point.temperatureK);
                } else if (mapMode_ == SurfaceMapMode::Height) {
                    color = heightToColor(point.heightKm);
                } else {
                    color = windToColor(point.windSpeedMps);
                }
                colorPainter.setBrush(QColor::fromRgb(color));
                for (const auto &path : colorPaths) {
                    colorPainter.drawPath(path);
                }
            }

            idPainter.setBrush(QColor::fromRgb(encodeId(pointIndex)));
            for (const auto &path : idPaths) {
                idPainter.drawPath(path);
            }
        }
    } else {
        for (int i = 0; i < grid_->pointCount(); ++i) {
            const SurfacePoint &point = grid_->points()[i];
            const QPoint idPixel = mapPointToPixel(point.latitudeDeg, point.longitudeDeg, size());
            const QRectF idEllipse(idPixel.x() - pointRadius,
                                   idPixel.y() - pointRadius,
                                   pointRadius * 2.0,
                                   pointRadius * 2.0);
            if (!interpolationEnabled_) {
                const QPoint pixel =
                    mapPointToPixel(point.latitudeDeg, point.longitudeDeg, scaledSize);
                if (colorImage_.rect().contains(pixel)) {
                    const QRectF ellipseRect(pixel.x() - scaledPointRadius,
                                             pixel.y() - scaledPointRadius,
                                             scaledPointRadius * 2.0,
                                             scaledPointRadius * 2.0);
                    QRgb color = 0;
                    if (mapMode_ == SurfaceMapMode::Temperature) {
                        color = temperatureToColor(point.temperatureK);
                    } else if (mapMode_ == SurfaceMapMode::Height) {
                        color = heightToColor(point.heightKm);
                    } else {
                        color = windToColor(point.windSpeedMps);
                    }
                    colorPainter.setBrush(QColor::fromRgb(color));
                    colorPainter.drawEllipse(ellipseRect);
                }
            }

            idPainter.setBrush(QColor::fromRgb(encodeId(i)));
            idPainter.drawEllipse(idEllipse);
        }
    }

    update();
}

QPoint SurfaceMapWidget::mapPointToPixel(double latitudeDeg,
                                         double longitudeDeg,
                                         const QSize &imageSize) const {
    if (imageSize.isEmpty()) {
        return QPoint();
    }
    const QPointF projected = projection_.project(latitudeDeg, longitudeDeg);
    const double normalizedX = (projected.x() + kMaxX) / (2.0 * kMaxX);
    const double normalizedY = (kMaxY - projected.y()) / (2.0 * kMaxY);
    const int x = qBound(0,
                         static_cast<int>(normalizedX * (imageSize.width() - 1)),
                         imageSize.width() - 1);
    const int y = qBound(0,
                         static_cast<int>(normalizedY * (imageSize.height() - 1)),
                         imageSize.height() - 1);
    return QPoint(x, y);
}

QRgb SurfaceMapWidget::temperatureToColor(double temperatureK) const {
    if (qFuzzyCompare(minTemperatureK_, maxTemperatureK_)) {
        return temperatureColorForRatio(0.5).rgb();
    }
    const double t = qBound(0.0,
                            (temperatureK - minTemperatureK_) / (maxTemperatureK_ - minTemperatureK_),
                            1.0);
    return temperatureColorForRatio(t).rgb();
}

QRgb SurfaceMapWidget::heightToColor(double heightKm) const {
    if (qFuzzyCompare(minHeightKm_, maxHeightKm_)) {
        return heightColorForRatio(0.5).rgb();
    }
    if (minHeightKm_ < 0.0 && maxHeightKm_ > 0.0) {
        // Делим диапазон по уровню моря: отрицательные высоты идут в нижнюю половину шкалы,
        // положительные - в верхнюю, чтобы береговая линия всегда попадала в середину.
        if (heightKm <= 0.0) {
            const double t = qBound(0.0, heightKm / minHeightKm_, 1.0);
            return heightColorForRatio(0.5 - 0.5 * t).rgb();
        }
        const double t = qBound(0.0, heightKm / maxHeightKm_, 1.0);
        return heightColorForRatio(0.5 + 0.5 * t).rgb();
    }

    const double t = qBound(0.0,
                            (heightKm - minHeightKm_) / (maxHeightKm_ - minHeightKm_),
                            1.0);
    return heightColorForRatio(t).rgb();
}

QRgb SurfaceMapWidget::windToColor(double speedMps) const {
    if (qFuzzyCompare(minWindSpeedMps_, maxWindSpeedMps_)) {
        return temperatureColorForRatio(0.5).rgb();
    }
    const double t = qBound(0.0,
                            (speedMps - minWindSpeedMps_) /
                                (maxWindSpeedMps_ - minWindSpeedMps_),
                            1.0);
    return temperatureColorForRatio(t).rgb();
}

int SurfaceMapWidget::pointIdAt(const QPoint &pixel) const {
    if (!idImage_.rect().contains(pixel)) {
        return -1;
    }
    return decodeId(idImage_.pixel(pixel));
}

QString SurfaceMapWidget::formatPointTooltip(const SurfacePoint &point) const {
    return QStringLiteral("Широта: %1°\nДолгота: %2°\nТемпература: %3 K\nВысота: %4 км\n"
                          "Давление: %5 атм\nВетер: %6 м/с")
        .arg(point.latitudeDeg, 0, 'f', 2)
        .arg(point.longitudeDeg, 0, 'f', 2)
        .arg(point.temperatureK, 0, 'f', 2)
        .arg(point.heightKm, 0, 'f', 2)
        .arg(point.pressureAtm, 0, 'f', 3)
        .arg(point.windSpeedMps, 0, 'f', 2);
}

double SurfaceMapWidget::pointRadiusPx(int pointCount, const QSize &imageSize) const {
    if (pointCount <= 0 || imageSize.isEmpty()) {
        return 1.0;
    }
    const double widgetArea =
        static_cast<double>(imageSize.width()) * static_cast<double>(imageSize.height());
    const double cellArea = widgetArea / static_cast<double>(pointCount);
    const double spacing = qSqrt(cellArea);
    // Радиус основан на среднем расстоянии между точками: при больших сетках он
    // уменьшается, чтобы точки не слипались, а при малых ограничивается сверху,
    // сохраняя читаемость без превращения точек в одиночные пиксели.
    return qBound(1.0, spacing * 0.45, 6.0);
}

QSize SurfaceMapWidget::scaledRenderSize() const {
    const int width = qMax(1, qRound(size().width() * renderScale_));
    const int height = qMax(1, qRound(size().height() * renderScale_));
    return QSize(width, height);
}

SurfaceMapCacheKey SurfaceMapWidget::currentCacheKey(int neighborCount) const {
    SurfaceMapCacheKey key;
    key.widgetSize = size();
    key.renderScale = renderScale_;
    key.neighborCount = neighborCount;
    key.power = interpolationPower_;
    return key;
}
