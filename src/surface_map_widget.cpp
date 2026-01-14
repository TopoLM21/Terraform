#include "surface_map_widget.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPalette>
#include <QPainterPath>
#include <QToolTip>
#include <QTransform>
#include <QVector3D>
#include <QtMath>

#include "height_color_scale.h"
#include "surface_realistic_color.h"
#include "temperature_color_scale.h"

namespace {
const double kMaxX = 2.0 * qSqrt(2.0);
const double kMaxY = qSqrt(2.0);
QVector3D latLonToUnitVector(double latitudeDeg, double longitudeDeg) {
    const double latRad = qDegreesToRadians(latitudeDeg);
    const double lonRad = qDegreesToRadians(longitudeDeg);
    const double cosLat = qCos(latRad);
    return QVector3D(static_cast<float>(cosLat * qCos(lonRad)),
                     static_cast<float>(qSin(latRad)),
                     static_cast<float>(cosLat * qSin(lonRad)));
}

double centralAngleRadians(const QVector3D &a, const QVector3D &b) {
    const QVector3D cross = QVector3D::crossProduct(a, b);
    const double crossLength = static_cast<double>(cross.length());
    const double dot = static_cast<double>(QVector3D::dotProduct(a, b));
    return qAtan2(crossLength, dot);
}

double sphericalTriangleAreaSteradians(const QVector3D &a,
                                       const QVector3D &b,
                                       const QVector3D &c) {
    // Формула сферического избытка для треугольника на единичной сфере:
    // E = 2 * atan2(|a · (b × c)|, 1 + a·b + b·c + c·a).
    const QVector3D cross = QVector3D::crossProduct(b, c);
    const double numerator =
        qAbs(static_cast<double>(QVector3D::dotProduct(a, cross)));
    const double denominator = 1.0 + static_cast<double>(QVector3D::dotProduct(a, b)) +
        static_cast<double>(QVector3D::dotProduct(b, c)) +
        static_cast<double>(QVector3D::dotProduct(c, a));
    return 2.0 * qAtan2(numerator, denominator);
}

double polygonAreaSteradians(const QVector<SurfaceVertex> &polygon) {
    if (polygon.size() < 3) {
        return 0.0;
    }
    const QVector3D origin = latLonToUnitVector(polygon[0].latitudeDeg,
                                                polygon[0].longitudeDeg);
    double area = 0.0;
    for (int i = 1; i + 1 < polygon.size(); ++i) {
        const QVector3D a = latLonToUnitVector(polygon[i].latitudeDeg,
                                               polygon[i].longitudeDeg);
        const QVector3D b = latLonToUnitVector(polygon[i + 1].latitudeDeg,
                                               polygon[i + 1].longitudeDeg);
        area += sphericalTriangleAreaSteradians(origin, a, b);
    }
    return area;
}

double polygonPerimeterRadians(const QVector<SurfaceVertex> &polygon) {
    if (polygon.size() < 2) {
        return 0.0;
    }
    double perimeter = 0.0;
    for (int i = 0; i < polygon.size(); ++i) {
        const int nextIndex = (i + 1) % polygon.size();
        const QVector3D a = latLonToUnitVector(polygon[i].latitudeDeg,
                                               polygon[i].longitudeDeg);
        const QVector3D b = latLonToUnitVector(polygon[nextIndex].latitudeDeg,
                                               polygon[nextIndex].longitudeDeg);
        perimeter += centralAngleRadians(a, b);
    }
    return perimeter;
}

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
    setAutoFillBackground(true);
    QPalette palette = this->palette();
    palette.setColor(QPalette::Window, Qt::black);
    setPalette(palette);
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

void SurfaceMapWidget::setPressureRange(double minAtm, double maxAtm) {
    minPressureAtm_ = minAtm;
    maxPressureAtm_ = maxAtm;
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
    QToolTip::showText(event->globalPos(), formatPointTooltip(*point, id), this);
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
        double minPressure = grid_->points().first().pressureAtm;
        double maxPressure = minPressure;
        for (const auto &point : grid_->points()) {
            minTemp = qMin(minTemp, point.temperatureK);
            maxTemp = qMax(maxTemp, point.temperatureK);
            minHeight = qMin(minHeight, point.heightKm);
            maxHeight = qMax(maxHeight, point.heightKm);
            minWind = qMin(minWind, point.windSpeedMps);
            maxWind = qMax(maxWind, point.windSpeedMps);
            minPressure = qMin(minPressure, point.pressureAtm);
            maxPressure = qMax(maxPressure, point.pressureAtm);
        }
        if (minTemp != maxTemp) {
            minTemperatureK_ = minTemp;
            maxTemperatureK_ = maxTemp;
        }
        minHeightKm_ = minHeight;
        maxHeightKm_ = maxHeight;
        minWindSpeedMps_ = minWind;
        maxWindSpeedMps_ = maxWind;
        minPressureAtm_ = minPressure;
        maxPressureAtm_ = maxPressure;
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

                if (mapMode_ == SurfaceMapMode::Realistic) {
                    double red = 0.0;
                    double green = 0.0;
                    double blue = 0.0;
                    for (int i = 0; i < pixelWeights.indices.size(); ++i) {
                        const int pointIndex = pixelWeights.indices[i];
                        if (pointIndex < 0 || pointIndex >= points.size()) {
                            continue;
                        }
                        const QColor color = realisticSurfaceColor(points[pointIndex],
                                                                   minHeightKm_,
                                                                   maxHeightKm_);
                        const double weight = pixelWeights.weights[i];
                        red += weight * color.red();
                        green += weight * color.green();
                        blue += weight * color.blue();
                    }
                    scanLine[x] = qRgba(qBound(0, static_cast<int>(qRound(red)), 255),
                                        qBound(0, static_cast<int>(qRound(green)), 255),
                                        qBound(0, static_cast<int>(qRound(blue)), 255),
                                        255);
                } else {
                    double value = 0.0;
                    for (int i = 0; i < pixelWeights.indices.size(); ++i) {
                        const int pointIndex = pixelWeights.indices[i];
                        if (pointIndex < 0 || pointIndex >= points.size()) {
                            continue;
                        }
                        double sample = 0.0;
                        if (mapMode_ == SurfaceMapMode::Temperature) {
                            sample = points[pointIndex].temperatureK;
                        } else if (mapMode_ == SurfaceMapMode::AirTemperature) {
                            sample = points[pointIndex].airTemperatureK;
                        } else if (mapMode_ == SurfaceMapMode::Height) {
                            sample = points[pointIndex].heightKm;
                        } else if (mapMode_ == SurfaceMapMode::Pressure) {
                            sample = points[pointIndex].pressureAtm;
                        } else {
                            sample = points[pointIndex].windSpeedMps;
                        }
                        value += pixelWeights.weights[i] * sample;
                    }

                    if (mapMode_ == SurfaceMapMode::Temperature ||
                        mapMode_ == SurfaceMapMode::AirTemperature) {
                        scanLine[x] = temperatureToColor(value);
                    } else if (mapMode_ == SurfaceMapMode::Height) {
                        scanLine[x] = heightToColor(value);
                    } else if (mapMode_ == SurfaceMapMode::Pressure) {
                        scanLine[x] = pressureToColor(value);
                    } else {
                        scanLine[x] = windToColor(value);
                    }
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
                } else if (mapMode_ == SurfaceMapMode::AirTemperature) {
                    color = temperatureToColor(point.airTemperatureK);
                } else if (mapMode_ == SurfaceMapMode::Height) {
                    color = heightToColor(point.heightKm);
                } else if (mapMode_ == SurfaceMapMode::Pressure) {
                    color = pressureToColor(point.pressureAtm);
                } else if (mapMode_ == SurfaceMapMode::Realistic) {
                    color = realisticSurfaceColor(point, minHeightKm_, maxHeightKm_).rgb();
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
                    } else if (mapMode_ == SurfaceMapMode::AirTemperature) {
                        color = temperatureToColor(point.airTemperatureK);
                    } else if (mapMode_ == SurfaceMapMode::Height) {
                        color = heightToColor(point.heightKm);
                    } else if (mapMode_ == SurfaceMapMode::Pressure) {
                        color = pressureToColor(point.pressureAtm);
                    } else if (mapMode_ == SurfaceMapMode::Realistic) {
                        color = realisticSurfaceColor(point, minHeightKm_, maxHeightKm_).rgb();
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

QRgb SurfaceMapWidget::pressureToColor(double pressureAtm) const {
    if (qFuzzyCompare(minPressureAtm_, maxPressureAtm_)) {
        return temperatureColorForRatio(0.5).rgb();
    }
    const double t = qBound(0.0,
                            (pressureAtm - minPressureAtm_) /
                                (maxPressureAtm_ - minPressureAtm_),
                            1.0);
    return temperatureColorForRatio(t).rgb();
}

int SurfaceMapWidget::pointIdAt(const QPoint &pixel) const {
    if (!idImage_.rect().contains(pixel)) {
        return -1;
    }
    return decodeId(idImage_.pixel(pixel));
}

QString SurfaceMapWidget::formatPointTooltip(const SurfacePoint &point, int pointIndex) const {
    const double areaKm2 = tileAreaKm2(pointIndex);
    const double edgeLengthKm = tileEdgeLengthKm(pointIndex);
    return QStringLiteral("Широта: %1°\nДолгота: %2°\nТемпература поверхности: %3 K\n"
                          "Температура воздуха: %4 K\nВысота: %5 км\nДавление: %6 атм\n"
                          "Солнечный поток: %7 Вт/м²\nВетер: %8 м/с\n"
                          "Площадь тайла: %9 км²\nСредняя длина ребра: %10 км")
        .arg(point.latitudeDeg, 0, 'f', 2)
        .arg(point.longitudeDeg, 0, 'f', 2)
        .arg(point.temperatureK, 0, 'f', 2)
        .arg(point.airTemperatureK, 0, 'f', 2)
        .arg(point.heightKm, 0, 'f', 2)
        .arg(point.pressureAtm, 0, 'f', 3)
        .arg(point.solarFluxWPerM2, 0, 'f', 2)
        .arg(point.windSpeedMps, 0, 'f', 2)
        .arg(areaKm2, 0, 'f', 2)
        .arg(edgeLengthKm, 0, 'f', 2);
}

double SurfaceMapWidget::tileAreaKm2(int pointIndex) const {
    if (!grid_) {
        return 0.0;
    }
    for (const auto &cell : grid_->cells()) {
        if (cell.pointIndex != pointIndex) {
            continue;
        }
        const double areaSteradians = polygonAreaSteradians(cell.polygon);
        return areaSteradians * grid_->radiusKm() * grid_->radiusKm();
    }
    return grid_->pointAreaKm2();
}

double SurfaceMapWidget::tileEdgeLengthKm(int pointIndex) const {
    if (!grid_) {
        return 0.0;
    }
    for (const auto &cell : grid_->cells()) {
        if (cell.pointIndex != pointIndex) {
            continue;
        }
        const double perimeterRadians = polygonPerimeterRadians(cell.polygon);
        if (cell.polygon.isEmpty()) {
            return 0.0;
        }
        return perimeterRadians * grid_->radiusKm() / cell.polygon.size();
    }
    const double areaKm2 = grid_->pointAreaKm2();
    if (areaKm2 <= 0.0) {
        return 0.0;
    }
    // Для равномерных точек без явного полигона используем эквивалентный
    // правильный шестиугольник, чтобы оценить длину ребра по площади.
    const double hexFactor = 3.0 * qSqrt(3.0) / 2.0;
    return qSqrt(areaKm2 / hexFactor);
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
