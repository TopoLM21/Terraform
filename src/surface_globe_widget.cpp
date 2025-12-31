#include "surface_globe_widget.h"

#include <algorithm>

#include <QPainter>
#include <QMouseEvent>
#include <QMatrix4x4>
#include <QPainterPath>
#include <QVector3D>
#include <QtMath>

#include "height_color_scale.h"
#include "temperature_color_scale.h"

namespace {
struct GlobePoint {
    QPointF position;
    double z = 0.0;
    QColor color;
    int pointIndex = -1;
};

struct GlobeCell {
    QPainterPath path;
    double depth = 0.0;
    QColor color;
};

struct ClippedSegment {
    QVector3D start;
    QVector3D end;
};

QVector3D latLonToCartesian(double latitudeDeg, double longitudeDeg) {
    const double latRad = qDegreesToRadians(latitudeDeg);
    const double lonRad = qDegreesToRadians(longitudeDeg);
    const double cosLat = qCos(latRad);
    return QVector3D(static_cast<float>(cosLat * qSin(lonRad)),
                     static_cast<float>(qSin(latRad)),
                     static_cast<float>(cosLat * qCos(lonRad)));
}

QVector<QVector3D> clipPolygonAgainstZ(const QVector<QVector3D> &input) {
    QVector<QVector3D> output;
    if (input.isEmpty()) {
        return output;
    }

    // Отсекаем по плоскости z=0 (видимая полусфера). Используем
    // Sutherland-Hodgman, чтобы сохранить форму многогранника после отсечения.
    QVector3D prev = input.last();
    bool prevInside = prev.z() > 0.0f;
    for (const auto &current : input) {
        const bool currentInside = current.z() > 0.0f;
        if (currentInside != prevInside) {
            const float t = (0.0f - prev.z()) / (current.z() - prev.z());
            const QVector3D intersection = prev + t * (current - prev);
            output.push_back(intersection);
        }
        if (currentInside) {
            output.push_back(current);
        }
        prev = current;
        prevInside = currentInside;
    }
    return output;
}

bool clipLineSegmentAgainstZ(const QVector3D &a, const QVector3D &b, ClippedSegment *segment) {
    const bool aInside = a.z() > 0.0f;
    const bool bInside = b.z() > 0.0f;
    if (!aInside && !bInside) {
        return false;
    }
    QVector3D start = a;
    QVector3D end = b;
    if (aInside != bInside) {
        const float t = (0.0f - a.z()) / (b.z() - a.z());
        const QVector3D intersection = a + t * (b - a);
        if (aInside) {
            end = intersection;
        } else {
            start = intersection;
        }
    }
    segment->start = start;
    segment->end = end;
    return true;
}

QVector<QVector3D> buildLongitudeLine(double longitudeDeg, int segments) {
    QVector<QVector3D> points;
    points.reserve(segments + 1);
    for (int i = 0; i <= segments; ++i) {
        const double latitude = -90.0 + 180.0 * (static_cast<double>(i) / segments);
        points.push_back(latLonToCartesian(latitude, longitudeDeg));
    }
    return points;
}

QVector<QVector3D> buildLatitudeLine(double latitudeDeg, int segments) {
    QVector<QVector3D> points;
    points.reserve(segments + 1);
    for (int i = 0; i <= segments; ++i) {
        const double longitude = -180.0 + 360.0 * (static_cast<double>(i) / segments);
        points.push_back(latLonToCartesian(latitudeDeg, longitude));
    }
    return points;
}
} // namespace

SurfaceGlobeWidget::SurfaceGlobeWidget(QWidget *parent)
    : QWidget(parent) {}

void SurfaceGlobeWidget::setGrid(const PlanetSurfaceGrid *grid) {
    grid_ = grid;
    if (grid_ && !grid_->points().isEmpty()) {
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
        minTemperatureK_ = minTemp;
        maxTemperatureK_ = maxTemp;
        minHeightKm_ = minHeight;
        maxHeightKm_ = maxHeight;
        minWindSpeedMps_ = minWind;
        maxWindSpeedMps_ = maxWind;
        minPressureAtm_ = minPressure;
        maxPressureAtm_ = maxPressure;
    }
    update();
}

void SurfaceGlobeWidget::setMapMode(SurfaceMapMode mode) {
    if (mapMode_ == mode) {
        return;
    }
    mapMode_ = mode;
    update();
}

void SurfaceGlobeWidget::setTemperatureRange(double minK, double maxK) {
    minTemperatureK_ = minK;
    maxTemperatureK_ = maxK;
    update();
}

void SurfaceGlobeWidget::setWindRange(double minMps, double maxMps) {
    minWindSpeedMps_ = minMps;
    maxWindSpeedMps_ = maxMps;
    update();
}

void SurfaceGlobeWidget::setPressureRange(double minAtm, double maxAtm) {
    minPressureAtm_ = minAtm;
    maxPressureAtm_ = maxAtm;
    update();
}

void SurfaceGlobeWidget::setMarkupVisible(bool visible) {
    if (markupVisible_ == visible) {
        return;
    }
    markupVisible_ = visible;
    update();
}

void SurfaceGlobeWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.fillRect(rect(), palette().window());

    projectedPoints_.clear();
    lastPointRadiusPx_ = 0.0;
    if (!grid_ || grid_->points().isEmpty()) {
        return;
    }

    const QPointF center(width() * 0.5, height() * 0.5);
    const double sphereRadius = 0.45 * qMin(width(), height());
    if (sphereRadius <= 1.0) {
        return;
    }

    QVector<GlobePoint> visiblePoints;
    visiblePoints.reserve(grid_->pointCount());
    QVector<GlobeCell> visibleCells;
    visibleCells.reserve(grid_->cells().size());

    for (int pointIndex = 0; pointIndex < grid_->points().size(); ++pointIndex) {
        const auto &point = grid_->points().at(pointIndex);
        const QVector3D normal = applyRotation(latLonToCartesian(point.latitudeDeg, point.longitudeDeg));
        if (normal.z() <= 0.0f) {
            continue;
        }

        const QPointF projected(center.x() + normal.x() * sphereRadius,
                                center.y() - normal.y() * sphereRadius);

        GlobePoint globePoint;
        globePoint.position = projected;
        globePoint.z = normal.z();
        globePoint.pointIndex = pointIndex;
        // Освещение отключено по требованию отображения без теней и подсветок.
        if (mapMode_ == SurfaceMapMode::Temperature) {
            globePoint.color = temperatureToColor(point.temperatureK);
        } else if (mapMode_ == SurfaceMapMode::Height) {
            globePoint.color = heightToColor(point.heightKm);
        } else if (mapMode_ == SurfaceMapMode::Pressure) {
            globePoint.color = pressureToColor(point.pressureAtm);
        } else {
            globePoint.color = windToColor(point.windSpeedMps);
        }
        visiblePoints.push_back(globePoint);
        projectedPoints_.push_back(ProjectedPoint{projected, pointIndex});
    }

    std::sort(visiblePoints.begin(), visiblePoints.end(), [](const GlobePoint &a, const GlobePoint &b) {
        return a.z < b.z;
    });

    if (!grid_->cells().isEmpty()) {
        for (const SurfaceCell &cell : grid_->cells()) {
            if (cell.polygon.size() < 3 || cell.pointIndex < 0 ||
                cell.pointIndex >= grid_->points().size()) {
                continue;
            }
            QVector<QVector3D> polygon3d;
            polygon3d.reserve(cell.polygon.size());
            for (const auto &vertex : cell.polygon) {
                const QVector3D rotated = applyRotation(latLonToCartesian(vertex.latitudeDeg,
                                                                          vertex.longitudeDeg));
                polygon3d.push_back(rotated);
            }

            QVector<QVector3D> clipped = clipPolygonAgainstZ(polygon3d);
            if (clipped.size() < 3) {
                continue;
            }

            QPainterPath path;
            const QPointF start(center.x() + clipped.first().x() * sphereRadius,
                                center.y() - clipped.first().y() * sphereRadius);
            path.moveTo(start);
            double depthSum = 0.0;
            for (int i = 0; i < clipped.size(); ++i) {
                const QVector3D &v = clipped[i];
                depthSum += v.z();
                const QPointF projected(center.x() + v.x() * sphereRadius,
                                         center.y() - v.y() * sphereRadius);
                if (i > 0) {
                    path.lineTo(projected);
                }
            }
            path.closeSubpath();

            GlobeCell cellDraw;
            cellDraw.path = path;
            cellDraw.depth = depthSum / static_cast<double>(clipped.size());
            const SurfacePoint &cellPoint = grid_->points().at(cell.pointIndex);
            if (mapMode_ == SurfaceMapMode::Temperature) {
                cellDraw.color = temperatureToColor(cellPoint.temperatureK);
            } else if (mapMode_ == SurfaceMapMode::Height) {
                cellDraw.color = heightToColor(cellPoint.heightKm);
            } else if (mapMode_ == SurfaceMapMode::Pressure) {
                cellDraw.color = pressureToColor(cellPoint.pressureAtm);
            } else {
                cellDraw.color = windToColor(cellPoint.windSpeedMps);
            }
            visibleCells.push_back(cellDraw);
        }

        std::sort(visibleCells.begin(), visibleCells.end(), [](const GlobeCell &a, const GlobeCell &b) {
            return a.depth < b.depth;
        });

        painter.setPen(Qt::NoPen);
        for (const auto &cell : visibleCells) {
            painter.setBrush(cell.color);
            painter.drawPath(cell.path);
        }
    } else {
        const double dotRadius = pointRadiusPx(visiblePoints.size(), sphereRadius);
        lastPointRadiusPx_ = dotRadius;
        painter.setPen(Qt::NoPen);

        for (const auto &point : visiblePoints) {
            painter.setBrush(point.color);
            painter.drawEllipse(point.position, dotRadius, dotRadius);
        }
    }

    const double dotRadius = pointRadiusPx(visiblePoints.size(), sphereRadius);
    lastPointRadiusPx_ = dotRadius;

    if (markupVisible_) {
        const QColor markupColor(255, 255, 255, 160);
        QPen markupPen(markupColor, qMax(1.0, sphereRadius * 0.004));
        markupPen.setCapStyle(Qt::RoundCap);
        painter.setPen(markupPen);
        painter.setBrush(Qt::NoBrush);

        const int lineSegments = 72;
        const QVector<QVector3D> equator = buildLatitudeLine(0.0, lineSegments);
        const QVector<QVector3D> meridianA = buildLongitudeLine(0.0, lineSegments);
        const QVector<QVector3D> meridianB = buildLongitudeLine(90.0, lineSegments);

        const auto drawLineStrip = [&](const QVector<QVector3D> &line) {
            if (line.size() < 2) {
                return;
            }
            for (int i = 1; i < line.size(); ++i) {
                const QVector3D start = applyRotation(line[i - 1]);
                const QVector3D end = applyRotation(line[i]);
                ClippedSegment segment;
                if (!clipLineSegmentAgainstZ(start, end, &segment)) {
                    continue;
                }
                const QPointF p1(center.x() + segment.start.x() * sphereRadius,
                                 center.y() - segment.start.y() * sphereRadius);
                const QPointF p2(center.x() + segment.end.x() * sphereRadius,
                                 center.y() - segment.end.y() * sphereRadius);
                painter.drawLine(p1, p2);
            }
        };

        drawLineStrip(equator);
        drawLineStrip(meridianB);

        // Линия через полюса выделена отдельно, чтобы отличаться от прочей разметки.
        QPen polarPen(QColor(180, 220, 255, 200), qMax(1.0, sphereRadius * 0.0045));
        polarPen.setCapStyle(Qt::RoundCap);
        painter.setPen(polarPen);
        drawLineStrip(meridianA);

        // Ось вращения рисуем как прямую через центр сферы.
        const QVector3D rotationAxis = applyRotation(QVector3D(0.0f, 1.0f, 0.0f));
        QPen axisPen(QColor(255, 220, 128, 200), qMax(1.0, sphereRadius * 0.005));
        axisPen.setCapStyle(Qt::RoundCap);
        painter.setPen(axisPen);
        const QVector3D axisStart = -rotationAxis;
        const QVector3D axisEnd = rotationAxis;
        ClippedSegment axisSegment;
        if (clipLineSegmentAgainstZ(axisStart, axisEnd, &axisSegment)) {
            const QPointF p1(center.x() + axisSegment.start.x() * sphereRadius,
                             center.y() - axisSegment.start.y() * sphereRadius);
            const QPointF p2(center.x() + axisSegment.end.x() * sphereRadius,
                             center.y() - axisSegment.end.y() * sphereRadius);
            painter.drawLine(p1, p2);
        }
    }
}

QColor SurfaceGlobeWidget::windToColor(double speedMps) const {
    if (qFuzzyCompare(minWindSpeedMps_, maxWindSpeedMps_)) {
        return temperatureColorForRatio(0.5);
    }
    const double t = qBound(0.0,
                            (speedMps - minWindSpeedMps_) /
                                (maxWindSpeedMps_ - minWindSpeedMps_),
                            1.0);
    return temperatureColorForRatio(t);
}

QColor SurfaceGlobeWidget::pressureToColor(double pressureAtm) const {
    if (qFuzzyCompare(minPressureAtm_, maxPressureAtm_)) {
        return temperatureColorForRatio(0.5);
    }
    const double t = qBound(0.0,
                            (pressureAtm - minPressureAtm_) /
                                (maxPressureAtm_ - minPressureAtm_),
                            1.0);
    return temperatureColorForRatio(t);
}

void SurfaceGlobeWidget::mousePressEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton && !projectedPoints_.isEmpty()) {
        const QPointF clickPos = event->pos();
        const double hitRadius = lastPointRadiusPx_;
        // Хит-тест делаем в 2D-проекции: точки рисуются как круги на экране,
        // поэтому используем расстояние до проекций и порог радиуса точки.
        double closestDistanceSq = hitRadius * hitRadius;
        int closestIndex = -1;
        for (const auto &point : projectedPoints_) {
            const double dx = point.position.x() - clickPos.x();
            const double dy = point.position.y() - clickPos.y();
            const double distanceSq = dx * dx + dy * dy;
            if (distanceSq <= closestDistanceSq) {
                closestDistanceSq = distanceSq;
                closestIndex = point.pointIndex;
            }
        }
        if (closestIndex >= 0) {
            emit pointClicked(closestIndex);
        }
    }
    if (event->button() == Qt::LeftButton) {
        isDragging_ = true;
        lastMousePos_ = event->pos();
    }
    QWidget::mousePressEvent(event);
}

void SurfaceGlobeWidget::mouseMoveEvent(QMouseEvent *event) {
    if (isDragging_ && (event->buttons() & Qt::LeftButton)) {
        const QPoint delta = event->pos() - lastMousePos_;
        lastMousePos_ = event->pos();

        const float sensitivity = 0.4f;
        yawDeg_ += delta.x() * sensitivity;
        pitchDeg_ += delta.y() * sensitivity;
        // Yaw - поворот вокруг вертикальной оси Y, pitch - вокруг оси X.
        // Pitch ограничен, чтобы избежать переворота «камеры» при взгляде через полюса.
        pitchDeg_ = qBound(-89.0f, pitchDeg_, 89.0f);
        update();
    }
    QWidget::mouseMoveEvent(event);
}

QColor SurfaceGlobeWidget::temperatureToColor(double temperatureK) const {
    if (qFuzzyCompare(minTemperatureK_, maxTemperatureK_)) {
        return temperatureColorForRatio(0.5);
    }
    const double t = qBound(0.0,
                            (temperatureK - minTemperatureK_) / (maxTemperatureK_ - minTemperatureK_),
                            1.0);
    return temperatureColorForRatio(t);
}

QColor SurfaceGlobeWidget::heightToColor(double heightKm) const {
    if (qFuzzyCompare(minHeightKm_, maxHeightKm_)) {
        return heightColorForRatio(0.5);
    }
    if (minHeightKm_ < 0.0 && maxHeightKm_ > 0.0) {
        // Разделяем океан и сушу на шкале: 0 км находится в центре цветовой палитры.
        if (heightKm <= 0.0) {
            const double t = qBound(0.0, heightKm / minHeightKm_, 1.0);
            return heightColorForRatio(0.5 - 0.5 * t);
        }
        const double t = qBound(0.0, heightKm / maxHeightKm_, 1.0);
        return heightColorForRatio(0.5 + 0.5 * t);
    }

    const double t = qBound(0.0,
                            (heightKm - minHeightKm_) / (maxHeightKm_ - minHeightKm_),
                            1.0);
    return heightColorForRatio(t);
}

double SurfaceGlobeWidget::pointRadiusPx(int pointCount, double sphereRadiusPx) const {
    if (pointCount <= 0 || sphereRadiusPx <= 0.0) {
        return 1.0;
    }
    const double sphereArea = M_PI * sphereRadiusPx * sphereRadiusPx;
    const double cellArea = sphereArea / static_cast<double>(pointCount);
    const double spacing = qSqrt(cellArea);
    return qBound(1.0, spacing * 0.45, 6.0);
}

QVector3D SurfaceGlobeWidget::applyRotation(const QVector3D &v) const {
    QMatrix4x4 rotation;
    rotation.rotate(yawDeg_, 0.0f, 1.0f, 0.0f);
    rotation.rotate(pitchDeg_, 1.0f, 0.0f, 0.0f);
    return rotation.mapVector(v);
}
