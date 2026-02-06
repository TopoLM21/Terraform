#include "surface_globe_widget.h"

#include <algorithm>
#include <cmath>

#include <QPainter>
#include <QPalette>
#include <QMouseEvent>
#include <QMatrix4x4>
#include <QPainterPath>
#include <QVector3D>
#include <QtMath>

#include "height_color_scale.h"
#include "precipitation_color_scale.h"
#include "star_color.h"
#include "surface_realistic_color.h"
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

QColor mixCloudColor(const QColor &baseColor, double opacity) {
    const double t = qBound(0.0, opacity, 1.0);
    const QColor cloudColor(240, 240, 240);
    const auto blend = [t](int baseChannel, int cloudChannel) {
        return qBound(0,
                      static_cast<int>(qRound(baseChannel * (1.0 - t) + cloudChannel * t)),
                      255);
    };
    return QColor(blend(baseColor.red(), cloudColor.red()),
                  blend(baseColor.green(), cloudColor.green()),
                  blend(baseColor.blue(), cloudColor.blue()));
}

QVector3D latLonToCartesian(double latitudeDeg, double longitudeDeg) {
    const double latRad = qDegreesToRadians(latitudeDeg);
    const double lonRad = qDegreesToRadians(longitudeDeg);
    const double cosLat = qCos(latRad);
    return QVector3D(static_cast<float>(cosLat * qSin(lonRad)),
                     static_cast<float>(qSin(latRad)),
                     static_cast<float>(cosLat * qCos(lonRad)));
}

constexpr double kSolarRadiusMeters = 6.957e8;
constexpr double kAstronomicalUnitMeters = 1.496e11;

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
    : QWidget(parent) {
    setAutoFillBackground(true);
    QPalette palette = this->palette();
    palette.setColor(QPalette::Window, Qt::black);
    setPalette(palette);
    updateStarBackgroundColors();
}

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
        // Используем сглаженную интенсивность осадков (EMA), кг/м² за интервал.
        double minPrecipitation = grid_->points().first().precipitationKgPerM2;
        double maxPrecipitation = minPrecipitation;
        for (const auto &point : grid_->points()) {
            minTemp = qMin(minTemp, point.temperatureK);
            maxTemp = qMax(maxTemp, point.temperatureK);
            minHeight = qMin(minHeight, point.heightKm);
            maxHeight = qMax(maxHeight, point.heightKm);
            minWind = qMin(minWind, point.windSpeedMps);
            maxWind = qMax(maxWind, point.windSpeedMps);
            minPressure = qMin(minPressure, point.pressureAtm);
            maxPressure = qMax(maxPressure, point.pressureAtm);
            const double precipitationIntensity = point.precipitationKgPerM2;
            minPrecipitation = qMin(minPrecipitation, precipitationIntensity);
            maxPrecipitation = qMax(maxPrecipitation, precipitationIntensity);
        }
        minTemperatureK_ = minTemp;
        maxTemperatureK_ = maxTemp;
        minHeightKm_ = minHeight;
        maxHeightKm_ = maxHeight;
        minWindSpeedMps_ = minWind;
        maxWindSpeedMps_ = maxWind;
        minPressureAtm_ = minPressure;
        maxPressureAtm_ = maxPressure;
        minPrecipitationKgPerM2_ = minPrecipitation;
        maxPrecipitationKgPerM2_ = maxPrecipitation;
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

void SurfaceGlobeWidget::setPrecipitationRange(double minKgPerM2, double maxKgPerM2) {
    minPrecipitationKgPerM2_ = minKgPerM2;
    maxPrecipitationKgPerM2_ = maxKgPerM2;
    update();
}

void SurfaceGlobeWidget::setMarkupVisible(bool visible) {
    if (markupVisible_ == visible) {
        return;
    }
    markupVisible_ = visible;
    update();
}

void SurfaceGlobeWidget::setAxisTiltDegrees(double tiltDegrees) {
    if (qFuzzyCompare(axisTiltDegrees_, tiltDegrees)) {
        return;
    }
    axisTiltDegrees_ = tiltDegrees;
    update();
}

void SurfaceGlobeWidget::setStarDirection(const QVector3D &direction) {
    setStarLightDirection(direction);
}

void SurfaceGlobeWidget::setStarTemperature(double temperatureK) {
    if (temperatureK <= 0.0) {
        return;
    }
    setStarColor(starColorFromTemperature(temperatureK));
}

void SurfaceGlobeWidget::setStarDistanceAu(double distanceAu) {
    if (qFuzzyCompare(starDistanceAu_, distanceAu)) {
        return;
    }
    starDistanceAu_ = qMax(0.0, distanceAu);
    updateStarAngularDiameter();
}

void SurfaceGlobeWidget::setStarRadiusSolar(double radiusSolar) {
    if (qFuzzyCompare(starRadiusSolar_, radiusSolar)) {
        return;
    }
    starRadiusSolar_ = qMax(0.0, radiusSolar);
    updateStarAngularDiameter();
}

void SurfaceGlobeWidget::setStarLightDirection(const QVector3D &direction) {
    starBackgroundRenderer_.setLightDirection(direction);
    update();
}

void SurfaceGlobeWidget::setStarColor(const QColor &color) {
    starColor_ = color;
    if (starBackgroundAutoColorsEnabled_) {
        updateStarBackgroundColors();
    }
    update();
}

void SurfaceGlobeWidget::setStarAngularDiameterDegrees(double angularDiameterDeg) {
    if (qFuzzyCompare(starBackgroundAngularDiameterDeg_, angularDiameterDeg)) {
        return;
    }
    starBackgroundAngularDiameterDeg_ = qMax(0.0, angularDiameterDeg);
    starBackgroundAutoDiameterEnabled_ = false;
    starBackgroundRenderer_.setAngularDiameterDegrees(starBackgroundAngularDiameterDeg_);
    update();
}

void SurfaceGlobeWidget::setStarBackgroundAutoDiameterEnabled(bool enabled) {
    if (starBackgroundAutoDiameterEnabled_ == enabled) {
        return;
    }
    starBackgroundAutoDiameterEnabled_ = enabled;
    updateStarAngularDiameter();
}

void SurfaceGlobeWidget::setStarBackgroundAutoColorsEnabled(bool enabled) {
    if (starBackgroundAutoColorsEnabled_ == enabled) {
        return;
    }
    starBackgroundAutoColorsEnabled_ = enabled;
    updateStarBackgroundColors();
}

void SurfaceGlobeWidget::setStarBackgroundGradientColors(const QColor &coreColor,
                                                         const QColor &edgeColor) {
    if (starBackgroundCoreColor_ == coreColor && starBackgroundEdgeColor_ == edgeColor &&
        !starBackgroundAutoColorsEnabled_) {
        return;
    }
    starBackgroundCoreColor_ = coreColor;
    starBackgroundEdgeColor_ = edgeColor;
    starBackgroundAutoColorsEnabled_ = false;
    updateStarBackgroundColors();
}

void SurfaceGlobeWidget::setStarBackgroundIntensity(double intensity) {
    const double clamped = qMax(0.0, intensity);
    if (qFuzzyCompare(starBackgroundIntensity_ + 1.0, clamped + 1.0)) {
        return;
    }
    starBackgroundIntensity_ = clamped;
    starBackgroundRenderer_.setIntensity(starBackgroundIntensity_);
    update();
}

void SurfaceGlobeWidget::setCloudOpacityBoost(double boost) {
    const double clampedBoost = qBound(0.0, boost, 2.0);
    if (qFuzzyCompare(cloudOpacityBoost_ + 1.0, clampedBoost + 1.0)) {
        return;
    }
    cloudOpacityBoost_ = clampedBoost;
    update();
}

void SurfaceGlobeWidget::updateStarAngularDiameter() {
    if (!starBackgroundAutoDiameterEnabled_) {
        starBackgroundRenderer_.setAngularDiameterDegrees(starBackgroundAngularDiameterDeg_);
        update();
        return;
    }
    if (starRadiusSolar_ <= 0.0 || starDistanceAu_ <= 0.0) {
        starBackgroundRenderer_.setAngularDiameterDegrees(0.0);
        update();
        return;
    }

    const double radiusMeters = starRadiusSolar_ * kSolarRadiusMeters;
    const double distanceMeters = starDistanceAu_ * kAstronomicalUnitMeters;
    // Угловой диаметр: θ = 2 * atan(R / d).
    const double angularDiameterRad = 2.0 * std::atan(radiusMeters / distanceMeters);
    starBackgroundRenderer_.setAngularDiameterDegrees(qRadiansToDegrees(angularDiameterRad));
    update();
}

void SurfaceGlobeWidget::updateStarBackgroundColors() {
    if (starBackgroundAutoColorsEnabled_) {
        starBackgroundCoreColor_ = starColor_;
        starBackgroundEdgeColor_ = starColor_;
    }
    starBackgroundRenderer_.setGradientColors(starBackgroundCoreColor_, starBackgroundEdgeColor_);
    starBackgroundRenderer_.setIntensity(starBackgroundIntensity_);
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

    const QVector3D baseLightDir = starBackgroundRenderer_.lightDirection();
    // Световое направление хранится в мировых координатах, поэтому поворачиваем его
    // тем же поворотом камеры, что и геометрию, чтобы освещение и диск звезды
    // оставались синхронизированы при вращении.
    const QVector3D lightDir = applyRotation(baseLightDir).normalized();
    const double ambientFactor = 0.18;

    starBackgroundRenderer_.draw(painter, center, sphereRadius, yawDeg_, pitchDeg_);

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
        // Освещённость lambert: max(0, dot(normal, lightDir)), с мягкой ambient-добавкой.
        const double diffuse = qMax(0.0, static_cast<double>(QVector3D::dotProduct(normal, lightDir)));
        const double lightFactor = ambientFactor + (1.0 - ambientFactor) * diffuse;
        if (mapMode_ == SurfaceMapMode::Temperature) {
            globePoint.color = temperatureToColor(point.temperatureK);
        } else if (mapMode_ == SurfaceMapMode::AirTemperature) {
            globePoint.color = temperatureToColor(point.airTemperatureK);
        } else if (mapMode_ == SurfaceMapMode::Height) {
            globePoint.color = heightToColor(point.heightKm);
        } else if (mapMode_ == SurfaceMapMode::Pressure) {
            globePoint.color = pressureToColor(point.pressureAtm);
        } else if (mapMode_ == SurfaceMapMode::Precipitation) {
            globePoint.color = precipitationToColor(point.precipitationKgPerM2);
        } else if (mapMode_ == SurfaceMapMode::Realistic) {
            const QColor baseColor =
                realisticSurfaceColor(point, minHeightKm_, maxHeightKm_);
            const double cloudOpacity =
                qBound(0.0, point.cloudOpacity * cloudOpacityBoost_, 1.0);
            const QColor cloudyColor = mixCloudColor(baseColor, cloudOpacity);
            globePoint.color = applyLighting(cloudyColor, lightFactor);
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
            const QVector3D cellNormal =
                applyRotation(latLonToCartesian(cellPoint.latitudeDeg, cellPoint.longitudeDeg));
            const double diffuse =
                qMax(0.0, static_cast<double>(QVector3D::dotProduct(cellNormal, lightDir)));
            const double lightFactor = ambientFactor + (1.0 - ambientFactor) * diffuse;
            if (mapMode_ == SurfaceMapMode::Temperature) {
                cellDraw.color = temperatureToColor(cellPoint.temperatureK);
            } else if (mapMode_ == SurfaceMapMode::AirTemperature) {
                cellDraw.color = temperatureToColor(cellPoint.airTemperatureK);
            } else if (mapMode_ == SurfaceMapMode::Height) {
                cellDraw.color = heightToColor(cellPoint.heightKm);
            } else if (mapMode_ == SurfaceMapMode::Pressure) {
                cellDraw.color = pressureToColor(cellPoint.pressureAtm);
            } else if (mapMode_ == SurfaceMapMode::Precipitation) {
                cellDraw.color = precipitationToColor(cellPoint.precipitationKgPerM2);
            } else if (mapMode_ == SurfaceMapMode::Realistic) {
                const QColor baseColor =
                    realisticSurfaceColor(cellPoint, minHeightKm_, maxHeightKm_);
                const double cloudOpacity =
                    qBound(0.0, cellPoint.cloudOpacity * cloudOpacityBoost_, 1.0);
                const QColor cloudyColor = mixCloudColor(baseColor, cloudOpacity);
                cellDraw.color = applyLighting(cloudyColor, lightFactor);
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

        const int gridStepDeg = 30;
        for (int latitude = -90 + gridStepDeg; latitude < 90; latitude += gridStepDeg) {
            drawLineStrip(buildLatitudeLine(latitude, lineSegments));
        }
        for (int longitude = 0; longitude < 360; longitude += gridStepDeg) {
            drawLineStrip(buildLongitudeLine(longitude, lineSegments));
        }

        // Линия через полюса выделена отдельно, чтобы отличаться от прочей разметки.
        QPen polarPen(QColor(180, 220, 255, 200), qMax(1.0, sphereRadius * 0.0045));
        polarPen.setCapStyle(Qt::RoundCap);
        painter.setPen(polarPen);
        drawLineStrip(buildLongitudeLine(0.0, lineSegments));

        // Географическая ось (истинный полюс) — через широты ±90°.
        const QVector3D geographicAxis = applyRotation(latLonToCartesian(90.0, 0.0));
        QPen geographicAxisPen(QColor(120, 200, 255, 200), qMax(1.0, sphereRadius * 0.0045));
        geographicAxisPen.setCapStyle(Qt::RoundCap);
        painter.setPen(geographicAxisPen);
        const QVector3D geographicStart = -geographicAxis;
        const QVector3D geographicEnd = geographicAxis;
        ClippedSegment geographicSegment;
        if (clipLineSegmentAgainstZ(geographicStart, geographicEnd, &geographicSegment)) {
            const QPointF p1(center.x() + geographicSegment.start.x() * sphereRadius,
                             center.y() - geographicSegment.start.y() * sphereRadius);
            const QPointF p2(center.x() + geographicSegment.end.x() * sphereRadius,
                             center.y() - geographicSegment.end.y() * sphereRadius);
            painter.drawLine(p1, p2);
        }

        // Ось вращения: поворачиваем географическую ось на угол обликости вокруг оси X,
        // моделируя физический наклон оси вращения относительно истинного полюса.
        QMatrix4x4 obliquityRotation;
        obliquityRotation.rotate(static_cast<float>(axisTiltDegrees_), 1.0f, 0.0f, 0.0f);
        const QVector3D tiltedAxis =
            applyRotation(obliquityRotation.mapVector(latLonToCartesian(90.0, 0.0)));
        QPen rotationAxisPen(QColor(255, 180, 120, 220), qMax(1.0, sphereRadius * 0.0055));
        rotationAxisPen.setCapStyle(Qt::RoundCap);
        painter.setPen(rotationAxisPen);
        const QVector3D axisStart = -tiltedAxis;
        const QVector3D axisEnd = tiltedAxis;
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

QColor SurfaceGlobeWidget::precipitationToColor(double precipitationKgPerM2) const {
    if (qFuzzyCompare(minPrecipitationKgPerM2_, maxPrecipitationKgPerM2_)) {
        return precipitationColorForRatio(0.5);
    }
    const double t = qBound(0.0,
                            (precipitationKgPerM2 - minPrecipitationKgPerM2_) /
                                (maxPrecipitationKgPerM2_ - minPrecipitationKgPerM2_),
                            1.0);
    return precipitationColorForRatio(t);
}

QColor SurfaceGlobeWidget::applyLighting(const QColor &baseColor, double lightFactor) const {
    const double factor = qBound(0.0, lightFactor, 1.0);
    return QColor(static_cast<int>(baseColor.red() * factor),
                  static_cast<int>(baseColor.green() * factor),
                  static_cast<int>(baseColor.blue() * factor),
                  baseColor.alpha());
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
