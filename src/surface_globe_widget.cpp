#include "surface_globe_widget.h"

#include <algorithm>

#include <QPainter>
#include <QMouseEvent>
#include <QMatrix4x4>
#include <QVector3D>
#include <QtMath>

#include "temperature_color_scale.h"

namespace {
struct GlobePoint {
    QPointF position;
    double z = 0.0;
    QColor color;
};

QVector3D latLonToCartesian(double latitudeDeg, double longitudeDeg) {
    const double latRad = qDegreesToRadians(latitudeDeg);
    const double lonRad = qDegreesToRadians(longitudeDeg);
    const double cosLat = qCos(latRad);
    return QVector3D(static_cast<float>(cosLat * qSin(lonRad)),
                     static_cast<float>(qSin(latRad)),
                     static_cast<float>(cosLat * qCos(lonRad)));
}
} // namespace

SurfaceGlobeWidget::SurfaceGlobeWidget(QWidget *parent)
    : QWidget(parent) {}

void SurfaceGlobeWidget::setGrid(const PlanetSurfaceGrid *grid) {
    grid_ = grid;
    if (grid_ && !grid_->points().isEmpty()) {
        double minTemp = grid_->points().first().temperatureK;
        double maxTemp = minTemp;
        for (const auto &point : grid_->points()) {
            minTemp = qMin(minTemp, point.temperatureK);
            maxTemp = qMax(maxTemp, point.temperatureK);
        }
        minTemperatureK_ = minTemp;
        maxTemperatureK_ = maxTemp;
    }
    update();
}

void SurfaceGlobeWidget::setTemperatureRange(double minK, double maxK) {
    minTemperatureK_ = minK;
    maxTemperatureK_ = maxK;
    update();
}

void SurfaceGlobeWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.fillRect(rect(), palette().window());

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

    for (const auto &point : grid_->points()) {
        const QVector3D normal = applyRotation(latLonToCartesian(point.latitudeDeg, point.longitudeDeg));
        if (normal.z() <= 0.0f) {
            continue;
        }

        const QPointF projected(center.x() + normal.x() * sphereRadius,
                                center.y() - normal.y() * sphereRadius);

        GlobePoint globePoint;
        globePoint.position = projected;
        globePoint.z = normal.z();
        // Освещение отключено по требованию отображения без теней и подсветок.
        globePoint.color = temperatureToColor(point.temperatureK);
        visiblePoints.push_back(globePoint);
    }

    std::sort(visiblePoints.begin(), visiblePoints.end(), [](const GlobePoint &a, const GlobePoint &b) {
        return a.z < b.z;
    });

    const double dotRadius = pointRadiusPx(visiblePoints.size(), sphereRadius);
    painter.setPen(Qt::NoPen);

    for (const auto &point : visiblePoints) {
        painter.setBrush(point.color);
        painter.drawEllipse(point.position, dotRadius, dotRadius);
    }

}

void SurfaceGlobeWidget::mousePressEvent(QMouseEvent *event) {
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
