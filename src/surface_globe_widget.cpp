#include "surface_globe_widget.h"

#include <algorithm>

#include <QPainter>
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

    // Направление света выбрано слегка сверху-слева и с фронта, чтобы подчеркнуть объем
    // без резкого пересвета на контуре сферы.
    const QVector3D lightDir = QVector3D(-0.35f, 0.4f, 0.85f).normalized();

    QVector<GlobePoint> visiblePoints;
    visiblePoints.reserve(grid_->pointCount());

    for (const auto &point : grid_->points()) {
        const QVector3D normal = latLonToCartesian(point.latitudeDeg, point.longitudeDeg);
        if (normal.z() <= 0.0f) {
            continue;
        }

        const QPointF projected(center.x() + normal.x() * sphereRadius,
                                center.y() - normal.y() * sphereRadius);

        GlobePoint globePoint;
        globePoint.position = projected;
        globePoint.z = normal.z();
        globePoint.color =
            applyLambertLighting(temperatureToColor(point.temperatureK), normal, lightDir);
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

    const QColor limbColor = palette().windowText().color();
    QRadialGradient limbGradient(center, sphereRadius);
    limbGradient.setColorAt(0.0, QColor(0, 0, 0, 0));
    limbGradient.setColorAt(0.82, QColor(0, 0, 0, 0));
    limbGradient.setColorAt(1.0, QColor(limbColor.red(), limbColor.green(), limbColor.blue(), 70));
    painter.setBrush(limbGradient);
    painter.setPen(Qt::NoPen);
    painter.drawEllipse(center, sphereRadius, sphereRadius);
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

QColor SurfaceGlobeWidget::applyLambertLighting(const QColor &color,
                                                const QVector3D &normal,
                                                const QVector3D &lightDir) const {
    // Ламбертово освещение: I = max(0, dot(n, L)).
    const float intensity = qMax(0.0f, QVector3D::dotProduct(normal, lightDir));
    const float r = static_cast<float>(color.redF()) * intensity;
    const float g = static_cast<float>(color.greenF()) * intensity;
    const float b = static_cast<float>(color.blueF()) * intensity;
    return QColor::fromRgbF(r, g, b, color.alphaF());
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
