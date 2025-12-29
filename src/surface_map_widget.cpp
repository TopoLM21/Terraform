#include "surface_map_widget.h"

#include <QMouseEvent>
#include <QPainter>
#include <QToolTip>
#include <QtMath>

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
} // namespace

SurfaceMapWidget::SurfaceMapWidget(QWidget *parent)
    : QWidget(parent) {
    setMouseTracking(true);
}

void SurfaceMapWidget::setGrid(const PlanetSurfaceGrid *grid) {
    grid_ = grid;
    rebuildImages();
}

void SurfaceMapWidget::setTemperatureRange(double minK, double maxK) {
    minTemperatureK_ = minK;
    maxTemperatureK_ = maxK;
    rebuildImages();
}

void SurfaceMapWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.fillRect(rect(), palette().window());
    if (!colorImage_.isNull()) {
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
    rebuildImages();
}

void SurfaceMapWidget::rebuildImages() {
    if (!grid_ || grid_->points().isEmpty() || width() <= 0 || height() <= 0) {
        colorImage_ = QImage();
        idImage_ = QImage();
        update();
        return;
    }

    colorImage_ = QImage(size(), QImage::Format_ARGB32);
    colorImage_.fill(Qt::transparent);
    idImage_ = QImage(size(), QImage::Format_ARGB32);
    idImage_.fill(Qt::transparent);

    if (grid_->pointCount() > 0) {
        double minTemp = grid_->points().first().temperatureK;
        double maxTemp = minTemp;
        for (const auto &point : grid_->points()) {
            minTemp = qMin(minTemp, point.temperatureK);
            maxTemp = qMax(maxTemp, point.temperatureK);
        }
        if (minTemp != maxTemp) {
            minTemperatureK_ = minTemp;
            maxTemperatureK_ = maxTemp;
        }
    }

    pointRadiusPx_ = pointRadiusPx(grid_->pointCount());

    QPainter colorPainter(&colorImage_);
    colorPainter.setRenderHint(QPainter::Antialiasing, true);
    colorPainter.setPen(Qt::NoPen);

    QPainter idPainter(&idImage_);
    idPainter.setRenderHint(QPainter::Antialiasing, false);
    idPainter.setPen(Qt::NoPen);

    for (int i = 0; i < grid_->pointCount(); ++i) {
        const SurfacePoint &point = grid_->points()[i];
        const QPoint pixel = mapPointToPixel(point.latitudeDeg, point.longitudeDeg);
        if (!colorImage_.rect().contains(pixel)) {
            continue;
        }
        const QRectF ellipseRect(pixel.x() - pointRadiusPx_,
                                 pixel.y() - pointRadiusPx_,
                                 pointRadiusPx_ * 2.0,
                                 pointRadiusPx_ * 2.0);
        colorPainter.setBrush(QColor::fromRgb(temperatureToColor(point.temperatureK)));
        colorPainter.drawEllipse(ellipseRect);

        idPainter.setBrush(QColor::fromRgb(encodeId(i)));
        idPainter.drawEllipse(ellipseRect);
    }

    update();
}

QPoint SurfaceMapWidget::mapPointToPixel(double latitudeDeg, double longitudeDeg) const {
    const QPointF projected = projection_.project(latitudeDeg, longitudeDeg);
    const double normalizedX = (projected.x() + kMaxX) / (2.0 * kMaxX);
    const double normalizedY = (kMaxY - projected.y()) / (2.0 * kMaxY);
    const int x = qBound(0, static_cast<int>(normalizedX * (width() - 1)), width() - 1);
    const int y = qBound(0, static_cast<int>(normalizedY * (height() - 1)), height() - 1);
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

int SurfaceMapWidget::pointIdAt(const QPoint &pixel) const {
    if (!idImage_.rect().contains(pixel)) {
        return -1;
    }
    return decodeId(idImage_.pixel(pixel));
}

QString SurfaceMapWidget::formatPointTooltip(const SurfacePoint &point) const {
    return QStringLiteral("Широта: %1°\nДолгота: %2°\nТемпература: %3 K\nВысота: %4 км")
        .arg(point.latitudeDeg, 0, 'f', 2)
        .arg(point.longitudeDeg, 0, 'f', 2)
        .arg(point.temperatureK, 0, 'f', 2)
        .arg(point.heightKm, 0, 'f', 2);
}

double SurfaceMapWidget::pointRadiusPx(int pointCount) const {
    if (pointCount <= 0) {
        return 1.0;
    }
    const double widgetArea = static_cast<double>(width()) * static_cast<double>(height());
    const double cellArea = widgetArea / static_cast<double>(pointCount);
    const double spacing = qSqrt(cellArea);
    // Радиус основан на среднем расстоянии между точками: при больших сетках он
    // уменьшается, чтобы точки не слипались, а при малых ограничивается сверху,
    // сохраняя читаемость без превращения точек в одиночные пиксели.
    return qBound(1.0, spacing * 0.45, 6.0);
}
