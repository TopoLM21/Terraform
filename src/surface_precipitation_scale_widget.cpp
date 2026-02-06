#include "surface_precipitation_scale_widget.h"

#include <QPainter>
#include <QtMath>

#include "precipitation_color_scale.h"

SurfacePrecipitationScaleWidget::SurfacePrecipitationScaleWidget(QWidget *parent)
    : QWidget(parent) {}

void SurfacePrecipitationScaleWidget::setPrecipitationRange(double minKgPerM2, double maxKgPerM2) {
    minPrecipitationKgPerM2_ = minKgPerM2;
    maxPrecipitationKgPerM2_ = maxKgPerM2;
    hasRange_ = minKgPerM2 <= maxKgPerM2;
    update();
}

void SurfacePrecipitationScaleWidget::clearRange() {
    hasRange_ = false;
    update();
}

void SurfacePrecipitationScaleWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.fillRect(rect(), palette().window());

    const QRectF barRect = rect().adjusted(4, 4, -4, -4);
    if (barRect.width() <= 0.0 || barRect.height() <= 0.0) {
        return;
    }

    painter.setPen(Qt::NoPen);
    if (!hasRange_ || qFuzzyCompare(minPrecipitationKgPerM2_, maxPrecipitationKgPerM2_)) {
        painter.setBrush(palette().mid());
        painter.drawRoundedRect(barRect, 3.0, 3.0);
        return;
    }

    QLinearGradient gradient(barRect.topLeft(), barRect.topRight());
    for (const auto &stop : precipitationColorStops()) {
        gradient.setColorAt(stop.position, stop.color);
    }
    painter.setBrush(gradient);
    painter.drawRoundedRect(barRect, 3.0, 3.0);
}
