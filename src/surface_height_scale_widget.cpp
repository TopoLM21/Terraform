#include "surface_height_scale_widget.h"

#include <QPainter>
#include <QtMath>

#include "height_color_scale.h"

SurfaceHeightScaleWidget::SurfaceHeightScaleWidget(QWidget *parent)
    : QWidget(parent) {}

void SurfaceHeightScaleWidget::setHeightRange(double minKm, double maxKm) {
    minHeightKm_ = minKm;
    maxHeightKm_ = maxKm;
    hasRange_ = minKm <= maxKm;
    update();
}

void SurfaceHeightScaleWidget::clearRange() {
    hasRange_ = false;
    update();
}

void SurfaceHeightScaleWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.fillRect(rect(), palette().window());

    const QRectF barRect = rect().adjusted(4, 4, -4, -4);
    if (barRect.width() <= 0.0 || barRect.height() <= 0.0) {
        return;
    }

    painter.setPen(Qt::NoPen);
    if (!hasRange_ || qFuzzyCompare(minHeightKm_, maxHeightKm_)) {
        painter.setBrush(palette().mid());
        painter.drawRoundedRect(barRect, 3.0, 3.0);
        return;
    }

    QLinearGradient gradient(barRect.topLeft(), barRect.topRight());
    // Нулевая отметка высоты визуально соответствует береговой линии,
    // поэтому цветовая шкала построена симметрично относительно середины.
    for (const auto &stop : heightColorStops()) {
        gradient.setColorAt(stop.position, stop.color);
    }
    painter.setBrush(gradient);
    painter.drawRoundedRect(barRect, 3.0, 3.0);
}
