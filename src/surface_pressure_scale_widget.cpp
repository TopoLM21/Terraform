#include "surface_pressure_scale_widget.h"

#include <QPainter>
#include <QtMath>

#include "temperature_color_scale.h"

SurfacePressureScaleWidget::SurfacePressureScaleWidget(QWidget *parent)
    : QWidget(parent) {}

void SurfacePressureScaleWidget::setPressureRange(double minAtm, double maxAtm) {
    minPressureAtm_ = minAtm;
    maxPressureAtm_ = maxAtm;
    hasRange_ = minAtm <= maxAtm;
    update();
}

void SurfacePressureScaleWidget::clearRange() {
    hasRange_ = false;
    update();
}

void SurfacePressureScaleWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.fillRect(rect(), palette().window());

    const QRectF barRect = rect().adjusted(4, 4, -4, -4);
    if (barRect.width() <= 0.0 || barRect.height() <= 0.0) {
        return;
    }

    painter.setPen(Qt::NoPen);
    if (!hasRange_ || qFuzzyCompare(minPressureAtm_, maxPressureAtm_)) {
        painter.setBrush(palette().mid());
        painter.drawRoundedRect(barRect, 3.0, 3.0);
        return;
    }

    QLinearGradient gradient(barRect.topLeft(), barRect.topRight());
    for (const auto &stop : temperatureColorStops()) {
        gradient.setColorAt(stop.position, stop.color);
    }
    painter.setBrush(gradient);
    painter.drawRoundedRect(barRect, 3.0, 3.0);
}
