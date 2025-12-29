#include "surface_temperature_scale_widget.h"

#include <QPainter>
#include <QtMath>

#include "temperature_color_scale.h"

SurfaceTemperatureScaleWidget::SurfaceTemperatureScaleWidget(QWidget *parent)
    : QWidget(parent) {}

void SurfaceTemperatureScaleWidget::setTemperatureRange(double minK, double maxK) {
    minTemperatureK_ = minK;
    maxTemperatureK_ = maxK;
    hasRange_ = minK <= maxK;
    update();
}

void SurfaceTemperatureScaleWidget::clearRange() {
    hasRange_ = false;
    update();
}

void SurfaceTemperatureScaleWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.fillRect(rect(), palette().window());

    const QRectF barRect = rect().adjusted(4, 4, -4, -4);
    if (barRect.width() <= 0.0 || barRect.height() <= 0.0) {
        return;
    }

    painter.setPen(Qt::NoPen);
    if (!hasRange_ || qFuzzyCompare(minTemperatureK_, maxTemperatureK_)) {
        painter.setBrush(palette().mid());
        painter.drawRoundedRect(barRect, 3.0, 3.0);
        return;
    }

    QLinearGradient gradient(barRect.topLeft(), barRect.topRight());
    // Остановки размещены плотнее в холодной части диапазона, чтобы небольшие перепады
    // температуры читались лучше и не сливались в один оттенок.
    for (const auto &stop : temperatureColorStops()) {
        gradient.setColorAt(stop.position, stop.color);
    }
    painter.setBrush(gradient);
    painter.drawRoundedRect(barRect, 3.0, 3.0);
}
