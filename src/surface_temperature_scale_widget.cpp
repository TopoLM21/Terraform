#include "surface_temperature_scale_widget.h"

#include <QPainter>
#include <QtMath>

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
    // Градиент повторяет схему карты температур: от холодного синего к тёплому красному.
    gradient.setColorAt(0.0, colorForRatio(0.0));
    gradient.setColorAt(0.5, colorForRatio(0.5));
    gradient.setColorAt(1.0, colorForRatio(1.0));
    painter.setBrush(gradient);
    painter.drawRoundedRect(barRect, 3.0, 3.0);
}

QColor SurfaceTemperatureScaleWidget::colorForRatio(double ratio) const {
    const double t = qBound(0.0, ratio, 1.0);
    const int r = static_cast<int>(255.0 * t);
    const int g = static_cast<int>(80.0 * (1.0 - t));
    const int b = static_cast<int>(255.0 * (1.0 - t));
    return QColor(r, g, b);
}
