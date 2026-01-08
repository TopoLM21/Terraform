#include "surface_realistic_scale_widget.h"

#include <QPainter>

SurfaceRealisticScaleWidget::SurfaceRealisticScaleWidget(QWidget *parent)
    : QWidget(parent) {}

void SurfaceRealisticScaleWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event)
    QPainter painter(this);
    painter.fillRect(rect(), palette().window());

    const QRectF barRect = rect().adjusted(4, 4, -4, -4);
    if (barRect.width() <= 0.0 || barRect.height() <= 0.0) {
        return;
    }

    painter.setPen(Qt::NoPen);
    QLinearGradient gradient(barRect.topLeft(), barRect.topRight());
    // Шкала собрана из характерных цветов: глубокая вода -> берег -> равнины -> горы -> снег.
    gradient.setColorAt(0.0, QColor(8, 28, 70));
    gradient.setColorAt(0.25, QColor(40, 140, 200));
    gradient.setColorAt(0.45, QColor(72, 132, 80));
    gradient.setColorAt(0.7, QColor(190, 160, 110));
    gradient.setColorAt(0.85, QColor(230, 230, 230));
    gradient.setColorAt(1.0, QColor(240, 245, 250));
    painter.setBrush(gradient);
    painter.drawRoundedRect(barRect, 3.0, 3.0);
}
