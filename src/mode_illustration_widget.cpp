#include "mode_illustration_widget.h"

#include <QtGui/QPainter>
#include <QtGui/QPaintEvent>
#include <QtGui/QFontMetrics>

namespace {
constexpr int kOuterMargin = 8;
constexpr int kLabelSpacing = 4;
constexpr int kColumnSpacing = 14;
}

ModeIllustrationWidget::ModeIllustrationWidget(QWidget *parent)
    : QWidget(parent) {
    setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    setMinimumHeight(90);
}

void ModeIllustrationWidget::setRotationMode(RotationMode mode) {
    if (rotationMode_ == mode) {
        return;
    }
    rotationMode_ = mode;
    update();
}

RotationMode ModeIllustrationWidget::rotationMode() const {
    return rotationMode_;
}

QSize ModeIllustrationWidget::sizeHint() const {
    return QSize(200, 100);
}

void ModeIllustrationWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event);

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);

    const QRect content = rect().adjusted(kOuterMargin, kOuterMargin,
                                          -kOuterMargin, -kOuterMargin);
    if (content.width() <= 0 || content.height() <= 0) {
        return;
    }

    const QFontMetrics metrics(font());
    const int labelHeight = metrics.height();
    const int availableHeight = content.height() - labelHeight - kLabelSpacing;
    if (availableHeight <= 0) {
        return;
    }

    const int circleDiameter = qMin(availableHeight,
                                    (content.width() - kColumnSpacing) / 2);
    if (circleDiameter <= 0) {
        return;
    }

    const int circlesTop = content.top() + (availableHeight - circleDiameter) / 2;
    const int leftX = content.left();
    const int rightX = content.left() + circleDiameter + kColumnSpacing;

    const QRectF starRect(leftX, circlesTop, circleDiameter, circleDiameter);
    const QRectF planetRect(rightX, circlesTop, circleDiameter, circleDiameter);

    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(255, 215, 0));
    painter.drawEllipse(starRect);

    painter.setBrush(QColor(100, 105, 110));
    painter.drawEllipse(planetRect);

    painter.save();
    // Подсветка показывает условно освещённую часть диска планеты:
    // при обычном вращении освещена левая половина (сторона к звезде),
    // а при приливной синхронизации — верхняя половина как «фиксированный» день.
    QRectF highlightRect = planetRect;
    if (rotationMode_ == RotationMode::Normal) {
        highlightRect.setWidth(planetRect.width() / 2.0);
    } else {
        highlightRect.setHeight(planetRect.height() / 2.0);
    }
    painter.setClipRect(highlightRect);
    painter.setBrush(QColor(180, 190, 200));
    painter.drawEllipse(planetRect);
    painter.restore();

    painter.setPen(QColor(40, 40, 40));
    painter.setBrush(Qt::NoBrush);
    painter.drawEllipse(starRect);
    painter.drawEllipse(planetRect);

    const int labelTop = circlesTop + circleDiameter + kLabelSpacing + metrics.ascent();
    painter.drawText(QPointF(starRect.center().x() - metrics.horizontalAdvance(QStringLiteral("Звезда")) / 2.0,
                             labelTop),
                     QStringLiteral("Звезда"));
    painter.drawText(QPointF(planetRect.center().x() - metrics.horizontalAdvance(QStringLiteral("Планета")) / 2.0,
                             labelTop),
                     QStringLiteral("Планета"));
}
