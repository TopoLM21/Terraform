#include "segment_selector_widget.h"

#include <QtCore/QLocale>
#include <QtCore/QtMath>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QSizePolicy>
#include <QtWidgets/QStyle>
#include <QtWidgets/QToolButton>

namespace {
constexpr int kSegmentCount = 12;
}

SegmentSelectorWidget::SegmentSelectorWidget(QWidget *parent)
    : QWidget(parent), segmentButtonGroup_(new QButtonGroup(this)) {
    auto *layout = new QHBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    previousButton_ = new QToolButton(this);
    previousButton_->setIcon(style()->standardIcon(QStyle::SP_ArrowLeft));
    previousButton_->setToolTip(QStringLiteral("Предыдущий сегмент"));
    layout->addWidget(previousButton_);

    auto *segmentsContainer = new QWidget(this);
    auto *segmentsLayout = new QGridLayout(segmentsContainer);
    segmentsLayout->setContentsMargins(0, 0, 0, 0);
    segmentsLayout->setSpacing(2);

    for (int index = 0; index < kSegmentCount; ++index) {
        auto *button = new QToolButton(segmentsContainer);
        button->setCheckable(true);
        button->setProperty("segmentButton", true);
        button->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        button->setEnabled(false);
        segmentsLayout->addWidget(button, 0, index);
        segmentButtons_.push_back(button);
        segmentButtonGroup_->addButton(button, index);
    }

    // Используем checked-состояние для подсветки выбора, чтобы не вводить ручной стиль отрисовки.
    setStyleSheet(QStringLiteral(
        "QToolButton[segmentButton=\"true\"]:checked {"
        "  background-color: palette(highlight);"
        "  color: palette(highlighted-text);"
        "}"
    ));

    layout->addWidget(segmentsContainer, 1);

    nextButton_ = new QToolButton(this);
    nextButton_->setIcon(style()->standardIcon(QStyle::SP_ArrowRight));
    nextButton_->setToolTip(QStringLiteral("Следующий сегмент"));
    layout->addWidget(nextButton_);

    segmentButtonGroup_->setExclusive(true);

    connect(segmentButtonGroup_, QOverload<int>::of(&QButtonGroup::buttonClicked),
            this, [this](int index) { setCurrentIndex(index); });
    connect(previousButton_, &QToolButton::clicked, this, [this]() {
        if (segments_.isEmpty()) {
            return;
        }
        const int nextIndex = (currentIndex_ <= 0)
            ? segments_.size() - 1
            : currentIndex_ - 1;
        setCurrentIndex(nextIndex);
    });
    connect(nextButton_, &QToolButton::clicked, this, [this]() {
        if (segments_.isEmpty()) {
            return;
        }
        const int nextIndex = (currentIndex_ + 1) % segments_.size();
        setCurrentIndex(nextIndex);
    });

    updateSegmentButtons();
}

void SegmentSelectorWidget::setSegments(const QVector<OrbitSegment> &segments) {
    segments_ = segments;
    updateSegmentButtons();

    if (segments_.isEmpty()) {
        updateSelection(-1, true);
        return;
    }

    if (currentIndex_ < 0 || currentIndex_ >= segments_.size()) {
        updateSelection(0, true);
    }
}

void SegmentSelectorWidget::setCurrentIndex(int index) {
    if (index == currentIndex_) {
        return;
    }

    if (index < -1 || (!segments_.isEmpty() && index >= segments_.size())) {
        return;
    }

    updateSelection(index, true);
}

int SegmentSelectorWidget::currentIndex() const {
    return currentIndex_;
}

void SegmentSelectorWidget::updateSelection(int index, bool emitSignal) {
    currentIndex_ = index;

    for (int i = 0; i < segmentButtons_.size(); ++i) {
        segmentButtons_[i]->setChecked(i == currentIndex_);
    }

    if (emitSignal) {
        emit currentIndexChanged(currentIndex_);
    }
}

void SegmentSelectorWidget::updateSegmentButtons() {
    const bool hasSegments = !segments_.isEmpty();
    const int availableSegments = segments_.size();

    for (int i = 0; i < segmentButtons_.size(); ++i) {
        auto *button = segmentButtons_[i];
        const bool inRange = i < availableSegments;
        button->setEnabled(hasSegments && inRange);
        if (inRange) {
            const auto &segment = segments_.at(i);
            button->setText(QString::number(segment.index + 1));
            button->setToolTip(formatSegmentLabel(segment));
        } else {
            button->setText(QString());
            button->setToolTip(QString());
        }
    }

    previousButton_->setEnabled(hasSegments);
    nextButton_->setEnabled(hasSegments);
    setEnabled(hasSegments);
}

QString SegmentSelectorWidget::formatSegmentLabel(const OrbitSegment &segment) const {
    const double meanAnomalyDegrees = qRadiansToDegrees(segment.meanAnomalyRadians);
    return QStringLiteral("Сегмент %1 (M=%2°, r=%3 а.е.)")
        .arg(segment.index + 1)
        .arg(QLocale().toString(meanAnomalyDegrees, 'f', 0))
        .arg(QLocale().toString(segment.distanceAU, 'f', 3));
}
