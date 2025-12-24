#pragma once

#include "orbit_segment_calculator.h"

#include <QtCore/QVector>
#include <QtWidgets/QWidget>

class QButtonGroup;
class QToolButton;

class SegmentSelectorWidget : public QWidget {
    Q_OBJECT

public:
    explicit SegmentSelectorWidget(QWidget *parent = nullptr);

    void setSegments(const QVector<OrbitSegment> &segments);
    void setCurrentIndex(int index);
    int currentIndex() const;

signals:
    void currentIndexChanged(int index);

private:
    void updateSelection(int index, bool emitSignal);
    void updateSegmentButtons();
    QString formatSegmentLabel(const OrbitSegment &segment) const;

    QVector<OrbitSegment> segments_;
    QVector<QToolButton *> segmentButtons_;
    QButtonGroup *segmentButtonGroup_ = nullptr;
    QToolButton *previousButton_ = nullptr;
    QToolButton *nextButton_ = nullptr;
    int currentIndex_ = -1;
};
