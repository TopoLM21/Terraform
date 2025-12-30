#pragma once

#include <QWidget>

class SurfaceHeightScaleWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfaceHeightScaleWidget(QWidget *parent = nullptr);

    void setHeightRange(double minKm, double maxKm);
    void clearRange();

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    double minHeightKm_ = 0.0;
    double maxHeightKm_ = 0.0;
    bool hasRange_ = false;
};
