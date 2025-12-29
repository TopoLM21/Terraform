#pragma once

#include <QWidget>

class SurfaceTemperatureScaleWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfaceTemperatureScaleWidget(QWidget *parent = nullptr);

    void setTemperatureRange(double minK, double maxK);
    void clearRange();

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    double minTemperatureK_ = 0.0;
    double maxTemperatureK_ = 0.0;
    bool hasRange_ = false;
};
