#pragma once

#include <QWidget>

class SurfaceWindScaleWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfaceWindScaleWidget(QWidget *parent = nullptr);

    void setWindRange(double minMps, double maxMps);
    void clearRange();

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    double minWindSpeedMps_ = 0.0;
    double maxWindSpeedMps_ = 0.0;
    bool hasRange_ = false;
};
