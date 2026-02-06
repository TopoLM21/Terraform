#pragma once

#include <QWidget>

class SurfacePrecipitationScaleWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfacePrecipitationScaleWidget(QWidget *parent = nullptr);

    void setPrecipitationRange(double minKgPerM2, double maxKgPerM2);
    void clearRange();

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    double minPrecipitationKgPerM2_ = 0.0;
    double maxPrecipitationKgPerM2_ = 0.0;
    bool hasRange_ = false;
};
