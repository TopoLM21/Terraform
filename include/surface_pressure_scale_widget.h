#pragma once

#include <QWidget>

class SurfacePressureScaleWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfacePressureScaleWidget(QWidget *parent = nullptr);

    void setPressureRange(double minAtm, double maxAtm);
    void clearRange();

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    double minPressureAtm_ = 0.0;
    double maxPressureAtm_ = 0.0;
    bool hasRange_ = false;
};
