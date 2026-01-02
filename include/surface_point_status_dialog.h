#pragma once

#include "planet_surface_point.h"

#include <QDialog>

class QLabel;

class SurfacePointStatusDialog : public QDialog {
    Q_OBJECT

public:
    explicit SurfacePointStatusDialog(QWidget *parent = nullptr);

    void setPoint(const SurfacePoint &point);
    void clearPoint();

private:
    QLabel *latitudeValueLabel_ = nullptr;
    QLabel *longitudeValueLabel_ = nullptr;
    QLabel *temperatureValueLabel_ = nullptr;
    QLabel *airTemperatureValueLabel_ = nullptr;
    QLabel *heightValueLabel_ = nullptr;
    QLabel *pressureValueLabel_ = nullptr;
    QLabel *windValueLabel_ = nullptr;
    QLabel *materialValueLabel_ = nullptr;
};
