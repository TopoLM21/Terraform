#pragma once

#include "planet_surface_point.h"

#include <QDialog>

class QLabel;

class SurfacePointStatusDialog : public QDialog {
    Q_OBJECT

public:
    explicit SurfacePointStatusDialog(QWidget *parent = nullptr);

    void setPoint(const SurfacePoint &point, double tileAreaKm2, double tileEdgeLengthKm);
    void clearPoint();

private:
    QLabel *latitudeValueLabel_ = nullptr;
    QLabel *longitudeValueLabel_ = nullptr;
    QLabel *temperatureValueLabel_ = nullptr;
    QLabel *airTemperatureValueLabel_ = nullptr;
    QLabel *heightValueLabel_ = nullptr;
    QLabel *pressureValueLabel_ = nullptr;
    QLabel *windValueLabel_ = nullptr;
    QLabel *insolationValueLabel_ = nullptr;
    QLabel *materialValueLabel_ = nullptr;
    QLabel *tileAreaValueLabel_ = nullptr;
    QLabel *tileEdgeLengthValueLabel_ = nullptr;
};
