#pragma once

#include "planet_surface_grid.h"

#include <QWidget>

class SurfaceGlobeWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfaceGlobeWidget(QWidget *parent = nullptr);

    void setGrid(const PlanetSurfaceGrid *grid);
    void setTemperatureRange(double minK, double maxK);

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    QColor temperatureToColor(double temperatureK) const;
    double pointRadiusPx(int pointCount, double sphereRadiusPx) const;

    const PlanetSurfaceGrid *grid_ = nullptr;
    double minTemperatureK_ = 200.0;
    double maxTemperatureK_ = 320.0;
};
