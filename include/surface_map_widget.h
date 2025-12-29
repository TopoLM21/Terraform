#pragma once

#include "mollweide_projection.h"
#include "planet_surface_grid.h"

#include <QImage>
#include <QWidget>

class SurfaceMapWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfaceMapWidget(QWidget *parent = nullptr);

    void setGrid(const PlanetSurfaceGrid *grid);
    void setTemperatureRange(double minK, double maxK);
    void setInterpolationEnabled(bool enabled);

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;

private:
    void rebuildImages();
    QPoint mapPointToPixel(double latitudeDeg, double longitudeDeg) const;
    QRgb temperatureToColor(double temperatureK) const;
    int pointIdAt(const QPoint &pixel) const;
    QString formatPointTooltip(const SurfacePoint &point) const;
    double pointRadiusPx(int pointCount) const;

    const PlanetSurfaceGrid *grid_ = nullptr;
    MollweideProjection projection_;
    QImage colorImage_;
    QImage idImage_;
    double minTemperatureK_ = 200.0;
    double maxTemperatureK_ = 320.0;
    double pointRadiusPx_ = 1.0;
    bool interpolationEnabled_ = false;
};
