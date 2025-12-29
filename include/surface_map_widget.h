#pragma once

#include "mollweide_projection.h"
#include "planet_surface_grid.h"
#include "surface_map_cache.h"

#include <QImage>
#include <QSize>
#include <QWidget>

class SurfaceMapWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfaceMapWidget(QWidget *parent = nullptr);

    void setGrid(const PlanetSurfaceGrid *grid);
    void setTemperatureRange(double minK, double maxK);
    void setInterpolationEnabled(bool enabled);
    void setRenderScale(double scale);
    void setInterpolationNeighborCount(int neighborCount);
    void setInterpolationPower(double power);

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;

private:
    void rebuildImages();
    QPoint mapPointToPixel(double latitudeDeg,
                           double longitudeDeg,
                           const QSize &imageSize) const;
    QRgb temperatureToColor(double temperatureK) const;
    int pointIdAt(const QPoint &pixel) const;
    QString formatPointTooltip(const SurfacePoint &point) const;
    double pointRadiusPx(int pointCount, const QSize &imageSize) const;
    QSize scaledRenderSize() const;
    SurfaceMapCacheKey currentCacheKey(int neighborCount) const;

    const PlanetSurfaceGrid *grid_ = nullptr;
    MollweideProjection projection_;
    QImage colorImage_;
    QImage idImage_;
    SurfaceMapCache geometryCache_;
    double minTemperatureK_ = 200.0;
    double maxTemperatureK_ = 320.0;
    bool interpolationEnabled_ = false;
    double renderScale_ = 1.0;
    int neighborCount_ = 8;
    double interpolationPower_ = 2.0;
};
