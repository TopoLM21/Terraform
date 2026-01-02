#pragma once

#include "planet_surface_grid.h"
#include "surface_map_mode.h"

#include <QWidget>
#include <QPointF>
#include <QVector>

class SurfaceGlobeWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfaceGlobeWidget(QWidget *parent = nullptr);

    void setGrid(const PlanetSurfaceGrid *grid);
    void setMapMode(SurfaceMapMode mode);
    void setTemperatureRange(double minK, double maxK);
    void setWindRange(double minMps, double maxMps);
    void setPressureRange(double minAtm, double maxAtm);
    void setMarkupVisible(bool visible);
    void setAxisTiltDegrees(double tiltDegrees);

signals:
    void pointClicked(int pointIndex);

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;

private:
    struct ProjectedPoint {
        QPointF position;
        int pointIndex = -1;
    };

    QColor temperatureToColor(double temperatureK) const;
    QColor heightToColor(double heightKm) const;
    QColor windToColor(double speedMps) const;
    QColor pressureToColor(double pressureAtm) const;
    double pointRadiusPx(int pointCount, double sphereRadiusPx) const;
    QVector3D applyRotation(const QVector3D &v) const;

    const PlanetSurfaceGrid *grid_ = nullptr;
    SurfaceMapMode mapMode_ = SurfaceMapMode::Temperature;
    double minTemperatureK_ = 200.0;
    double maxTemperatureK_ = 320.0;
    double minHeightKm_ = -5.0;
    double maxHeightKm_ = 5.0;
    double minWindSpeedMps_ = 0.0;
    double maxWindSpeedMps_ = 50.0;
    double minPressureAtm_ = 0.0;
    double maxPressureAtm_ = 2.0;
    float yawDeg_ = 0.0f;
    float pitchDeg_ = 0.0f;
    QPoint lastMousePos_;
    bool isDragging_ = false;
    QVector<ProjectedPoint> projectedPoints_;
    double lastPointRadiusPx_ = 0.0;
    bool markupVisible_ = false;
    double axisTiltDegrees_ = 0.0;
};
