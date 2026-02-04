#pragma once

#include <QColor>
#include <QPointF>
#include <QVector3D>

#include "StarOcclusionGeometry.h"

class QPainter;

class StarBackgroundRenderer {
public:
    void setAngularDiameterDegrees(double angularDiameterDeg);
    void setLightDirection(const QVector3D &direction);

    double angularDiameterDegrees() const;
    const QVector3D &lightDirection() const;

    void draw(QPainter &painter,
              const QPointF &screenCenter,
              double sphereRadiusPx,
              float yawDeg,
              float pitchDeg,
              const QColor &starColor) const;

private:
    QVector3D rotatedDirection(float yawDeg, float pitchDeg) const;

    double angularDiameterDeg_ = 0.5;
    QVector3D lightDirection_ = QVector3D(0.0f, 0.0f, 1.0f);
    StarOcclusionGeometry occlusionGeometry_;
};
