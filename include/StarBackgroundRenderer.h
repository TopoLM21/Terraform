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
    void setGradientColors(const QColor &coreColor, const QColor &edgeColor);
    void setIntensity(double intensity);

    double angularDiameterDegrees() const;
    const QVector3D &lightDirection() const;
    QColor coreColor() const;
    QColor edgeColor() const;
    double intensity() const;

    void draw(QPainter &painter,
              const QPointF &screenCenter,
              double sphereRadiusPx,
              float yawDeg,
              float pitchDeg) const;

private:
    QVector3D rotatedDirection(float yawDeg, float pitchDeg) const;

    double angularDiameterDeg_ = 0.5;
    QVector3D lightDirection_ = QVector3D(0.0f, 0.0f, 1.0f);
    QColor coreColor_ = QColor(255, 244, 234);
    QColor edgeColor_ = QColor(255, 244, 234);
    double intensity_ = 1.0;
    StarOcclusionGeometry occlusionGeometry_;
};
