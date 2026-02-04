#include "StarBackgroundRenderer.h"

#include <QMatrix4x4>
#include <QPainter>
#include <QRadialGradient>
#include <QtMath>

void StarBackgroundRenderer::setAngularDiameterDegrees(double angularDiameterDeg) {
    angularDiameterDeg_ = qMax(0.0, angularDiameterDeg);
}

void StarBackgroundRenderer::setLightDirection(const QVector3D &direction) {
    if (direction.lengthSquared() < 1e-6f) {
        lightDirection_ = QVector3D(0.0f, 0.0f, 1.0f);
    } else {
        lightDirection_ = direction.normalized();
    }
}

double StarBackgroundRenderer::angularDiameterDegrees() const {
    return angularDiameterDeg_;
}

const QVector3D &StarBackgroundRenderer::lightDirection() const {
    return lightDirection_;
}

void StarBackgroundRenderer::draw(QPainter &painter,
                                  const QPointF &screenCenter,
                                  double sphereRadiusPx,
                                  float yawDeg,
                                  float pitchDeg,
                                  const QColor &starColor) const {
    if (angularDiameterDeg_ <= 0.0 || sphereRadiusPx <= 0.0) {
        return;
    }

    const QVector3D lightDir = rotatedDirection(yawDeg, pitchDeg).normalized();

    // Угловой радиус = половина диаметра; линейный размер на экране r ≈ R * tan(angularRadius).
    const double angularRadiusRad = qDegreesToRadians(angularDiameterDeg_ * 0.5);
    const double starRadius = sphereRadiusPx * qTan(angularRadiusRad);
    if (starRadius <= 0.0) {
        return;
    }

    const StarOcclusionGeometry::Result occlusion =
        occlusionGeometry_.evaluate(screenCenter, sphereRadiusPx, lightDir, starRadius);
    if (occlusion.hidden) {
        return;
    }

    QRadialGradient gradient(occlusion.starCenter, starRadius);
    QColor coreColor = starColor;
    coreColor.setAlpha(220);
    QColor edgeColor = starColor;
    edgeColor.setAlpha(0);
    gradient.setColorAt(0.0, coreColor);
    gradient.setColorAt(1.0, edgeColor);

    painter.setPen(Qt::NoPen);
    painter.setBrush(gradient);
    if (occlusion.hasClip) {
        painter.save();
        painter.setClipPath(occlusion.clipPath);
        painter.drawEllipse(occlusion.starCenter, starRadius, starRadius);
        painter.restore();
    } else {
        painter.drawEllipse(occlusion.starCenter, starRadius, starRadius);
    }
}

QVector3D StarBackgroundRenderer::rotatedDirection(float yawDeg, float pitchDeg) const {
    QMatrix4x4 rotation;
    rotation.rotate(yawDeg, 0.0f, 1.0f, 0.0f);
    rotation.rotate(pitchDeg, 1.0f, 0.0f, 0.0f);
    // Направление света хранится в мировых координатах, поэтому поворачиваем
    // его тем же вращением камеры, что и сцену.
    return rotation.mapVector(lightDirection_);
}
