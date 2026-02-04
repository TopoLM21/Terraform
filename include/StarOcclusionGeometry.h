#pragma once

#include <QPainterPath>
#include <QPointF>
#include <QVector3D>

class StarOcclusionGeometry {
public:
    struct Result {
        QPointF starCenter;
        double starRadiusPx = 0.0;
        bool hidden = false;
        bool hasClip = false;
        QPainterPath clipPath;
    };

    Result evaluate(const QPointF &screenCenter,
                    double sphereRadiusPx,
                    const QVector3D &screenDirection,
                    double starRadiusPx) const;
};
