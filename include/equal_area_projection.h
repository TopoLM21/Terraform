#pragma once

#include <QPointF>

class EqualAreaProjection {
public:
    virtual ~EqualAreaProjection() = default;

    virtual QPointF project(double latitudeDeg, double longitudeDeg) const = 0;
    virtual QPointF unproject(double x, double y) const = 0;
};
