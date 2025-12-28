#pragma once

#include "equal_area_projection.h"

class MollweideProjection final : public EqualAreaProjection {
public:
    QPointF project(double latitudeDeg, double longitudeDeg) const override;
    QPointF unproject(double x, double y) const override;
};
