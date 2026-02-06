#pragma once

#include "atmospheric_column.h"

class VerticalMoistureMixingSolver {
public:
    explicit VerticalMoistureMixingSolver(double mixingCoefficientKz = 1.0);

    void setMixingCoefficient(double mixingCoefficientKz);
    double mixingCoefficient() const;

    void mix(AtmosphericColumn &column, double dtSeconds) const;

private:
    double mixingCoefficientKz_ = 1.0;
};
