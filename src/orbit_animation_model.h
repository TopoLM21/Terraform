#pragma once

#include "orbit_segment_calculator.h"

class OrbitAnimationModel {
public:
    OrbitAnimationModel();

    void reset(double semiMajorAxisAU,
               double eccentricity,
               double obliquityDegrees,
               double perihelionArgumentDegrees,
               int segmentsPerOrbit = 12);
    void advanceSegment();

    double declinationDegrees() const;
    double distanceAU() const;

private:
    double semiMajorAxisAU_ = 1.0;
    double eccentricity_ = 0.0;
    double obliquityDegrees_ = 0.0;
    double perihelionArgumentDegrees_ = 0.0;
    int segmentsPerOrbit_ = 12;
    int segmentIndex_ = 0;
};
