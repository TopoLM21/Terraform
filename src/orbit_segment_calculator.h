#pragma once

#include <QtCore/QVector>

struct OrbitSegment {
    int index = 0;
    double meanAnomalyRadians = 0.0;
    double trueAnomalyRadians = 0.0;
    double distanceAU = 0.0;
};

class OrbitSegmentCalculator {
public:
    OrbitSegmentCalculator(double semiMajorAxisAU, double eccentricity);

    QVector<OrbitSegment> segments(int count = 12) const;

private:
    double semiMajorAxisAU_ = 0.0;
    double eccentricity_ = 0.0;
};
