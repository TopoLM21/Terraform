#include "orbit_segment_calculator.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr int kMaxKeplerIterations = 30;
constexpr double kKeplerTolerance = 1e-10;

double solveKeplerEquation(double meanAnomaly, double eccentricity) {
    double eccentricAnomaly = meanAnomaly;
    for (int i = 0; i < kMaxKeplerIterations; ++i) {
        const double delta =
            (eccentricAnomaly - eccentricity * std::sin(eccentricAnomaly) - meanAnomaly) /
            (1.0 - eccentricity * std::cos(eccentricAnomaly));
        eccentricAnomaly -= delta;
        if (std::abs(delta) < kKeplerTolerance) {
            break;
        }
    }
    return eccentricAnomaly;
}
}  // namespace

OrbitSegmentCalculator::OrbitSegmentCalculator(double semiMajorAxisAU, double eccentricity)
    : semiMajorAxisAU_(semiMajorAxisAU), eccentricity_(eccentricity) {}

QVector<OrbitSegment> OrbitSegmentCalculator::segments(int count) const {
    QVector<OrbitSegment> result;
    if (count <= 0) {
        return result;
    }

    result.reserve(count);
    const double step = 2.0 * kPi / static_cast<double>(count);

    for (int i = 0; i < count; ++i) {
        const double meanAnomaly = step * static_cast<double>(i);
        OrbitSegment segment = orbitAtMeanAnomaly(meanAnomaly);
        segment.index = i;
        result.push_back(segment);
    }

    return result;
}

OrbitSegment OrbitSegmentCalculator::orbitAtMeanAnomaly(double meanAnomalyRadians) const {
    const double eccentricAnomaly = solveKeplerEquation(meanAnomalyRadians, eccentricity_);
    const double trueAnomaly = 2.0 * std::atan2(std::sqrt(1.0 + eccentricity_) *
                                                   std::sin(eccentricAnomaly / 2.0),
                                               std::sqrt(1.0 - eccentricity_) *
                                                   std::cos(eccentricAnomaly / 2.0));
    const double distance = semiMajorAxisAU_ * (1.0 - eccentricity_ * std::cos(eccentricAnomaly));

    OrbitSegment segment;
    segment.meanAnomalyRadians = meanAnomalyRadians;
    segment.trueAnomalyRadians = trueAnomaly;
    segment.distanceAU = distance;
    return segment;
}
