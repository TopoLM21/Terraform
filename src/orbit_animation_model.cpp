#include "orbit_animation_model.h"

#include <QtCore/QtMath>

namespace {
constexpr double kTwoPi = 6.28318530717958647692;
}

OrbitAnimationModel::OrbitAnimationModel() = default;

void OrbitAnimationModel::reset(double semiMajorAxisAU,
                                double eccentricity,
                                double obliquityDegrees,
                                double perihelionArgumentDegrees,
                                int segmentsPerOrbit) {
    semiMajorAxisAU_ = semiMajorAxisAU;
    eccentricity_ = eccentricity;
    obliquityDegrees_ = obliquityDegrees;
    perihelionArgumentDegrees_ = perihelionArgumentDegrees;
    segmentsPerOrbit_ = qMax(1, segmentsPerOrbit);
    segmentIndex_ = 0;
}

void OrbitAnimationModel::advanceSegment() {
    if (segmentsPerOrbit_ <= 0) {
        segmentIndex_ = 0;
        return;
    }
    segmentIndex_ = (segmentIndex_ + 1) % segmentsPerOrbit_;
}

double OrbitAnimationModel::declinationDegrees() const {
    const double meanAnomaly =
        kTwoPi * static_cast<double>(segmentIndex_) / static_cast<double>(segmentsPerOrbit_);
    OrbitSegmentCalculator calculator(semiMajorAxisAU_, eccentricity_);
    const OrbitSegment segment = calculator.orbitAtMeanAnomaly(meanAnomaly);
    const double obliquityRadians = qDegreesToRadians(obliquityDegrees_);
    const double perihelionArgumentRadians = qDegreesToRadians(perihelionArgumentDegrees_);
    const double solarLongitude = segment.trueAnomalyRadians + perihelionArgumentRadians;
    // Сезонная деклинация: δ = asin(sin(наклон оси) * sin(истинная долгота звезды)).
    return qRadiansToDegrees(
        std::asin(std::sin(obliquityRadians) * std::sin(solarLongitude)));
}

double OrbitAnimationModel::distanceAU() const {
    const double meanAnomaly =
        kTwoPi * static_cast<double>(segmentIndex_) / static_cast<double>(segmentsPerOrbit_);
    OrbitSegmentCalculator calculator(semiMajorAxisAU_, eccentricity_);
    const OrbitSegment segment = calculator.orbitAtMeanAnomaly(meanAnomaly);
    return segment.distanceAU;
}
