#include "mollweide_projection.h"

#include <QtMath>

namespace {
constexpr double kPi = 3.14159265358979323846;

double degreesToRadians(double degrees) {
    return degrees * kPi / 180.0;
}

double radiansToDegrees(double radians) {
    return radians * 180.0 / kPi;
}

double solveTheta(double latitudeRad) {
    // Решаем 2θ + sin(2θ) = π sin(φ) методом Ньютона.
    const double target = kPi * qSin(latitudeRad);
    double theta = latitudeRad;
    for (int i = 0; i < 8; ++i) {
        const double f = 2.0 * theta + qSin(2.0 * theta) - target;
        const double df = 2.0 + 2.0 * qCos(2.0 * theta);
        if (qFuzzyIsNull(df)) {
            break;
        }
        const double next = theta - f / df;
        if (qAbs(next - theta) < 1e-9) {
            theta = next;
            break;
        }
        theta = next;
    }
    return theta;
}
} // namespace

QPointF MollweideProjection::project(double latitudeDeg, double longitudeDeg) const {
    const double lat = degreesToRadians(latitudeDeg);
    const double lon = degreesToRadians(longitudeDeg);
    const double theta = solveTheta(lat);

    const double x = (2.0 * qSqrt(2.0) / kPi) * lon * qCos(theta);
    const double y = qSqrt(2.0) * qSin(theta);
    return QPointF(x, y);
}

QPointF MollweideProjection::unproject(double x, double y) const {
    // Обратное преобразование Mollweide.
    const double theta = qAsin(qBound(-1.0, y / qSqrt(2.0), 1.0));
    const double lat = qAsin(qBound(-1.0, (2.0 * theta + qSin(2.0 * theta)) / kPi, 1.0));
    const double cosTheta = qCos(theta);
    double lon = 0.0;
    if (!qFuzzyIsNull(cosTheta)) {
        lon = (kPi * x) / (2.0 * qSqrt(2.0) * cosTheta);
    }
    return QPointF(radiansToDegrees(lat), radiansToDegrees(lon));
}
