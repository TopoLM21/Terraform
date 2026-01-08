#include "star_color.h"

#include <QtMath>

namespace {
constexpr double kMinKelvin = 1000.0;
constexpr double kMaxKelvin = 40000.0;

int clampChannel(double value) {
    return static_cast<int>(qBound(0.0, value, 255.0));
}
} // namespace

QColor starColorFromTemperature(double temperatureK) {
    const double clampedK = qBound(kMinKelvin, temperatureK, kMaxKelvin);
    const double t = clampedK / 100.0;

    double red = 0.0;
    double green = 0.0;
    double blue = 0.0;

    if (t <= 66.0) {
        red = 255.0;
        green = 99.4708025861 * qLn(t) - 161.1195681661;
        if (t <= 19.0) {
            blue = 0.0;
        } else {
            blue = 138.5177312231 * qLn(t - 10.0) - 305.0447927307;
        }
    } else {
        red = 329.698727446 * qPow(t - 60.0, -0.1332047592);
        green = 288.1221695283 * qPow(t - 60.0, -0.0755148492);
        blue = 255.0;
    }

    return QColor(clampChannel(red), clampChannel(green), clampChannel(blue));
}
