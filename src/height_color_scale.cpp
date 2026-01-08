#include "height_color_scale.h"

#include <QtMath>

namespace {
const QVector<HeightColorStop> &defaultStops() {
    static const QVector<HeightColorStop> kStops = {
        // Hypsometric tints: greens for lowlands, yellow/brown for hills, white for high peaks.
        {0.0, QColor(34, 90, 38)},
        {0.25, QColor(72, 140, 70)},
        {0.45, QColor(120, 180, 95)},
        {0.58, QColor(200, 190, 120)},
        {0.72, QColor(190, 150, 95)},
        {0.85, QColor(150, 110, 75)},
        {0.94, QColor(205, 200, 195)},
        {1.0, QColor(250, 250, 250)}
    };
    return kStops;
}
} // namespace

const QVector<HeightColorStop> &heightColorStops() {
    return defaultStops();
}

QColor heightColorForRatio(double ratio) {
    const auto &stops = heightColorStops();
    if (stops.isEmpty()) {
        return QColor();
    }

    const double t = qBound(0.0, ratio, 1.0);
    for (int i = 1; i < stops.size(); ++i) {
        if (t <= stops[i].position) {
            const HeightColorStop &left = stops[i - 1];
            const HeightColorStop &right = stops[i];
            const double span = right.position - left.position;
            const double local = qFuzzyIsNull(span) ? 0.0 : (t - left.position) / span;

            auto mixChannel = [local](int from, int to) {
                return static_cast<int>(qRound(from + (to - from) * local));
            };

            return QColor(mixChannel(left.color.red(), right.color.red()),
                          mixChannel(left.color.green(), right.color.green()),
                          mixChannel(left.color.blue(), right.color.blue()));
        }
    }

    return stops.back().color;
}
