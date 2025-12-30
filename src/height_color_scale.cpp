#include "height_color_scale.h"

#include <QtMath>

namespace {
const QVector<HeightColorStop> &defaultStops() {
    static const QVector<HeightColorStop> kStops = {
        {0.0, QColor(6, 22, 68)},
        {0.22, QColor(0, 64, 140)},
        {0.42, QColor(40, 120, 190)},
        {0.5, QColor(210, 200, 140)},
        {0.62, QColor(60, 160, 90)},
        {0.78, QColor(130, 110, 70)},
        {0.9, QColor(180, 170, 160)},
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
