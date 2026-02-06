#include "precipitation_color_scale.h"

#include <QtMath>

namespace {
const QVector<PrecipitationColorStop> &defaultStops() {
    static const QVector<PrecipitationColorStop> kStops = {
        // Шкала осадков: от почти сухого состояния к ливням.
        // Светлые оттенки оставляют фон читаемым, а насыщенные сдвигаются в синие тона.
        {0.0, QColor(245, 247, 250)},
        {0.2, QColor(200, 230, 255)},
        {0.4, QColor(120, 190, 255)},
        {0.6, QColor(40, 120, 220)},
        {0.8, QColor(20, 80, 180)},
        {1.0, QColor(15, 40, 120)}
    };
    return kStops;
}

QColor mixColors(const QColor &left, const QColor &right, double t) {
    auto mixChannel = [t](int from, int to) {
        return static_cast<int>(qRound(from + (to - from) * t));
    };

    return QColor(mixChannel(left.red(), right.red()),
                  mixChannel(left.green(), right.green()),
                  mixChannel(left.blue(), right.blue()));
}
} // namespace

const QVector<PrecipitationColorStop> &precipitationColorStops() {
    return defaultStops();
}

QColor precipitationColorForRatio(double ratio) {
    const auto &stops = precipitationColorStops();
    if (stops.isEmpty()) {
        return QColor();
    }

    const double t = qBound(0.0, ratio, 1.0);
    for (int i = 1; i < stops.size(); ++i) {
        if (t <= stops[i].position) {
            const PrecipitationColorStop &left = stops[i - 1];
            const PrecipitationColorStop &right = stops[i];
            const double span = right.position - left.position;
            const double local = qFuzzyIsNull(span) ? 0.0 : (t - left.position) / span;
            return mixColors(left.color, right.color, local);
        }
    }

    return stops.back().color;
}
