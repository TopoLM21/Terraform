#include "temperature_color_scale.h"

#include <QtMath>

namespace {
const QVector<TemperatureColorStop> &defaultStops() {
    static const QVector<TemperatureColorStop> kStops = {
        // Холодная часть диапазона разбита более мелкими шагами,
        // чтобы тонкие различия температур были заметнее на карте и шкале.
        // У максимума добавлен более плотный stop для акцента на суточных колебаниях.
        {0.0, QColor(12, 24, 96)},
        {0.08, QColor(0, 64, 164)},
        {0.18, QColor(0, 120, 210)},
        {0.32, QColor(0, 180, 220)},
        {0.48, QColor(40, 210, 140)},
        {0.62, QColor(190, 220, 80)},
        {0.76, QColor(250, 180, 60)},
        {0.9, QColor(240, 80, 40)},
        {0.97, QColor(250, 100, 50)},
        {1.0, QColor(170, 0, 20)}
    };
    return kStops;
}
} // namespace

const QVector<TemperatureColorStop> &temperatureColorStops() {
    return defaultStops();
}

QColor temperatureColorForRatio(double ratio) {
    const auto &stops = temperatureColorStops();
    if (stops.isEmpty()) {
        return QColor();
    }

    const double t = qBound(0.0, ratio, 1.0);
    for (int i = 1; i < stops.size(); ++i) {
        if (t <= stops[i].position) {
            const TemperatureColorStop &left = stops[i - 1];
            const TemperatureColorStop &right = stops[i];
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
