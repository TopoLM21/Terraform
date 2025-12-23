#include "temperature_plot_tracker.h"

#include <qwt/qwt_text.h>

#include <QtGui/QBrush>
#include <QtGui/QColor>

#include <algorithm>

TemperaturePlotTracker::TemperaturePlotTracker(QWidget *canvas)
    : QwtPlotPicker(canvas) {
    setTrackerMode(QwtPicker::AlwaysOn);
    setRubberBand(QwtPicker::VLineRubberBand);
}

void TemperaturePlotTracker::setTemperatureSeries(const QVector<TemperatureRangePoint> &points) {
    points_ = points;
}

void TemperaturePlotTracker::clearSeries() {
    points_.clear();
}

QwtText TemperaturePlotTracker::trackerTextF(const QPointF &pos) const {
    if (points_.isEmpty()) {
        return QwtText();
    }

    const int index = nearestPointIndex(pos.x());
    if (index < 0) {
        return QwtText();
    }

    const auto &point = points_.at(index);
    const QString text = QStringLiteral("Широта: %1°\nМин: %2 K (%3 °C)\nМакс: %4 K (%5 °C)")
                             .arg(point.latitudeDegrees, 0, 'f', 0)
                             .arg(point.minimumKelvin, 0, 'f', 1)
                             .arg(point.minimumCelsius, 0, 'f', 1)
                             .arg(point.maximumKelvin, 0, 'f', 1)
                             .arg(point.maximumCelsius, 0, 'f', 1);

    QwtText trackerText(text);
    trackerText.setBackgroundBrush(QBrush(QColor(255, 255, 255, 220)));
    return trackerText;
}

int TemperaturePlotTracker::nearestPointIndex(double latitudeDegrees) const {
    if (points_.isEmpty()) {
        return -1;
    }

    auto it = std::lower_bound(points_.cbegin(), points_.cend(), latitudeDegrees,
                               [](const TemperatureRangePoint &point, double latitude) {
                                   return point.latitudeDegrees < latitude;
                               });

    if (it == points_.cbegin()) {
        return 0;
    }
    if (it == points_.cend()) {
        return points_.size() - 1;
    }

    const int rightIndex = static_cast<int>(std::distance(points_.cbegin(), it));
    const int leftIndex = rightIndex - 1;
    const double leftDiff =
        std::abs(points_.at(leftIndex).latitudeDegrees - latitudeDegrees);
    const double rightDiff =
        std::abs(points_.at(rightIndex).latitudeDegrees - latitudeDegrees);

    return (leftDiff <= rightDiff) ? leftIndex : rightIndex;
}
