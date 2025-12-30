#include "temperature_plot_tracker.h"

#include <qwt-qt5/qwt_text.h>

#include <QtGui/QBrush>
#include <QtGui/QColor>

#include <algorithm>

TemperaturePlotTracker::TemperaturePlotTracker(QWidget *canvas)
    : QwtPlotPicker(canvas) {
    setTrackerMode(QwtPicker::AlwaysOn);
    setRubberBand(QwtPicker::VLineRubberBand);
}

void TemperaturePlotTracker::setTemperatureSeries(const QVector<TemperatureRangePoint> &points,
                                                  const QVector<TemperatureSummaryPoint> &summaryPoints,
                                                  const QString &segmentLabel,
                                                  RotationMode rotationMode) {
    points_ = points;
    summaryPoints_ = summaryPoints;
    segmentLabel_ = segmentLabel;
    rotationMode_ = rotationMode;
}

void TemperaturePlotTracker::clearSeries() {
    points_.clear();
    summaryPoints_.clear();
    segmentLabel_.clear();
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
    const QString axisLabel =
        (rotationMode_ == RotationMode::Normal)
            ? QStringLiteral("Широта")
            : QStringLiteral("Угол от подсолнечной точки");
    QString text = QStringLiteral("%1: %2°\nСегмент: мин %3 K (%4 °C)\n"
                                  "Сегмент: макс %5 K (%6 °C)")
                       .arg(axisLabel)
                       .arg(point.latitudeDegrees, 0, 'f', 0)
                       .arg(point.minimumKelvin, 0, 'f', 1)
                       .arg(point.minimumCelsius, 0, 'f', 1)
                       .arg(point.maximumKelvin, 0, 'f', 1)
                       .arg(point.maximumCelsius, 0, 'f', 1);

    if (index >= 0 && index < summaryPoints_.size()) {
        const auto &summaryPoint = summaryPoints_.at(index);
        text.append(QStringLiteral("\nСредняя за год: %1 K (%2 °C)\n"
                                   "Средняя за год (день): %3 K (%4 °C)\n"
                                   "Средняя за год (ночь): %5 K (%6 °C)\n"
                                   "Год: мин %7 K (%8 °C)\nГод: макс %9 K (%10 °C)")
                        .arg(summaryPoint.meanAnnualKelvin, 0, 'f', 1)
                        .arg(summaryPoint.meanAnnualCelsius, 0, 'f', 1)
                        .arg(summaryPoint.meanAnnualDayKelvin, 0, 'f', 1)
                        .arg(summaryPoint.meanAnnualDayCelsius, 0, 'f', 1)
                        .arg(summaryPoint.meanAnnualNightKelvin, 0, 'f', 1)
                        .arg(summaryPoint.meanAnnualNightCelsius, 0, 'f', 1)
                        .arg(summaryPoint.minimumKelvin, 0, 'f', 1)
                        .arg(summaryPoint.minimumCelsius, 0, 'f', 1)
                        .arg(summaryPoint.maximumKelvin, 0, 'f', 1)
                        .arg(summaryPoint.maximumCelsius, 0, 'f', 1));
    }
    if (!segmentLabel_.isEmpty()) {
        text.append(QStringLiteral("\n%1").arg(segmentLabel_));
    }

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
