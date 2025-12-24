#pragma once

#include "surface_temperature_calculator.h"

#include <qwt/qwt_plot_picker.h>
#include <QtCore/QString>

class TemperaturePlotTracker : public QwtPlotPicker {
public:
    explicit TemperaturePlotTracker(QWidget *canvas);

    void setTemperatureSeries(const QVector<TemperatureRangePoint> &points,
                              const QVector<TemperatureSummaryPoint> &summaryPoints,
                              const QString &segmentLabel);
    void clearSeries();

protected:
    QwtText trackerTextF(const QPointF &pos) const override;

private:
    int nearestPointIndex(double latitudeDegrees) const;

    QVector<TemperatureRangePoint> points_;
    QVector<TemperatureSummaryPoint> summaryPoints_;
    QString segmentLabel_;
};
