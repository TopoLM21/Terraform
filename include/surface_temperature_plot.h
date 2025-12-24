#pragma once

#include "surface_temperature_calculator.h"

#include <qwt/qwt_plot.h>
#include <QtCore/QString>

class QwtPlotCurve;
class QwtPlotGrid;
class QwtPlotMarker;
class TemperaturePlotTracker;

class SurfaceTemperaturePlot : public QwtPlot {
public:
    explicit SurfaceTemperaturePlot(QWidget *parent = nullptr);

    void setTemperatureSeries(const QVector<TemperatureRangePoint> &points,
                              const QVector<TemperatureSummaryPoint> &summaryPoints,
                              const QString &segmentLabel);
    void clearSeries();

private:
    QwtPlotCurve *minimumCurve_;
    QwtPlotCurve *maximumCurve_;
    QwtPlotCurve *meanAnnualCurve_;
    QwtPlotGrid *grid_;
    QwtPlotMarker *freezingMarker_;
    TemperaturePlotTracker *tracker_;
    QVector<TemperatureRangePoint> points_;
    QVector<TemperatureSummaryPoint> summaryPoints_;
    QString segmentLabel_;
};
