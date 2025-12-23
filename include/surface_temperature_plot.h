#pragma once

#include "surface_temperature_calculator.h"

#include <qwt/qwt_plot.h>

class QwtPlotCurve;
class QwtPlotGrid;
class QwtPlotMarker;

class SurfaceTemperaturePlot : public QwtPlot {
public:
    explicit SurfaceTemperaturePlot(QWidget *parent = nullptr);

    void setTemperatureSeries(const QVector<TemperatureRangePoint> &points);
    void clearSeries();

private:
    QwtPlotCurve *minimumCurve_;
    QwtPlotCurve *maximumCurve_;
    QwtPlotGrid *grid_;
    QwtPlotMarker *freezingMarker_;
};
