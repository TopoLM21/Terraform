#pragma once

#include "surface_temperature_calculator.h"

#include <qwt/qwt_plot_picker.h>

class TemperaturePlotTracker : public QwtPlotPicker {
public:
    explicit TemperaturePlotTracker(QWidget *canvas);

    void setTemperatureSeries(const QVector<TemperatureRangePoint> &points);
    void clearSeries();

protected:
    QwtText trackerTextF(const QPointF &pos) const override;

private:
    int nearestPointIndex(double latitudeDegrees) const;

    QVector<TemperatureRangePoint> points_;
};
