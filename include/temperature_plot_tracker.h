#pragma once

#include "rotation_mode.h"
#include "surface_temperature_calculator.h"

#include <qwt/qwt_plot_picker.h>
#include <QtCore/QString>

class TemperaturePlotTracker : public QwtPlotPicker {
public:
    explicit TemperaturePlotTracker(QWidget *canvas);

    void setTemperatureSeries(const QVector<TemperatureRangePoint> &points,
                              const QVector<TemperatureSummaryPoint> &summaryPoints,
                              const QString &segmentLabel,
                              RotationMode rotationMode);
    void clearSeries();

protected:
    QwtText trackerTextF(const QPointF &pos) const override;

private:
    int nearestPointIndex(double latitudeDegrees) const;

    QVector<TemperatureRangePoint> points_;
    QVector<TemperatureSummaryPoint> summaryPoints_;
    QString segmentLabel_;
    RotationMode rotationMode_ = RotationMode::Normal;
};
