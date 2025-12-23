#include "surface_temperature_plot.h"

#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_marker.h>
#include <qwt/qwt_scale_draw.h>
#include <qwt/qwt_text.h>

#include <QtGui/QPalette>

namespace {
class TemperatureScaleDraw : public QwtScaleDraw {
public:
    QwtText label(double value) const override {
        const double celsius = value - 273.15;
        return QwtText(QStringLiteral("%1 (%2)")
                           .arg(value, 0, 'f', 0)
                           .arg(celsius, 0, 'f', 0));
    }
};
}  // namespace

SurfaceTemperaturePlot::SurfaceTemperaturePlot(QWidget *parent)
    : QwtPlot(parent),
      minimumCurve_(new QwtPlotCurve(QStringLiteral("Минимум за сутки"))),
      maximumCurve_(new QwtPlotCurve(QStringLiteral("Максимум за сутки"))),
      grid_(new QwtPlotGrid()), freezingMarker_(new QwtPlotMarker()) {
    setTitle(QStringLiteral("Температура поверхности по широтам"));
    setAxisTitle(QwtPlot::xBottom, QStringLiteral("Широта (°)"));
    setAxisTitle(QwtPlot::yLeft, QStringLiteral("Температура (K, °C)"));
    setAxisScale(QwtPlot::xBottom, -90.0, 90.0, 15.0);
    setAxisScaleDraw(QwtPlot::yLeft, new TemperatureScaleDraw());

    grid_->setMajorPen(QPen(palette().color(QPalette::Mid), 0.0, Qt::DashLine));
    grid_->attach(this);

    freezingMarker_->setLineStyle(QwtPlotMarker::HLine);
    // Линия замерзания воды: 0 °C соответствует 273.15 K.
    freezingMarker_->setYValue(273.15);
    freezingMarker_->setLinePen(QPen(palette().color(QPalette::Mid), 1.0, Qt::DashLine));
    freezingMarker_->setLabel(QwtText(QStringLiteral("0 °C")));
    freezingMarker_->setLabelAlignment(Qt::AlignRight | Qt::AlignTop);
    freezingMarker_->attach(this);

    minimumCurve_->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    minimumCurve_->setPen(QPen(palette().color(QPalette::Highlight), 2.0));
    minimumCurve_->attach(this);

    maximumCurve_->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    maximumCurve_->setPen(QPen(palette().color(QPalette::Highlight).lighter(130), 2.0));
    maximumCurve_->attach(this);
}

void SurfaceTemperaturePlot::setTemperatureSeries(const QVector<TemperatureRangePoint> &points) {
    QVector<QPointF> minimumSeries;
    QVector<QPointF> maximumSeries;
    minimumSeries.reserve(points.size());
    maximumSeries.reserve(points.size());

    for (const auto &point : points) {
        minimumSeries.push_back(QPointF(point.latitudeDegrees, point.minimumKelvin));
        maximumSeries.push_back(QPointF(point.latitudeDegrees, point.maximumKelvin));
    }

    minimumCurve_->setSamples(minimumSeries);
    maximumCurve_->setSamples(maximumSeries);
    replot();
}

void SurfaceTemperaturePlot::clearSeries() {
    minimumCurve_->setSamples(QVector<QPointF>{});
    maximumCurve_->setSamples(QVector<QPointF>{});
    replot();
}
