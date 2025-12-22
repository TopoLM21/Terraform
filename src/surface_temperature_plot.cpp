#include "surface_temperature_plot.h"

#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_scale_draw.h>
#include <qwt_text.h>

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
    : QwtPlot(parent), curve_(new QwtPlotCurve(QStringLiteral("Температура"))),
      grid_(new QwtPlotGrid()) {
    setTitle(QStringLiteral("Температура поверхности по широтам"));
    setAxisTitle(QwtPlot::xBottom, QStringLiteral("Широта (°)"));
    setAxisTitle(QwtPlot::yLeft, QStringLiteral("Температура (K, °C)"));
    setAxisScale(QwtPlot::xBottom, -90.0, 90.0, 15.0);
    setAxisScaleDraw(QwtPlot::yLeft, new TemperatureScaleDraw());

    grid_->setMajorPen(QPen(palette().color(QPalette::Mid), 0.0, Qt::DashLine));
    grid_->attach(this);

    curve_->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    curve_->setPen(QPen(palette().color(QPalette::Highlight), 2.0));
    curve_->attach(this);
}

void SurfaceTemperaturePlot::setTemperatureSeries(const QVector<TemperaturePoint> &points) {
    QVector<QPointF> series;
    series.reserve(points.size());

    for (const auto &point : points) {
        series.push_back(QPointF(point.latitudeDegrees, point.temperatureKelvin));
    }

    curve_->setSamples(series);
    replot();
}

void SurfaceTemperaturePlot::clearSeries() {
    curve_->setSamples(QVector<QPointF>{});
    replot();
}
