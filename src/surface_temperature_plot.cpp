#include "surface_temperature_plot.h"
#include "temperature_plot_tracker.h"

#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_marker.h>
#include <qwt/qwt_scale_draw.h>
// #include <qwt/qwt_spline_curve_fitter.h>
#include <qwt/qwt_legend.h>
#include <qwt/qwt_legend_data.h>
#include <qwt/qwt_text.h>

#include <QtCore/QPair>
#include <QtGui/QColor>
#include <QtGui/QFont>
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

QString axisTitleForMode(RotationMode mode) {
    return (mode == RotationMode::Normal)
               ? QStringLiteral("Широта (°)")
               : QStringLiteral("Угол от подсолнечной точки (°)");
}

QString plotTitleForMode(RotationMode mode) {
    return (mode == RotationMode::Normal)
               ? QStringLiteral("Температура поверхности по широтам")
               : QStringLiteral("Температура поверхности по углу от подсолнечной точки");
}

QString temperatureContextLabel(bool hasAtmosphere) {
    return hasAtmosphere
               ? QStringLiteral("с учетом атмосферы")
               : QStringLiteral("поверхность");
}

QPair<double, double> axisRangeForMode(RotationMode mode) {
    return (mode == RotationMode::Normal)
               ? QPair<double, double>{-90.0, 90.0}
               : QPair<double, double>{0.0, 180.0};
}

double axisStepForMode(RotationMode mode) {
    return (mode == RotationMode::Normal) ? 15.0 : 30.0;
}
}  // namespace

SurfaceTemperaturePlot::SurfaceTemperaturePlot(QWidget *parent)
    : QwtPlot(parent),
      minimumCurve_(new QwtPlotCurve(QStringLiteral("Минимум за сутки"))),
      maximumCurve_(new QwtPlotCurve(QStringLiteral("Максимум за сутки"))),
      meanAnnualCurve_(new QwtPlotCurve(QStringLiteral("Средняя за год"))),
      meanAnnualDayCurve_(new QwtPlotCurve(QStringLiteral("Средняя за год (день)"))),
      meanAnnualNightCurve_(new QwtPlotCurve(QStringLiteral("Средняя за год (ночь)"))),
      grid_(new QwtPlotGrid()),
      freezingMarker_(new QwtPlotMarker()),
      tracker_(new TemperaturePlotTracker(canvas())) {
    setTitle(plotTitleForMode(rotationMode_));
    setAxisTitle(QwtPlot::xBottom, axisTitleForMode(rotationMode_));
    setAxisTitle(QwtPlot::yLeft, QStringLiteral("Температура (K, °C)"));
    const auto axisRange = axisRangeForMode(rotationMode_);
    setAxisScale(QwtPlot::xBottom, axisRange.first, axisRange.second,
                 axisStepForMode(rotationMode_));
    setAxisScaleDraw(QwtPlot::yLeft, new TemperatureScaleDraw());

    canvas()->setMouseTracking(true);

    auto *legend = new QwtLegend(this);
    legend->setDefaultItemMode(QwtLegendData::Checkable);
    insertLegend(legend, QwtPlot::RightLegend);
    connect(legend, &QwtLegend::checked, this, [this](const QVariant &info, bool checked, int index) {
        Q_UNUSED(index);
        QwtPlotItem *item = infoToItem(info);
        if (!item) {
            return;
        }
        item->setVisible(checked);
        replot();
    });

    grid_->setMajorPen(QPen(palette().color(QPalette::Mid), 0.0, Qt::DashLine));
    grid_->attach(this);

    freezingMarker_->setLineStyle(QwtPlotMarker::HLine);
    // Линия замерзания воды: 0 °C соответствует 273.15 K.
    freezingMarker_->setYValue(273.15);
    freezingMarker_->setLinePen(QPen(QColor(220, 30, 30), 2.0, Qt::DashDotLine));
    QwtText freezingLabel(QStringLiteral("0 °C"));
    QFont freezingFont = freezingLabel.font();
    freezingFont.setBold(true);
    freezingLabel.setFont(freezingFont);
    freezingLabel.setColor(QColor(220, 30, 30));
    freezingMarker_->setLabel(freezingLabel);
    freezingMarker_->setLabelAlignment(Qt::AlignRight | Qt::AlignTop);
    freezingMarker_->attach(this);

    minimumCurve_->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    minimumCurve_->setPen(QPen(palette().color(QPalette::Highlight), 2.0));
    minimumCurve_->attach(this);

    maximumCurve_->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    maximumCurve_->setPen(QPen(palette().color(QPalette::Highlight).lighter(130), 2.0));
    maximumCurve_->attach(this);

    meanAnnualCurve_->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    meanAnnualCurve_->setPen(QPen(QColor(0, 150, 0), 1.8));
    meanAnnualCurve_->attach(this);

    meanAnnualDayCurve_->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    meanAnnualDayCurve_->setPen(QPen(QColor(220, 140, 0), 1.6, Qt::DashLine));
    meanAnnualDayCurve_->attach(this);

    meanAnnualNightCurve_->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    meanAnnualNightCurve_->setPen(QPen(QColor(70, 110, 200), 1.6, Qt::DashLine));
    meanAnnualNightCurve_->attach(this);

    updateCurveTitles();
}

void SurfaceTemperaturePlot::updateCurveTitles() {
    const QString context = temperatureContextLabel(hasAtmosphere_);
    minimumCurve_->setTitle(QStringLiteral("Минимум за сутки (%1)").arg(context));
    maximumCurve_->setTitle(QStringLiteral("Максимум за сутки (%1)").arg(context));
    meanAnnualCurve_->setTitle(QStringLiteral("Средняя за год (%1)").arg(context));
    meanAnnualDayCurve_->setTitle(QStringLiteral("Средняя за год (день) (%1)").arg(context));
    meanAnnualNightCurve_->setTitle(QStringLiteral("Средняя за год (ночь) (%1)").arg(context));
}

void SurfaceTemperaturePlot::setSmoothingEnabled(bool enabled) {
    if (smoothingEnabled_ == enabled) {
        return;
    }

    smoothingEnabled_ = enabled;
    // Сглаживание временно отключено: оставляем флаг, но не применяем аппроксиматор.
    // const auto applyFitter = [enabled](QwtPlotCurve *curve) {
    //     curve->setCurveFitter(enabled ? new QwtSplineCurveFitter() : nullptr);
    // };
    //
    // applyFitter(minimumCurve_);
    // applyFitter(maximumCurve_);
    // applyFitter(meanAnnualCurve_);
    // applyFitter(meanAnnualDayCurve_);
    // applyFitter(meanAnnualNightCurve_);
    replot();
}

void SurfaceTemperaturePlot::setTemperatureSeries(const QVector<TemperatureRangePoint> &points,
                                                  const QVector<TemperatureSummaryPoint> &summaryPoints,
                                                  const QString &segmentLabel,
                                                  RotationMode rotationMode,
                                                  bool hasAtmosphere) {
    points_ = points;
    summaryPoints_ = summaryPoints;
    segmentLabel_ = segmentLabel;
    rotationMode_ = rotationMode;
    hasAtmosphere_ = hasAtmosphere;
    updateCurveTitles();
    tracker_->setTemperatureSeries(points_, summaryPoints_, segmentLabel_, rotationMode_);

    const auto axisRange = axisRangeForMode(rotationMode_);
    setAxisTitle(QwtPlot::xBottom, axisTitleForMode(rotationMode_));
    setAxisScale(QwtPlot::xBottom, axisRange.first, axisRange.second,
                 axisStepForMode(rotationMode_));
    if (!segmentLabel_.isEmpty()) {
        setTitle(QStringLiteral("%1 (%2) — %3")
                     .arg(plotTitleForMode(rotationMode_),
                          temperatureContextLabel(hasAtmosphere_),
                          segmentLabel_));
    } else {
        setTitle(QStringLiteral("%1 (%2)")
                     .arg(plotTitleForMode(rotationMode_),
                          temperatureContextLabel(hasAtmosphere_)));
    }

    QVector<QPointF> minimumSeries;
    QVector<QPointF> maximumSeries;
    QVector<QPointF> meanAnnualSeries;
    QVector<QPointF> meanAnnualDaySeries;
    QVector<QPointF> meanAnnualNightSeries;
    minimumSeries.reserve(points.size());
    maximumSeries.reserve(points.size());
    meanAnnualSeries.reserve(summaryPoints.size());
    meanAnnualDaySeries.reserve(summaryPoints.size());
    meanAnnualNightSeries.reserve(summaryPoints.size());

    for (const auto &point : points) {
        minimumSeries.push_back(QPointF(point.latitudeDegrees, point.minimumKelvin));
        maximumSeries.push_back(QPointF(point.latitudeDegrees, point.maximumKelvin));
    }

    for (const auto &summaryPoint : summaryPoints) {
        meanAnnualSeries.push_back(QPointF(summaryPoint.latitudeDegrees, summaryPoint.meanAnnualKelvin));
        meanAnnualDaySeries.push_back(
            QPointF(summaryPoint.latitudeDegrees, summaryPoint.meanAnnualDayKelvin));
        meanAnnualNightSeries.push_back(
            QPointF(summaryPoint.latitudeDegrees, summaryPoint.meanAnnualNightKelvin));
    }

    minimumCurve_->setSamples(minimumSeries);
    maximumCurve_->setSamples(maximumSeries);
    meanAnnualCurve_->setSamples(meanAnnualSeries);
    meanAnnualDayCurve_->setSamples(meanAnnualDaySeries);
    meanAnnualNightCurve_->setSamples(meanAnnualNightSeries);
    replot();
}

void SurfaceTemperaturePlot::clearSeries() {
    points_.clear();
    summaryPoints_.clear();
    segmentLabel_.clear();
    hasAtmosphere_ = false;
    tracker_->clearSeries();
    updateCurveTitles();
    setTitle(QStringLiteral("%1 (%2)")
                 .arg(plotTitleForMode(rotationMode_),
                      temperatureContextLabel(hasAtmosphere_)));
    minimumCurve_->setSamples(QVector<QPointF>{});
    maximumCurve_->setSamples(QVector<QPointF>{});
    meanAnnualCurve_->setSamples(QVector<QPointF>{});
    meanAnnualDayCurve_->setSamples(QVector<QPointF>{});
    meanAnnualNightCurve_->setSamples(QVector<QPointF>{});
    replot();
}
