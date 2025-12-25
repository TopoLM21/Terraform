#include "atmosphere_chart_widget.h"

#include <QtCharts/QBarCategoryAxis>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QValueAxis>
#include <QtGui/QPainter>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QVBoxLayout>

AtmosphereChartWidget::AtmosphereChartWidget(const QVector<GasSpec> &gases, QWidget *parent)
    : QWidget(parent)
    , gases_(gases) {
    axisSelector_ = new QComboBox(this);
    axisSelector_->addItem(QStringLiteral("Масса (Гт)"),
                           static_cast<int>(AtmosphereChartWidget::AxisMode::MassGigatons));
    axisSelector_->addItem(QStringLiteral("Давление (атм)"),
                           static_cast<int>(AtmosphereChartWidget::AxisMode::PressureAtm));

    chart_ = new QtCharts::QChart();
    chart_->legend()->hide();
    chart_->setTitle(QStringLiteral("Состав атмосферы"));

    chartView_ = new QtCharts::QChartView(chart_, this);
    chartView_->setRenderHint(QPainter::Antialiasing);
    // Минимальная высота нужна, чтобы график не сжимался в диалоге добавления планеты.
    chartView_->setMinimumHeight(200);

    axisX_ = new QtCharts::QBarCategoryAxis();
    axisY_ = new QtCharts::QValueAxis();
    axisY_->setLabelFormat(QStringLiteral("%.3f"));

    series_ = new QtCharts::QBarSeries(chart_);
    chart_->addSeries(series_);
    chart_->addAxis(axisX_, Qt::AlignBottom);
    chart_->addAxis(axisY_, Qt::AlignLeft);
    series_->attachAxis(axisX_);
    series_->attachAxis(axisY_);

    auto *layout = new QVBoxLayout(this);
    layout->addWidget(new QLabel(QStringLiteral("Ось Y:"), this));
    layout->addWidget(axisSelector_);
    layout->addWidget(chartView_, 1);
    setLayout(layout);

    connect(axisSelector_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
            [this]() { rebuildChart(); });

    rebuildChart();
}

void AtmosphereChartWidget::setPlanetParameters(double massEarths, double radiusKm) {
    if (massEarths <= 0.0 || radiusKm <= 0.0) {
        clearPlanetParameters();
        return;
    }
    planetMassEarths_ = massEarths;
    planetRadiusKm_ = radiusKm;
    hasPlanetParameters_ = true;
    rebuildChart();
}

void AtmosphereChartWidget::clearPlanetParameters() {
    planetMassEarths_ = 0.0;
    planetRadiusKm_ = 0.0;
    hasPlanetParameters_ = false;
    rebuildChart();
}

void AtmosphereChartWidget::setComposition(const AtmosphereComposition &composition) {
    composition_ = composition;
    rebuildChart();
}

AtmosphereChartWidget::AxisMode AtmosphereChartWidget::currentAxisMode() const {
    if (!axisSelector_) {
        return AtmosphereChartWidget::AxisMode::MassGigatons;
    }
    return static_cast<AtmosphereChartWidget::AxisMode>(axisSelector_->currentData().toInt());
}

void AtmosphereChartWidget::rebuildChart() {
    series_->clear();
    axisX_->clear();

    const AxisMode mode = currentAxisMode();
    const bool canShowPressure = hasPlanetParameters_;
    const bool showPressure = (mode == AxisMode::PressureAtm && canShowPressure);
    const QString yTitle = showPressure ? QStringLiteral("Давление (атм)")
                                        : QStringLiteral("Масса (Гт)");
    axisY_->setTitleText(yTitle);

    auto *set = new QtCharts::QBarSet(QStringLiteral("Газы"));
    double maxValue = 0.0;

    QStringList categories;
    for (const auto &gas : gases_) {
        const double massGigatons = composition_.massGigatons(gas.id);
        if (massGigatons <= 0.0) {
            continue;
        }

        double value = massGigatons;
        if (showPressure) {
            // Переводим массу газа (Гт -> кг), чтобы оценить вклад в давление по общей формуле.
            const double massKg = massGigatons * 1.0e12;
            value = calculatePressureAtmFromKg(massKg, planetMassEarths_, planetRadiusKm_);
        } else if (mode == AxisMode::PressureAtm && !canShowPressure) {
            // Если параметры планеты не заданы, отображаем нулевые столбцы вместо подстановки массы.
            value = 0.0;
        }

        categories.append(gas.displayName);
        *set << value;
        maxValue = qMax(maxValue, value);
    }

    series_->append(set);
    axisX_->append(categories);
    axisY_->setRange(0.0, maxValue > 0.0 ? maxValue * 1.2 : 1.0);
    axisY_->setLabelFormat(QStringLiteral("%.3f"));

    if (mode == AxisMode::PressureAtm && !canShowPressure) {
        axisY_->setTitleText(QStringLiteral("Давление (атм, требуется масса и радиус)"));
    }

    chart_->update();
}
