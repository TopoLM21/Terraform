#pragma once

#include "atmosphere_model.h"

#include <QtWidgets/QWidget>

namespace QtCharts {
class QBarCategoryAxis;
class QBarSeries;
class QChart;
class QChartView;
class QValueAxis;
} // namespace QtCharts

class QComboBox;

class AtmosphereChartWidget : public QWidget {
public:
    explicit AtmosphereChartWidget(const QVector<GasSpec> &gases, QWidget *parent = nullptr);

    void setPlanetParameters(double massEarths, double radiusKm);
    void clearPlanetParameters();
    void setComposition(const AtmosphereComposition &composition);

private:
    enum class AxisMode {
        MassGigatons = 0,
        PressureAtm
    };

    AxisMode currentAxisMode() const;
    void rebuildChart();

    QVector<GasSpec> gases_;
    AtmosphereComposition composition_;
    bool hasPlanetParameters_ = false;
    double planetMassEarths_ = 0.0;
    double planetRadiusKm_ = 0.0;

    QComboBox *axisSelector_ = nullptr;
    QtCharts::QChart *chart_ = nullptr;
    QtCharts::QChartView *chartView_ = nullptr;
    QtCharts::QBarCategoryAxis *axisX_ = nullptr;
    QtCharts::QValueAxis *axisY_ = nullptr;
    QtCharts::QBarSeries *series_ = nullptr;
};
