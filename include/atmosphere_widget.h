#pragma once

#include "atmosphere_model.h"

#include <QtCore/QVector>
#include <QtWidgets/QGroupBox>

class QLabel;
class QTableWidget;
class AtmosphereChartWidget;
class QDoubleSpinBox;

class AtmosphereWidget : public QGroupBox {
    Q_OBJECT

public:
    explicit AtmosphereWidget(QWidget *parent = nullptr, bool showTable = true);

    void setPlanetParameters(double massEarths, double radiusKm);
    void setComposition(const AtmosphereComposition &composition);
    void clearPlanetParameters();
    AtmosphereComposition composition(bool includeZeroes = false) const;
    double minBottomLayerThicknessMeters() const;
    void setMinBottomLayerThicknessMeters(double meters);
    double minTopPressureAtm() const;
    void setMinTopPressureAtm(double pressureAtm);
    double verticalWindMixingCoefficient() const;
    void setVerticalWindMixingCoefficient(double coefficient);

signals:
    void compositionChanged(const AtmosphereComposition &composition);
    void minBottomLayerThicknessChanged(double meters);
    void minTopPressureChanged(double pressureAtm);
    void verticalWindMixingCoefficientChanged(double coefficient);

private:
    void populateTable();
    void updateAllShares();
    void updateAllPressures();
    void updateSummary();
    void normalizeMassItem(int row);
    void normalizePressureItem(int row);
    double parseMassText(const QString &text) const;
    double parsePressureText(const QString &text) const;
    double rowMassGigatons(int row) const;

    QTableWidget *table_ = nullptr;
    AtmosphereChartWidget *chartWidget_ = nullptr;
    QLabel *totalMassLabel_ = nullptr;
    QLabel *pressureLabel_ = nullptr;
    QLabel *meanMolarMassLabel_ = nullptr;
    QDoubleSpinBox *minBottomLayerThicknessSpinBox_ = nullptr;
    QDoubleSpinBox *minTopPressureSpinBox_ = nullptr;
    QDoubleSpinBox *verticalWindMixingSpinBox_ = nullptr;
    QVector<GasSpec> gases_;
    double planetMassEarths_ = 0.0;
    double planetRadiusKm_ = 0.0;
    bool hasPlanetParameters_ = false;
};
