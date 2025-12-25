#pragma once

#include "atmosphere_model.h"

#include <QtCore/QVector>
#include <QtWidgets/QGroupBox>

class QLabel;
class QTableWidget;
class AtmosphereChartWidget;

class AtmosphereWidget : public QGroupBox {
public:
    explicit AtmosphereWidget(QWidget *parent = nullptr, bool showTable = true);

    void setPlanetParameters(double massEarths, double radiusKm);
    void setComposition(const AtmosphereComposition &composition);
    void clearPlanetParameters();
    AtmosphereComposition composition(bool includeZeroes = false) const;

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
    QVector<GasSpec> gases_;
    double planetMassEarths_ = 0.0;
    double planetRadiusKm_ = 0.0;
    bool hasPlanetParameters_ = false;
};
