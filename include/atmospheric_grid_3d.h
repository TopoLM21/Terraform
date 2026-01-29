#pragma once

#include "atmosphere_model.h"
#include "atmospheric_column.h"

#include <QVector>

class AtmosphericGrid3D {
public:
    void resizeColumns(int columnCount, int layerCount = 0);
    int columnCount() const;
    int layerCount() const;
    double minTopPressureAtm() const;

    QVector<AtmosphericColumn> &columns();
    const QVector<AtmosphericColumn> &columns() const;

    void initialize(const AtmosphereComposition &composition,
                    double planetMassEarths,
                    double radiusKm,
                    double baseTemperatureKelvin,
                    int columnCount,
                    int layerCount,
                    double minBottomLayerThicknessMeters,
                    double minTopPressureAtm = 0.0);

    bool updateLayerCountForTopPressure(double surfacePressureAtm,
                                        double minTopPressureAtm,
                                        double surfaceTemperatureKelvin,
                                        double gravityMps2,
                                        double specificGasConstant,
                                        double changeThresholdFraction);
    bool updateColumnLayerCountForTopPressure(int columnIndex,
                                              double surfacePressureAtm,
                                              double minTopPressureAtm,
                                              double surfaceTemperatureKelvin,
                                              double gravityMps2,
                                              double specificGasConstant,
                                              double changeThresholdFraction);

private:
    void rebuildLayersWithInterpolation(int newLayerCount, double topHeightMeters);
    bool rebuildColumnLayersWithInterpolation(AtmosphericColumn &column,
                                              int newLayerCount,
                                              double topHeightMeters);
    void updateMaxLayerCount();
    void ensureHistorySize();

    QVector<AtmosphericColumn> columns_;
    int layerCount_ = 0;
    double minBottomLayerThicknessMeters_ = 0.0;
    double minTopPressureAtm_ = 0.0;
    QVector<double> lastScaleHeightMetersByColumn_;
    QVector<double> lastSurfacePressureAtmByColumn_;
    QVector<double> lastTopHeightMetersByColumn_;
};
