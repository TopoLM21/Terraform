#pragma once

#include "atmosphere_model.h"
#include "atmospheric_column.h"

#include <QVector>

class AtmosphericGrid3D {
public:
    void resizeColumns(int columnCount, int layerCount = 0);
    int columnCount() const;
    int layerCount() const;

    QVector<AtmosphericColumn> &columns();
    const QVector<AtmosphericColumn> &columns() const;

    void initialize(const AtmosphereComposition &composition,
                    double planetMassEarths,
                    double radiusKm,
                    int columnCount,
                    int layerCount);

private:
    QVector<AtmosphericColumn> columns_;
    int layerCount_ = 0;
};
