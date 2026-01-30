#pragma once

#include "atmospheric_column.h"

#include <QVector>

class LayeredRadiationSolver {
public:
    explicit LayeredRadiationSolver(double timeStepSeconds);

    // Возвращает изменение температуры по слоям за шаг dt.
    QVector<double> solve(const AtmosphericColumn &column,
                          double insolationWPerM2,
                          double albedo,
                          double cloudShortwaveTransmission,
                          double surfaceTemperatureKelvin) const;

private:
    double timeStepSeconds_ = 0.0;
};
