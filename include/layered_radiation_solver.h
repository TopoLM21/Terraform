#pragma once

#include "atmospheric_column.h"

#include <QVector>

class LayeredRadiationSolver {
public:
    // Потоки, которые радиационный солвер вычисляет для поверхности.
    struct SurfaceRadiativeFluxes {
        // Коротковолновый поток, поглощённый поверхностью (Вт/м²):
        // insolation * cloudTransmission * atmosphericTransmission * (1 - albedo).
        double shortwaveAbsorbedWPerM2 = 0.0;
        // Нисходящий длинноволновый поток от атмосферы к поверхности (backradiation, Вт/м²).
        double longwaveDownWPerM2 = 0.0;
    };

    explicit LayeredRadiationSolver(double timeStepSeconds);

    // Возвращает изменение температуры по слоям за шаг dt.
    QVector<double> solve(const AtmosphericColumn &column,
                          double insolationWPerM2,
                          double albedo,
                          double cloudShortwaveTransmission,
                          double surfaceTemperatureKelvin) const;

    // Версия solve, которая дополнительно возвращает потоки для поверхности.
    QVector<double> solve(const AtmosphericColumn &column,
                          double insolationWPerM2,
                          double albedo,
                          double cloudShortwaveTransmission,
                          double surfaceTemperatureKelvin,
                          SurfaceRadiativeFluxes &surfaceFluxes) const;

private:
    double timeStepSeconds_ = 0.0;
};
