#pragma once

#include "atmospheric_column.h"

#include <QVector>

// Количество ИК-полос в многополосной модели.
constexpr int kRadiationBandCount = 4;

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

    // Многополосные оптические толщины для одного слоя.
    struct BandOpticalDepths {
        double tauLw[kRadiationBandCount] = {};
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

    // Многополосная версия: для каждого слоя задаётся массив τ по ИК-полосам.
    // Двухпоточный перенос вычисляется независимо для каждой полосы,
    // суммируются потоки → корректное насыщение каждой полосы отдельно.
    QVector<double> solveMultiband(
        const AtmosphericColumn &column,
        const QVector<BandOpticalDepths> &bandTaus,
        double insolationWPerM2,
        double albedo,
        double cloudShortwaveTransmission,
        double surfaceTemperatureKelvin,
        SurfaceRadiativeFluxes &surfaceFluxes) const;

private:
    double timeStepSeconds_ = 0.0;
};
