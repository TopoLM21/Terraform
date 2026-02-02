#pragma once

#include <QtCore/QVector>

#include "FluidLayer.h"

class MixingModel {
public:
    using LayerStack = QVector<FluidLayer>;
    using FluidGrid = QVector<QVector<LayerStack>>;

    MixingModel() = default;
    MixingModel(double horizontalDiffusionM2PerS, double verticalDiffusionM2PerS);

    void setHorizontalDiffusionCoefficient(double horizontalDiffusionM2PerS);
    double horizontalDiffusionCoefficient() const;

    void setVerticalDiffusionCoefficient(double verticalDiffusionM2PerS);
    double verticalDiffusionCoefficient() const;

    // Горизонтальная диффузия между клетками сетки по уравнению:
    //   ∂σ/∂t = K_h ∇²σ,
    // где σ = ρ * h — масса на единицу площади. Обновление выполняется
    // по явной схеме с ограничением Courant: K_h * dt / L² ≤ 1/4,
    // чтобы гарантировать устойчивость и сохранение массы.
    void mixCells(FluidGrid &grid, double cellSizeMeters, double dtSeconds) const;

    // Вертикальное перетекание между слоями в колонке:
    //   ∂σ/∂t = ∂/∂z (K_v ∂σ/∂z),
    // где поток на границе слоев F = -K_v (σ_{k+1} - σ_k) / Δz.
    // Обновление консервативно по массе: Σσ остается неизменной.
    void mixLayers(QVector<FluidLayer> &layers, double dtSeconds) const;

private:
    double horizontalDiffusionM2PerS_ = 0.0;
    double verticalDiffusionM2PerS_ = 0.0;
};
