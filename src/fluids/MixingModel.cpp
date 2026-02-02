#include "MixingModel.h"

#include <QtCore/QtMath>

#include <limits>

namespace {
constexpr double kMinThicknessMeters = 1.0e-6;
constexpr double kMaxDiffusionCourant2D = 0.25;
constexpr double kMaxDiffusionCourant1D = 0.5;
} // namespace

MixingModel::MixingModel(double horizontalDiffusionM2PerS, double verticalDiffusionM2PerS)
    : horizontalDiffusionM2PerS_(horizontalDiffusionM2PerS)
    , verticalDiffusionM2PerS_(verticalDiffusionM2PerS) {}

void MixingModel::setHorizontalDiffusionCoefficient(double horizontalDiffusionM2PerS) {
    horizontalDiffusionM2PerS_ = horizontalDiffusionM2PerS;
}

double MixingModel::horizontalDiffusionCoefficient() const {
    return horizontalDiffusionM2PerS_;
}

void MixingModel::setVerticalDiffusionCoefficient(double verticalDiffusionM2PerS) {
    verticalDiffusionM2PerS_ = verticalDiffusionM2PerS;
}

double MixingModel::verticalDiffusionCoefficient() const {
    return verticalDiffusionM2PerS_;
}

void MixingModel::mixCells(FluidGrid &grid, double cellSizeMeters, double dtSeconds) const {
    if (dtSeconds <= 0.0 || horizontalDiffusionM2PerS_ <= 0.0 || cellSizeMeters <= 0.0) {
        return;
    }

    const int rows = grid.size();
    if (rows == 0) {
        return;
    }
    const int cols = grid.front().size();
    if (cols == 0) {
        return;
    }

    const int layerCount = grid.front().front().size();
    if (layerCount == 0) {
        return;
    }

    // Ограничение Courant для 2D диффузии:
    //   K_h * dt / L² ≤ 1/4.
    const double maxStableDt = kMaxDiffusionCourant2D
        * cellSizeMeters * cellSizeMeters / horizontalDiffusionM2PerS_;
    if (maxStableDt <= 0.0) {
        return;
    }

    const int substeps = qMax(1, static_cast<int>(qCeil(dtSeconds / maxStableDt)));
    const double stepDt = dtSeconds / static_cast<double>(substeps);
    const double coeff = horizontalDiffusionM2PerS_ * stepDt
        / (cellSizeMeters * cellSizeMeters);

    QVector<QVector<double>> sigma(rows, QVector<double>(cols));
    QVector<QVector<double>> sigmaNext(rows, QVector<double>(cols));

    for (int layerIndex = 0; layerIndex < layerCount; ++layerIndex) {
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                const auto &layer = grid[r][c].at(layerIndex);
                sigma[r][c] = layer.columnMassKgPerM2();
            }
        }

        for (int step = 0; step < substeps; ++step) {
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    const double center = sigma[r][c];
                    const double north = (r > 0) ? sigma[r - 1][c] : center;
                    const double south = (r + 1 < rows) ? sigma[r + 1][c] : center;
                    const double west = (c > 0) ? sigma[r][c - 1] : center;
                    const double east = (c + 1 < cols) ? sigma[r][c + 1] : center;

                    // Явная схема: σ_new = σ + K_h * dt / L² * (Σ соседей - 4σ).
                    sigmaNext[r][c] = center + coeff * (north + south + west + east - 4.0 * center);
                }
            }
            sigma.swap(sigmaNext);
        }

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                auto &layer = grid[r][c][layerIndex];
                const double density = layer.densityKgPerM3();
                if (density <= 0.0) {
                    continue;
                }
                // Масса на единицу площади σ = ρ * h => h = σ / ρ.
                layer.setThicknessMeters(sigma[r][c] / density);
            }
        }
    }
}

void MixingModel::mixLayers(QVector<FluidLayer> &layers, double dtSeconds) const {
    if (dtSeconds <= 0.0 || verticalDiffusionM2PerS_ <= 0.0) {
        return;
    }

    const int layerCount = layers.size();
    if (layerCount < 2) {
        return;
    }

    QVector<double> sigma(layerCount);
    QVector<double> thickness(layerCount);
    double minDz = std::numeric_limits<double>::max();

    for (int i = 0; i < layerCount; ++i) {
        const double dz = qMax(layers.at(i).thicknessMeters(), kMinThicknessMeters);
        thickness[i] = dz;
        sigma[i] = layers.at(i).columnMassKgPerM2();
        if (i > 0) {
            const double interfaceDz = 0.5 * (thickness.at(i - 1) + dz);
            minDz = qMin(minDz, interfaceDz);
        }
    }

    if (!qIsFinite(minDz) || minDz <= 0.0) {
        return;
    }

    // Ограничение Courant для 1D диффузии:
    //   K_v * dt / Δz² ≤ 1/2.
    const double maxStableDt = kMaxDiffusionCourant1D * minDz * minDz / verticalDiffusionM2PerS_;
    if (maxStableDt <= 0.0) {
        return;
    }

    const int substeps = qMax(1, static_cast<int>(qCeil(dtSeconds / maxStableDt)));
    const double stepDt = dtSeconds / static_cast<double>(substeps);
    const int interfaceCount = layerCount - 1;

    QVector<double> flux(interfaceCount);
    QVector<double> sigmaNext(layerCount);

    for (int step = 0; step < substeps; ++step) {
        for (int i = 0; i < interfaceCount; ++i) {
            const double dz = 0.5 * (thickness.at(i) + thickness.at(i + 1));
            // Поток массы через границу слоев:
            //   F = -K_v * (σ_{k+1} - σ_k) / Δz.
            flux[i] = -verticalDiffusionM2PerS_ * (sigma.at(i + 1) - sigma.at(i)) / dz;
        }

        for (int k = 0; k < layerCount; ++k) {
            const double fluxDown = (k > 0) ? flux.at(k - 1) : 0.0;
            const double fluxUp = (k < interfaceCount) ? flux.at(k) : 0.0;
            // Консервативное обновление: σ^{n+1} = σ^n + (F_{k-1/2} - F_{k+1/2}) * dt.
            sigmaNext[k] = sigma.at(k) + (fluxDown - fluxUp) * stepDt;
        }
        sigma.swap(sigmaNext);
    }

    for (int i = 0; i < layerCount; ++i) {
        const double density = layers.at(i).densityKgPerM3();
        if (density <= 0.0) {
            continue;
        }
        layers[i].setThicknessMeters(sigma.at(i) / density);
    }
}
