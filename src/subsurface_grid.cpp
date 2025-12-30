#include "subsurface_grid.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr int kMaxIterations = 64;
} // namespace

SubsurfaceGrid::SubsurfaceGrid(int layerCount,
                               double topLayerThicknessMeters,
                               double bottomDepthMeters) {
    rebuild(layerCount, topLayerThicknessMeters, bottomDepthMeters);
}

void SubsurfaceGrid::rebuild(int layerCount,
                             double topLayerThicknessMeters,
                             double bottomDepthMeters) {
    layerThicknessesMeters_.clear();
    layerDepthsMeters_.clear();

    const int safeLayerCount = qMax(0, layerCount);
    const double safeBottomDepth = qMax(0.0, bottomDepthMeters);
    bottomDepthMeters_ = safeBottomDepth;

    if (safeLayerCount == 0 || safeBottomDepth <= 0.0) {
        return;
    }

    layerThicknessesMeters_.resize(safeLayerCount);
    layerDepthsMeters_.resize(safeLayerCount);

    if (safeLayerCount == 1) {
        layerThicknessesMeters_[0] = safeBottomDepth;
        layerDepthsMeters_[0] = 0.5 * safeBottomDepth;
        return;
    }

    const double safeTopThickness = qMax(1e-6, topLayerThicknessMeters);
    const double uniformThickness = safeBottomDepth / static_cast<double>(safeLayerCount);

    if (safeTopThickness * safeLayerCount >= safeBottomDepth) {
        // Если заданная верхняя толщина слишком велика, используем равномерную сетку.
        for (int i = 0; i < safeLayerCount; ++i) {
            layerThicknessesMeters_[i] = uniformThickness;
        }
    } else {
        // Логарифмическая (геометрическая) сетка: толщина растет с глубиной,
        // чтобы слой около поверхности был тонким, а нижние слои — более грубые.
        auto sumForRatio = [safeLayerCount, safeTopThickness](double ratio) {
            if (qFuzzyCompare(ratio, 1.0)) {
                return safeTopThickness * static_cast<double>(safeLayerCount);
            }
            return safeTopThickness * (std::pow(ratio, safeLayerCount) - 1.0) / (ratio - 1.0);
        };

        double low = 1.0;
        double high = qMax(1.0, safeBottomDepth / safeTopThickness);
        for (int iter = 0; iter < kMaxIterations; ++iter) {
            const double mid = 0.5 * (low + high);
            const double sum = sumForRatio(mid);
            if (sum < safeBottomDepth) {
                low = mid;
            } else {
                high = mid;
            }
        }

        const double ratio = 0.5 * (low + high);
        for (int i = 0; i < safeLayerCount; ++i) {
            layerThicknessesMeters_[i] = safeTopThickness * std::pow(ratio, i);
        }
    }

    double accumulated = 0.0;
    for (int i = 0; i < safeLayerCount; ++i) {
        if (i == safeLayerCount - 1) {
            layerThicknessesMeters_[i] = qMax(0.0, safeBottomDepth - accumulated);
        }
        const double thickness = layerThicknessesMeters_[i];
        layerDepthsMeters_[i] = accumulated + 0.5 * thickness;
        accumulated += thickness;
    }
}

int SubsurfaceGrid::layerCount() const {
    return layerThicknessesMeters_.size();
}

double SubsurfaceGrid::bottomDepthMeters() const {
    return bottomDepthMeters_;
}

const QVector<double> &SubsurfaceGrid::layerThicknessesMeters() const {
    return layerThicknessesMeters_;
}

const QVector<double> &SubsurfaceGrid::layerDepthsMeters() const {
    return layerDepthsMeters_;
}
