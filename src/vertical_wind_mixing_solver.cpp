#include "vertical_wind_mixing_solver.h"

#include <QtCore/QLoggingCategory>
#include <QtCore/QtMath>

#include <limits>

namespace {
constexpr double kMinThicknessMeters = 1.0e-3;
constexpr double kMaxDiffusionCourant = 0.5;
}

Q_LOGGING_CATEGORY(verticalWindMixingLog, "solar.atmosphere.vertical_mixing")

VerticalWindMixingSolver::VerticalWindMixingSolver(double mixingCoefficientKz)
    : mixingCoefficientKz_(mixingCoefficientKz) {}

void VerticalWindMixingSolver::setMixingCoefficient(double mixingCoefficientKz) {
    mixingCoefficientKz_ = mixingCoefficientKz;
}

double VerticalWindMixingSolver::mixingCoefficient() const {
    return mixingCoefficientKz_;
}

void VerticalWindMixingSolver::mix(AtmosphericColumn &column, double dtSeconds) const {
    if (dtSeconds <= 0.0 || mixingCoefficientKz_ <= 0.0) {
        return;
    }

    auto &layers = column.layers();
    const int layerCount = layers.size();
    if (layerCount < 2) {
        return;
    }

    QVector<double> windU;
    QVector<double> windV;
    QVector<double> thickness;
    windU.reserve(layerCount);
    windV.reserve(layerCount);
    thickness.reserve(layerCount);

    const bool logEnabled = verticalWindMixingLog().isInfoEnabled();
    double energyBefore = 0.0;
    double minThickness = std::numeric_limits<double>::max();
    for (const auto &layer : layers) {
        const double u = layer.windUMps();
        const double v = layer.windVMps();
        windU.push_back(u);
        windV.push_back(v);
        if (logEnabled) {
            // Кинетическая энергия на единицу массы: 0.5 * (u^2 + v^2).
            energyBefore += 0.5 * (u * u + v * v);
        }
        const double layerThickness = qMax(layer.thicknessMeters(), kMinThicknessMeters);
        thickness.push_back(layerThickness);
        minThickness = qMin(minThickness, layerThickness);
    }
    minThickness = qMax(minThickness, kMinThicknessMeters);

    // Ограничение диффузионного числа Куранта: dt <= C * dz^2 / Kz
    // для устойчивости явной схемы. Используем минимальную толщину слоя.
    const double maxStableDt = kMaxDiffusionCourant * minThickness * minThickness / mixingCoefficientKz_;
    if (maxStableDt <= 0.0) {
        return;
    }

    // Сохраняем полный шаг модели, но делим на подшаги для соблюдения
    // ограничения Куранта по вертикальной диффузии (физическая устойчивость явной схемы).
    const int substeps = qMax(1, static_cast<int>(qCeil(dtSeconds / maxStableDt)));
    const double stepDt = dtSeconds / static_cast<double>(substeps);

    const int interfaceCount = layerCount - 1;
    QVector<double> initialWindU;
    QVector<double> initialWindV;
    if (logEnabled) {
        initialWindU = windU;
        initialWindV = windV;
    }

    QVector<double> fluxU(interfaceCount);
    QVector<double> fluxV(interfaceCount);
    QVector<double> nextWindU(layerCount);
    QVector<double> nextWindV(layerCount);

    // Нижняя граница: выбираем no-slip (u, v = 0 у поверхности), чтобы учесть
    // торможение ветра о шероховатую поверхность. Верхняя граница: нулевой поток
    // (∂u/∂z = 0), что эквивалентно отсутствию потока импульса через верх.
    const double bottomBoundaryU = 0.0;
    const double bottomBoundaryV = 0.0;

    for (int step = 0; step < substeps; ++step) {
        for (int i = 0; i < interfaceCount; ++i) {
            const double dz = qMax(0.5 * (thickness.at(i) + thickness.at(i + 1)), kMinThicknessMeters);
            // Вертикальная диффузия по уравнению
            //   ∂u/∂t = ∂/∂z (Kz ∂u/∂z)
            // в дискретной форме: поток на границе слоев
            //   F = -Kz (u_{k+1} - u_k) / Δz,
            // а изменение в слое: Δu = (F_{k-1/2} - F_{k+1/2}) * dt / h_k.
            fluxU[i] = -mixingCoefficientKz_ * (windU.at(i + 1) - windU.at(i)) / dz;
            fluxV[i] = -mixingCoefficientKz_ * (windV.at(i + 1) - windV.at(i)) / dz;
        }

        const double bottomDz = qMax(0.5 * thickness.at(0), kMinThicknessMeters);
        const double bottomFluxU = -mixingCoefficientKz_ * (windU.at(0) - bottomBoundaryU) / bottomDz;
        const double bottomFluxV = -mixingCoefficientKz_ * (windV.at(0) - bottomBoundaryV) / bottomDz;

        for (int k = 0; k < layerCount; ++k) {
            const double dzLayer = thickness.at(k);
            const double fluxDownU = (k > 0) ? fluxU.at(k - 1) : bottomFluxU;
            const double fluxUpU = (k < interfaceCount) ? fluxU.at(k) : 0.0;
            const double fluxDownV = (k > 0) ? fluxV.at(k - 1) : bottomFluxV;
            const double fluxUpV = (k < interfaceCount) ? fluxV.at(k) : 0.0;

            nextWindU[k] = windU.at(k) + (fluxDownU - fluxUpU) * stepDt / dzLayer;
            nextWindV[k] = windV.at(k) + (fluxDownV - fluxUpV) * stepDt / dzLayer;
        }

        windU.swap(nextWindU);
        windV.swap(nextWindV);
    }

    double energyAfter = 0.0;
    double maxDeltaU = 0.0;
    double maxDeltaV = 0.0;
    for (int k = 0; k < layerCount; ++k) {
        if (logEnabled) {
            energyAfter += 0.5 * (windU.at(k) * windU.at(k) + windV.at(k) * windV.at(k));
            maxDeltaU = qMax(maxDeltaU, qAbs(windU.at(k) - initialWindU.at(k)));
            maxDeltaV = qMax(maxDeltaV, qAbs(windV.at(k) - initialWindV.at(k)));
        }
        layers[k].setWindUMps(windU.at(k));
        layers[k].setWindVMps(windV.at(k));
    }

    if (logEnabled) {
        qCInfo(verticalWindMixingLog) << "Vertical wind mixing step"
                                      << "layers=" << layerCount
                                      << "dt=" << dtSeconds
                                      << "substeps=" << substeps
                                      << "energyBefore=" << energyBefore
                                      << "energyAfter=" << energyAfter
                                      << "maxDeltaU=" << maxDeltaU
                                      << "maxDeltaV=" << maxDeltaV;
    }
}
