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
    const double effectiveDt = qMin(dtSeconds, maxStableDt);
    if (effectiveDt <= 0.0) {
        return;
    }

    const int interfaceCount = layerCount - 1;
    QVector<double> fluxU(interfaceCount);
    QVector<double> fluxV(interfaceCount);

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

    // Нижняя граница: выбираем no-slip (u, v = 0 у поверхности), чтобы учесть
    // торможение ветра о шероховатую поверхность. Верхняя граница: нулевой поток
    // (∂u/∂z = 0), что эквивалентно отсутствию потока импульса через верх.
    const double bottomBoundaryU = 0.0;
    const double bottomBoundaryV = 0.0;
    const double bottomDz = qMax(0.5 * thickness.at(0), kMinThicknessMeters);
    const double bottomFluxU = -mixingCoefficientKz_ * (windU.at(0) - bottomBoundaryU) / bottomDz;
    const double bottomFluxV = -mixingCoefficientKz_ * (windV.at(0) - bottomBoundaryV) / bottomDz;

    double energyAfter = 0.0;
    double maxDeltaU = 0.0;
    double maxDeltaV = 0.0;
    for (int k = 0; k < layerCount; ++k) {
        const double dzLayer = thickness.at(k);
        const double fluxDownU = (k > 0) ? fluxU.at(k - 1) : bottomFluxU;
        const double fluxUpU = (k < interfaceCount) ? fluxU.at(k) : 0.0;
        const double fluxDownV = (k > 0) ? fluxV.at(k - 1) : bottomFluxV;
        const double fluxUpV = (k < interfaceCount) ? fluxV.at(k) : 0.0;

        const double updatedU = windU.at(k) + (fluxDownU - fluxUpU) * effectiveDt / dzLayer;
        const double updatedV = windV.at(k) + (fluxDownV - fluxUpV) * effectiveDt / dzLayer;
        if (logEnabled) {
            energyAfter += 0.5 * (updatedU * updatedU + updatedV * updatedV);
            maxDeltaU = qMax(maxDeltaU, qAbs(updatedU - windU.at(k)));
            maxDeltaV = qMax(maxDeltaV, qAbs(updatedV - windV.at(k)));
        }
        layers[k].setWindUMps(updatedU);
        layers[k].setWindVMps(updatedV);
    }

    if (logEnabled) {
        qCInfo(verticalWindMixingLog) << "Vertical wind mixing step"
                                      << "layers=" << layerCount
                                      << "dt=" << effectiveDt
                                      << "energyBefore=" << energyBefore
                                      << "energyAfter=" << energyAfter
                                      << "maxDeltaU=" << maxDeltaU
                                      << "maxDeltaV=" << maxDeltaV;
    }
}
