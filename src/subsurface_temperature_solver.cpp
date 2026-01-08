#include "subsurface_temperature_solver.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kMinHeatCapacity = 1e-6;
} // namespace

SubsurfaceTemperatureSolver::SubsurfaceTemperatureSolver(
    const SubsurfaceGrid &grid,
    const QVector<double> &thermalConductivityByLayer,
    const QVector<double> &densityByLayer,
    const QVector<double> &specificHeatByLayer,
    SubsurfaceBottomBoundaryCondition bottomBoundary,
    double bottomTemperatureKelvin,
    double geothermalFluxWPerM2) {
    reset(grid,
          thermalConductivity,
          density,
          specificHeat,
          bottomBoundary,
          bottomTemperatureKelvin,
          geothermalFluxWPerM2);
}

void SubsurfaceTemperatureSolver::reset(const SubsurfaceGrid &grid,
                                        const QVector<double> &thermalConductivityByLayer,
                                        const QVector<double> &densityByLayer,
                                        const QVector<double> &specificHeatByLayer,
                                        SubsurfaceBottomBoundaryCondition bottomBoundary,
                                        double bottomTemperatureKelvin,
                                        double geothermalFluxWPerM2) {
    grid_ = grid;
    bottomBoundary_ = bottomBoundary;
    bottomTemperatureKelvin_ = bottomTemperatureKelvin;
    geothermalFluxWPerM2_ = geothermalFluxWPerM2;

    const int layers = grid_.layerCount();
    const auto normalizeProfile = [layers](const QVector<double> &values, double fallback) {
        QVector<double> result;
        result.resize(layers);
        for (int i = 0; i < layers; ++i) {
            const double value = (i < values.size()) ? values[i] : fallback;
            result[i] = qMax(1e-6, value);
        }
        return result;
    };
    thermalConductivityByLayer_ = normalizeProfile(thermalConductivityByLayer, 1.0);
    densityByLayer_ = normalizeProfile(densityByLayer, 1.0);
    specificHeatByLayer_ = normalizeProfile(specificHeatByLayer, 1.0);
    temperatures_.resize(layers);
}

void SubsurfaceTemperatureSolver::setInitialTemperature(double temperatureKelvin) {
    for (double &value : temperatures_) {
        value = temperatureKelvin;
    }
    if (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::FixedTemperature) {
        bottomTemperatureKelvin_ = temperatureKelvin;
    }
}

void SubsurfaceTemperatureSolver::setTemperatures(const QVector<double> &temperatures) {
    temperatures_ = temperatures;
    if (temperatures_.size() != grid_.layerCount()) {
        temperatures_.resize(grid_.layerCount());
    }
}

const QVector<double> &SubsurfaceTemperatureSolver::temperatures() const {
    return temperatures_;
}

double SubsurfaceTemperatureSolver::surfaceTemperatureKelvin() const {
    if (temperatures_.isEmpty()) {
        return 0.0;
    }
    return temperatures_.first();
}

double SubsurfaceTemperatureSolver::topLayerHeatCapacity() const {
    const auto &dz = grid_.layerThicknessesMeters();
    if (dz.isEmpty()) {
        return 0.0;
    }
    const double volumetricHeatCapacity = qMax(kMinHeatCapacity, density_ * specificHeat_);
    return volumetricHeatCapacity * dz.first();
}

void SubsurfaceTemperatureSolver::stepImplicit(double netSurfaceFlux, double dtSeconds) {
    if (temperatures_.isEmpty() || dtSeconds <= 0.0) {
        return;
    }

    // Неявная схема Backward Euler для уравнения теплопроводности:
    // C_i (T_i^{n+1} - T_i^n) / dt = F_{i-1/2}^{n+1} - F_{i+1/2}^{n+1},
    // F = -lambda * dT/dz. Это дает трёхдиагональную матрицу и решается методом прогонки.
    const int n = temperatures_.size();
    const auto &dz = grid_.layerThicknessesMeters();
    QVector<double> a(n, 0.0);
    QVector<double> b(n, 0.0);
    QVector<double> c(n, 0.0);
    QVector<double> d(n, 0.0);

    const auto interfaceConductivity = [this](int upperIndex, int lowerIndex) {
        const double kUpper = thermalConductivityByLayer_.value(upperIndex, 1.0);
        const double kLower = thermalConductivityByLayer_.value(lowerIndex, kUpper);
        const double sum = kUpper + kLower;
        if (sum <= 0.0) {
            return 0.0;
        }
        return 2.0 * kUpper * kLower / sum;
    };

    if (n == 1) {
        const double volumetricHeatCapacity =
            qMax(kMinHeatCapacity, densityByLayer_[0] * specificHeatByLayer_[0]);
        const double capacity = volumetricHeatCapacity * dz[0];
        const double dtInv = capacity / dtSeconds;
        double kBottom = 0.0;
        // Геотермальный поток задаётся как положительный вверх, добавляя тепло нижнему слою.
        const double bottomFlux = (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::Insulating)
            ? geothermalFluxWPerM2_
            : 0.0;
        if (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::FixedTemperature) {
            const double distBottom = 0.5 * dz[0];
            kBottom = thermalConductivityByLayer_[0] / qMax(1e-6, distBottom);
        }
        a[0] = 0.0;
        b[0] = dtInv + kBottom;
        c[0] = 0.0;
        d[0] = dtInv * temperatures_[0] + netSurfaceFlux + bottomFlux +
            kBottom * bottomTemperatureKelvin_;
        solveTridiagonal(a, b, c, d);
        temperatures_[0] = d[0];
        return;
    }

    for (int i = 0; i < n; ++i) {
        const double volumetricHeatCapacity =
            qMax(kMinHeatCapacity, densityByLayer_[i] * specificHeatByLayer_[i]);
        const double capacity = volumetricHeatCapacity * dz[i];
        const double dtInv = capacity / dtSeconds;

        if (i == 0) {
            const double distDown = 0.5 * (dz[0] + dz[1]);
            const double kDown =
                interfaceConductivity(0, 1) / qMax(1e-6, distDown);
            a[i] = 0.0;
            b[i] = dtInv + kDown;
            c[i] = -kDown;
            d[i] = dtInv * temperatures_[i] + netSurfaceFlux;
            continue;
        }

        if (i == n - 1) {
            const double distUp = 0.5 * (dz[n - 2] + dz[n - 1]);
            const double kUp =
                interfaceConductivity(n - 2, n - 1) / qMax(1e-6, distUp);
            double kBottom = 0.0;
            // Геотермальный поток задаётся как положительный вверх, добавляя тепло нижнему слою.
            const double bottomFlux =
                (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::Insulating)
                ? geothermalFluxWPerM2_
                : 0.0;
            if (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::FixedTemperature) {
                const double distBottom = 0.5 * dz[n - 1];
                kBottom = thermalConductivityByLayer_[n - 1] / qMax(1e-6, distBottom);
            }
            a[i] = -kUp;
            b[i] = dtInv + kUp + kBottom;
            c[i] = 0.0;
            d[i] = dtInv * temperatures_[i] + bottomFlux + kBottom * bottomTemperatureKelvin_;
            continue;
        }

        const double distUp = 0.5 * (dz[i - 1] + dz[i]);
        const double distDown = 0.5 * (dz[i] + dz[i + 1]);
        const double kUp = interfaceConductivity(i - 1, i) / qMax(1e-6, distUp);
        const double kDown = interfaceConductivity(i, i + 1) / qMax(1e-6, distDown);
        a[i] = -kUp;
        b[i] = dtInv + kUp + kDown;
        c[i] = -kDown;
        d[i] = dtInv * temperatures_[i];
    }

    solveTridiagonal(a, b, c, d);
    temperatures_ = d;
}

void SubsurfaceTemperatureSolver::setBottomBoundary(
    SubsurfaceBottomBoundaryCondition bottomBoundary,
    double bottomTemperatureKelvin) {
    bottomBoundary_ = bottomBoundary;
    bottomTemperatureKelvin_ = bottomTemperatureKelvin;
}

void SubsurfaceTemperatureSolver::solveTridiagonal(QVector<double> &a,
                                                   QVector<double> &b,
                                                   QVector<double> &c,
                                                   QVector<double> &d) const {
    const int n = b.size();
    if (n == 0) {
        return;
    }

    for (int i = 1; i < n; ++i) {
        const double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    d[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
    }
}
