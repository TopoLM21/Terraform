#include "subsurface_temperature_solver.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kMinHeatCapacity = 1e-6;
} // namespace

SubsurfaceTemperatureSolver::SubsurfaceTemperatureSolver(
    const SubsurfaceGrid &grid,
    double thermalConductivity,
    const ThermalConductivityModel &thermalConductivityModel,
    double density,
    double specificHeat,
    SubsurfaceBottomBoundaryCondition bottomBoundary,
    double bottomTemperatureKelvin,
    double geothermalFluxWPerM2) {
    reset(grid,
          thermalConductivity,
          thermalConductivityModel,
          density,
          specificHeat,
          bottomBoundary,
          bottomTemperatureKelvin,
          geothermalFluxWPerM2);
}

void SubsurfaceTemperatureSolver::reset(const SubsurfaceGrid &grid,
                                        double thermalConductivity,
                                        const ThermalConductivityModel &thermalConductivityModel,
                                        double density,
                                        double specificHeat,
                                        SubsurfaceBottomBoundaryCondition bottomBoundary,
                                        double bottomTemperatureKelvin,
                                        double geothermalFluxWPerM2) {
    grid_ = grid;
    thermalConductivity_ = qMax(1e-6, thermalConductivity);
    thermalConductivityModel_ = thermalConductivityModel;
    density_ = qMax(1e-6, density);
    specificHeat_ = qMax(1e-6, specificHeat);
    layerProperties_.clear();
    layerPropertiesProvider_ = LayerPropertiesProvider{};
    // Коэффициент температуропроводности: alpha = lambda / (rho * c).
    // Большая теплопроводность увеличивает скорость выравнивания,
    // а высокая плотность/удельная теплоемкость повышают тепловую инерцию.
    alpha_ = thermalConductivityAtTemperature(thermalConductivity_,
                                              thermalConductivityModel_,
                                              thermalConductivityModel_.referenceTemperatureKelvin) /
        (density_ * specificHeat_);
    bottomBoundary_ = bottomBoundary;
    bottomTemperatureKelvin_ = bottomTemperatureKelvin;
    geothermalFluxWPerM2_ = geothermalFluxWPerM2;

    const int layers = grid_.layerCount();
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

void SubsurfaceTemperatureSolver::setSurfaceLayerTemperatureKelvin(double temperatureKelvin) {
    if (temperatures_.isEmpty()) {
        return;
    }
    temperatures_.first() = temperatureKelvin;
}

void SubsurfaceTemperatureSolver::setTemperatures(const QVector<double> &temperatures) {
    temperatures_ = temperatures;
    const int targetLayers = grid_.layerCount();
    if (temperatures_.size() < targetLayers) {
        const double fillTemperature = temperatures_.isEmpty() ? bottomTemperatureKelvin_
                                                               : temperatures_.last();
        const int oldSize = temperatures_.size();
        temperatures_.resize(targetLayers);
        // Продлеваем профиль последней известной температурой (или нижней границей),
        // чтобы не вводить физически нереалистичные нули и избегать всплесков градиента.
        for (int i = oldSize; i < targetLayers; ++i) {
            temperatures_[i] = fillTemperature;
        }
    } else if (temperatures_.size() > targetLayers) {
        temperatures_.resize(targetLayers);
    }
}

void SubsurfaceTemperatureSolver::setLayerProperties(
    const QVector<SubsurfaceLayerProperties> &layerProperties) {
    layerProperties_ = layerProperties;
    layerPropertiesProvider_ = LayerPropertiesProvider{};
}

void SubsurfaceTemperatureSolver::setLayerPropertiesProvider(LayerPropertiesProvider provider) {
    layerPropertiesProvider_ = std::move(provider);
    layerProperties_.clear();
}

const QVector<double> &SubsurfaceTemperatureSolver::temperatures() const {
    return temperatures_;
}

const QVector<double> &SubsurfaceTemperatureSolver::layerThicknessesMeters() const {
    return grid_.layerThicknessesMeters();
}

const QVector<double> &SubsurfaceTemperatureSolver::layerDepthsMeters() const {
    return grid_.layerDepthsMeters();
}

double SubsurfaceTemperatureSolver::surfaceTemperatureKelvin() const {
    if (temperatures_.isEmpty()) {
        return 0.0;
    }
    return temperatures_.first();
}

double SubsurfaceTemperatureSolver::topLayerHeatCapacity() const {
    const auto &dz = grid_.layerThicknessesMeters();
    if (dz.isEmpty() || temperatures_.isEmpty()) {
        return 0.0;
    }
    const SubsurfaceLayerProperties layerProperties = layerPropertiesFor(0, temperatures_.first());
    const double volumetricHeatCapacity =
        qMax(kMinHeatCapacity,
             qMax(1e-6, layerProperties.density) * qMax(1e-6, layerProperties.specificHeat));
    return volumetricHeatCapacity * dz.first();
}

SubsurfaceBottomBoundaryCondition SubsurfaceTemperatureSolver::bottomBoundary() const {
    return bottomBoundary_;
}

double SubsurfaceTemperatureSolver::bottomTemperatureKelvin() const {
    return bottomTemperatureKelvin_;
}

void SubsurfaceTemperatureSolver::stepImplicit(double netSurfaceFlux, double dtSeconds) {
    if (temperatures_.isEmpty() || dtSeconds <= 0.0) {
        return;
    }

    // Неявная схема Backward Euler для уравнения теплопроводности:
    // C_i (T_i^{n+1} - T_i^n) / dt = F_{i-1/2}^{n+1} - F_{i+1/2}^{n+1},
    // F = -lambda * dT/dz. Это дает трёхдиагональную матрицу и решается методом прогонки.
    // Для температурозависимой теплопроводности используем k(T) на предыдущем шаге,
    // оценивая её по межслойной средней температуре.
    const int n = temperatures_.size();
    const auto &dz = grid_.layerThicknessesMeters();
    QVector<double> a(n, 0.0);
    QVector<double> b(n, 0.0);
    QVector<double> c(n, 0.0);
    QVector<double> d(n, 0.0);

    if (n == 1) {
        const SubsurfaceLayerProperties layerProperties =
            layerPropertiesFor(0, temperatures_[0]);
        const double volumetricHeatCapacity =
            qMax(kMinHeatCapacity,
                 qMax(1e-6, layerProperties.density) * qMax(1e-6, layerProperties.specificHeat));
        const double capacity = volumetricHeatCapacity * dz[0];
        const double dtInv = capacity / dtSeconds;
        double kBottom = 0.0;
        // Геотермальный поток задаётся как положительный вверх, добавляя тепло нижнему слою.
        const double bottomFlux = (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::Insulating)
            ? geothermalFluxWPerM2_
            : 0.0;
        if (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::FixedTemperature) {
            const double distBottom = 0.5 * dz[0];
            const double kInterface = conductivityAtTemperature(
                layerPropertiesFor(0, temperatures_[0]),
                0.5 * (temperatures_[0] + bottomTemperatureKelvin_));
            kBottom = kInterface / qMax(1e-6, distBottom);
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
        const SubsurfaceLayerProperties layerProperties =
            layerPropertiesFor(i, temperatures_[i]);
        const double volumetricHeatCapacity =
            qMax(kMinHeatCapacity,
                 qMax(1e-6, layerProperties.density) * qMax(1e-6, layerProperties.specificHeat));
        const double capacity = volumetricHeatCapacity * dz[i];
        const double dtInv = capacity / dtSeconds;

        if (i == 0) {
            const double distDown = 0.5 * (dz[0] + dz[1]);
            const double interfaceTemperature = 0.5 * (temperatures_[0] + temperatures_[1]);
            // Свойства материала (лёд/вода) определяются по фактической температуре
            // каждого слоя, а не по средней межслойной: иначе при T_avg > 273.15 K
            // замёрзший слой получит теплопроводность воды вместо льда.
            const double kInterface = 0.5 * (conductivityAtTemperature(
                                                 layerPropertiesFor(0, temperatures_[0]),
                                                 interfaceTemperature) +
                                             conductivityAtTemperature(
                                                 layerPropertiesFor(1, temperatures_[1]),
                                                 interfaceTemperature));
            const double kDown = kInterface / qMax(1e-6, distDown);
            a[i] = 0.0;
            b[i] = dtInv + kDown;
            c[i] = -kDown;
            d[i] = dtInv * temperatures_[i] + netSurfaceFlux;
            continue;
        }

        if (i == n - 1) {
            const double distUp = 0.5 * (dz[n - 2] + dz[n - 1]);
            const double interfaceTemperature = 0.5 * (temperatures_[n - 2] + temperatures_[n - 1]);
            const double kInterface = 0.5 * (conductivityAtTemperature(
                                                 layerPropertiesFor(n - 2, temperatures_[n - 2]),
                                                 interfaceTemperature) +
                                             conductivityAtTemperature(
                                                 layerPropertiesFor(n - 1, temperatures_[n - 1]),
                                                 interfaceTemperature));
            const double kUp = kInterface / qMax(1e-6, distUp);
            double kBottom = 0.0;
            // Геотермальный поток задаётся как положительный вверх, добавляя тепло нижнему слою.
            const double bottomFlux =
                (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::Insulating)
                ? geothermalFluxWPerM2_
                : 0.0;
            if (bottomBoundary_ == SubsurfaceBottomBoundaryCondition::FixedTemperature) {
                const double distBottom = 0.5 * dz[n - 1];
                const double interfaceTemperatureBottom =
                    0.5 * (temperatures_[n - 1] + bottomTemperatureKelvin_);
                const double kInterfaceBottom =
                    conductivityAtTemperature(layerPropertiesFor(n - 1, temperatures_[n - 1]),
                                               interfaceTemperatureBottom);
                kBottom = kInterfaceBottom / qMax(1e-6, distBottom);
            }
            a[i] = -kUp;
            b[i] = dtInv + kUp + kBottom;
            c[i] = 0.0;
            d[i] = dtInv * temperatures_[i] + bottomFlux + kBottom * bottomTemperatureKelvin_;
            continue;
        }

        const double distUp = 0.5 * (dz[i - 1] + dz[i]);
        const double distDown = 0.5 * (dz[i] + dz[i + 1]);
        const double upInterfaceTemperature = 0.5 * (temperatures_[i - 1] + temperatures_[i]);
        const double downInterfaceTemperature = 0.5 * (temperatures_[i] + temperatures_[i + 1]);
        const double kUpInterface =
            0.5 * (conductivityAtTemperature(layerPropertiesFor(i - 1, temperatures_[i - 1]),
                                             upInterfaceTemperature) +
                   conductivityAtTemperature(layerPropertiesFor(i, temperatures_[i]),
                                             upInterfaceTemperature));
        const double kDownInterface =
            0.5 * (conductivityAtTemperature(layerPropertiesFor(i, temperatures_[i]),
                                             downInterfaceTemperature) +
                   conductivityAtTemperature(layerPropertiesFor(i + 1, temperatures_[i + 1]),
                                             downInterfaceTemperature));
        const double kUp = kUpInterface / qMax(1e-6, distUp);
        const double kDown = kDownInterface / qMax(1e-6, distDown);
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

SubsurfaceLayerProperties SubsurfaceTemperatureSolver::layerPropertiesFor(
    int layerIndex,
    double temperatureKelvin) const {
    if (layerPropertiesProvider_) {
        return layerPropertiesProvider_(layerIndex, temperatureKelvin);
    }
    if (!layerProperties_.isEmpty()) {
        const int safeIndex = qBound(0, layerIndex, layerProperties_.size() - 1);
        return layerProperties_.at(safeIndex);
    }
    return {thermalConductivity_, thermalConductivityModel_, density_, specificHeat_};
}

double SubsurfaceTemperatureSolver::conductivityAtTemperature(
    const SubsurfaceLayerProperties &properties,
    double temperatureKelvin) const {
    const double baseConductivity = qMax(1e-6, properties.thermalConductivity);
    return thermalConductivityAtTemperature(baseConductivity,
                                            properties.thermalConductivityModel,
                                            temperatureKelvin);
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
