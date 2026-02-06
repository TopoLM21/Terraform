#include "vertical_moisture_mixing_solver.h"

#include "atmosphere/EvaporationModel.h"

#include <QtCore/QLoggingCategory>
#include <QtCore/QtMath>

#include <limits>

namespace {
constexpr double kMinThicknessMeters = 1.0e-3;
constexpr double kMaxDiffusionCourant = 0.5;
}

Q_LOGGING_CATEGORY(verticalMoistureMixingLog, "solar.atmosphere.vertical_moisture_mixing")

VerticalMoistureMixingSolver::VerticalMoistureMixingSolver(double mixingCoefficientKz)
    : mixingCoefficientKz_(mixingCoefficientKz) {}

void VerticalMoistureMixingSolver::setMixingCoefficient(double mixingCoefficientKz) {
    mixingCoefficientKz_ = mixingCoefficientKz;
}

double VerticalMoistureMixingSolver::mixingCoefficient() const {
    return mixingCoefficientKz_;
}

void VerticalMoistureMixingSolver::mix(AtmosphericColumn &column,
                                       double dtSeconds,
                                       double minTopPressureAtm) const {
    if (dtSeconds <= 0.0 || mixingCoefficientKz_ <= 0.0) {
        return;
    }

    auto &layers = column.layers();
    const int layerCount = layers.size();
    if (layerCount < 2) {
        return;
    }

    QVector<double> vapor;
    QVector<double> liquid;
    QVector<double> ice;
    QVector<double> thickness;
    vapor.reserve(layerCount);
    liquid.reserve(layerCount);
    ice.reserve(layerCount);
    thickness.reserve(layerCount);

    double minThickness = std::numeric_limits<double>::max();
    for (const auto &layer : layers) {
        vapor.push_back(layer.waterVaporKgPerM2());
        liquid.push_back(layer.liquidWaterKgPerM2());
        ice.push_back(layer.iceWaterKgPerM2());
        const double layerThickness = qMax(layer.thicknessMeters(), kMinThicknessMeters);
        thickness.push_back(layerThickness);
        minThickness = qMin(minThickness, layerThickness);
    }
    minThickness = qMax(minThickness, kMinThicknessMeters);

    // Ограничение диффузионного числа Куранта: dt <= C * dz^2 / Kz
    // для устойчивости явной схемы вертикального переноса влаги.
    const double maxStableDt = kMaxDiffusionCourant * minThickness * minThickness / mixingCoefficientKz_;
    if (maxStableDt <= 0.0) {
        return;
    }

    // Делим шаг на подшаги, чтобы сохранить устойчивость по вертикали
    // при параметрическом обмене (Kz) и сохранить массу каждого компонента.
    const int substeps = qMax(1, static_cast<int>(qCeil(dtSeconds / maxStableDt)));
    const double stepDt = dtSeconds / static_cast<double>(substeps);

    const int interfaceCount = layerCount - 1;
    QVector<double> fluxVapor(interfaceCount);
    QVector<double> fluxLiquid(interfaceCount);
    QVector<double> fluxIce(interfaceCount);
    QVector<double> nextVapor(layerCount);
    QVector<double> nextLiquid(layerCount);
    QVector<double> nextIce(layerCount);

    const bool logEnabled = verticalMoistureMixingLog().isInfoEnabled();
    double totalBefore = 0.0;
    if (logEnabled) {
        for (int k = 0; k < layerCount; ++k) {
            totalBefore += qMax(0.0, vapor.at(k));
            totalBefore += qMax(0.0, liquid.at(k));
            totalBefore += qMax(0.0, ice.at(k));
        }
    }

    for (int step = 0; step < substeps; ++step) {
        for (int i = 0; i < interfaceCount; ++i) {
            const double dz = qMax(0.5 * (thickness.at(i) + thickness.at(i + 1)), kMinThicknessMeters);
            // Вертикальная диффузия влаги по уравнению
            //   ∂q/∂t = ∂/∂z (Kz ∂q/∂z),
            // где q — масса компонента на м² в слое.
            fluxVapor[i] = -mixingCoefficientKz_ * (vapor.at(i + 1) - vapor.at(i)) / dz;
            fluxLiquid[i] = -mixingCoefficientKz_ * (liquid.at(i + 1) - liquid.at(i)) / dz;
            fluxIce[i] = -mixingCoefficientKz_ * (ice.at(i + 1) - ice.at(i)) / dz;
        }

        if (minTopPressureAtm > 0.0 && interfaceCount > 0) {
            const double topPressureAtm = layers.last().pressureAtm();
            // Верхняя граница: не даём влаге уходить в разрежённые слои,
            // чтобы предотвратить неограниченный вынос влаги из модели.
            if (topPressureAtm > 0.0) {
                const double transitionScale = qMax(minTopPressureAtm, 1.0e-6);
                const double filter =
                    qBound(0.0, (topPressureAtm - minTopPressureAtm) / transitionScale, 1.0);
                const int topInterface = interfaceCount - 1;
                fluxVapor[topInterface] *= filter;
                fluxLiquid[topInterface] *= filter;
                fluxIce[topInterface] *= filter;
            }
        }

        for (int k = 0; k < layerCount; ++k) {
            const double dzLayer = thickness.at(k);
            const double fluxDownVapor = (k > 0) ? fluxVapor.at(k - 1) : 0.0;
            const double fluxUpVapor = (k < interfaceCount) ? fluxVapor.at(k) : 0.0;
            const double fluxDownLiquid = (k > 0) ? fluxLiquid.at(k - 1) : 0.0;
            const double fluxUpLiquid = (k < interfaceCount) ? fluxLiquid.at(k) : 0.0;
            const double fluxDownIce = (k > 0) ? fluxIce.at(k - 1) : 0.0;
            const double fluxUpIce = (k < interfaceCount) ? fluxIce.at(k) : 0.0;

            nextVapor[k] = vapor.at(k) + (fluxDownVapor - fluxUpVapor) * stepDt / dzLayer;
            nextLiquid[k] = liquid.at(k) + (fluxDownLiquid - fluxUpLiquid) * stepDt / dzLayer;
            nextIce[k] = ice.at(k) + (fluxDownIce - fluxUpIce) * stepDt / dzLayer;
        }

        vapor.swap(nextVapor);
        liquid.swap(nextLiquid);
        ice.swap(nextIce);
    }

    double totalAfter = 0.0;
    for (int k = 0; k < layerCount; ++k) {
        const double vaporValue = qMax(0.0, vapor.at(k));
        const double liquidValue = qMax(0.0, liquid.at(k));
        const double iceValue = qMax(0.0, ice.at(k));
        layers[k].setWaterVaporKgPerM2(vaporValue);
        layers[k].setLiquidWaterKgPerM2(liquidValue);
        layers[k].setIceWaterKgPerM2(iceValue);

        const double maxVapor = EvaporationModel::saturationWaterVaporKgPerM2(
            layers[k].temperatureKelvin(),
            layers[k].thicknessMeters());
        const double relativeHumidity =
            (maxVapor > 0.0) ? qBound(0.0, vaporValue / maxVapor, 1.0) : 0.0;
        layers[k].setRelativeHumidity(relativeHumidity);

        if (logEnabled) {
            totalAfter += vaporValue + liquidValue + iceValue;
        }
    }

    if (logEnabled) {
        qCInfo(verticalMoistureMixingLog) << "Vertical moisture mixing step"
                                          << "layers=" << layerCount
                                          << "dt=" << dtSeconds
                                          << "substeps=" << substeps
                                          << "massBefore=" << totalBefore
                                          << "massAfter=" << totalAfter;
    }
}
