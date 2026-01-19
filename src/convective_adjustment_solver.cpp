#include "convective_adjustment_solver.h"

#include "atmospheric_thermodynamics.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
bool hasWaterVapor(const AtmosphereComposition &composition) {
    const auto fractions = composition.fractions();
    for (const auto &fraction : fractions) {
        if (fraction.id == QLatin1String("h2o") && fraction.share > 0.0) {
            return true;
        }
    }
    return false;
}
}  // namespace

ConvectiveAdjustmentSolver::ConvectiveAdjustmentSolver(const AtmosphereComposition &composition,
                                                       double gravityMps2)
    : composition_(composition),
      hasWaterVapor_(hasWaterVapor(composition)),
      gravityMps2_(gravityMps2),
      specificHeatCp_(AtmosphericThermodynamics::specificHeatCp(composition)),
      specificGasConstant_(AtmosphericThermodynamics::specificGasConstant(composition)) {}

void ConvectiveAdjustmentSolver::adjust(AtmosphericColumn &column) const {
    auto &layers = column.layers();
    if (layers.size() < 2 || gravityMps2_ <= 0.0) {
        return;
    }

    // Проходим снизу вверх: при супер-адиабатическом градиенте выравниваем слой
    // до нейтральной сухой/влажной адиабаты, сохраняя энергию пары слоёв.
    for (int i = 0; i < layers.size() - 1; ++i) {
        auto &lower = layers[i];
        auto &upper = layers[i + 1];

        const double temperatureLower = lower.temperatureKelvin();
        const double temperatureUpper = upper.temperatureKelvin();
        const double pressureLowerAtm = lower.pressureAtm();
        const double pressureUpperAtm = upper.pressureAtm();

        if (temperatureLower <= 0.0 || temperatureUpper <= 0.0 ||
            pressureLowerAtm <= 0.0 || pressureUpperAtm <= 0.0) {
            continue;
        }

        const double dlnP = std::log(pressureUpperAtm / pressureLowerAtm);
        if (!std::isfinite(dlnP) || qFuzzyIsNull(dlnP)) {
            continue;
        }

        const double gradAd = adiabaticGradientDlnT(temperatureLower, pressureLowerAtm);
        const double adiabaticUpper = temperatureLower * std::exp(gradAd * dlnP);

        if (temperatureUpper >= adiabaticUpper) {
            continue;
        }

        const double heatCapacityLower = qMax(0.0, lower.heatCapacityJPerM2K());
        const double heatCapacityUpper = qMax(0.0, upper.heatCapacityJPerM2K());
        if (heatCapacityLower <= 0.0 || heatCapacityUpper <= 0.0) {
            // Если нет теплоёмкости, просто поднимаем верхний слой до нейтрали.
            upper.setTemperatureKelvin(adiabaticUpper);
            continue;
        }

        // Сохраняем энергию: C1*T1 + C2*T2 = C1*T1' + C2*T2',
        // при условии T2' = T1' * exp(∇_ad * ΔlnP).
        const double totalEnergy = heatCapacityLower * temperatureLower +
                                   heatCapacityUpper * temperatureUpper;
        const double factor = std::exp(gradAd * dlnP);
        const double denominator = heatCapacityLower + heatCapacityUpper * factor;
        if (denominator <= 0.0) {
            continue;
        }

        const double adjustedLower = totalEnergy / denominator;
        const double adjustedUpper = adjustedLower * factor;
        lower.setTemperatureKelvin(qMax(0.0, adjustedLower));
        upper.setTemperatureKelvin(qMax(0.0, adjustedUpper));
    }
}

double ConvectiveAdjustmentSolver::adiabaticGradientDlnT(double temperatureKelvin,
                                                         double pressureAtm) const {
    if (hasWaterVapor_) {
        return AtmosphericThermodynamics::adiabaticGradientDlnT(composition_,
                                                                pressureAtm,
                                                                temperatureKelvin,
                                                                gravityMps2_);
    }

    if (gravityMps2_ <= 0.0 || specificHeatCp_ <= 0.0) {
        return 0.0;
    }

    // Для сухой адиабаты: ∇ = Γ * R / g, Γ = g / c_p.
    const double dryLapseRate = gravityMps2_ / qMax(1.0, specificHeatCp_);
    return (dryLapseRate * specificGasConstant_) / gravityMps2_;
}
