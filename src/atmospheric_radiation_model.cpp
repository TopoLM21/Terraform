#include "atmospheric_radiation_model.h"

#include <QtCore/QHash>
#include <QtCore/QtMath>

#include <algorithm>
#include <array>
#include <cmath>

namespace {
struct GasOpticalSpec {
    const char *id;
    double baseOpticalDepth;
    double saturationPressureAtm;
    double shortwaveShare;
};

struct OverlapSpec {
    const char *firstId;
    const char *secondId;
    double overlapCoefficient;
};

constexpr std::array<GasOpticalSpec, 6> kGasOpticalSpecs = {{
    {"co2", 1.45, 0.03, 0.08},
    {"h2o", 1.80, 0.02, 0.15},
    {"ch4", 1.10, 0.005, 0.05},
    {"nh3", 1.05, 0.004, 0.04},
    {"sf6", 2.80, 0.002, 0.02},
    {"nf3", 2.20, 0.003, 0.02},
}};

constexpr std::array<OverlapSpec, 4> kOverlapSpecs = {{
    {"co2", "h2o", 0.22},
    {"co2", "ch4", 0.12},
    {"ch4", "h2o", 0.10},
    {"co2", "nh3", 0.08},
}};

constexpr double kReferenceTemperatureKelvin = 288.0;
constexpr double kSaturationExponent = 0.65;
}  // namespace

AtmosphericRadiationModel::AtmosphericRadiationModel(const AtmosphereComposition &composition,
                                                     double pressureAtm,
                                                     double baseTemperatureKelvin)
    : composition_(composition),
      pressureAtm_(pressureAtm),
      baseTemperatureKelvin_(baseTemperatureKelvin) {
    computeOpticalDepths();
}

void AtmosphericRadiationModel::computeOpticalDepths() {
    effectiveOpticalDepth_ = 0.0;
    shortwaveOpticalDepth_ = 0.0;
    if (pressureAtm_ <= 0.0) {
        return;
    }

    QHash<QString, double> gasOpticalDepths;
    gasOpticalDepths.reserve(static_cast<int>(kGasOpticalSpecs.size()));

    const double temperatureScale =
        std::sqrt(qMax(1.0, baseTemperatureKelvin_) / kReferenceTemperatureKelvin);

    const auto fractions = composition_.fractions();
    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }

        const auto it = std::find_if(kGasOpticalSpecs.begin(), kGasOpticalSpecs.end(),
                                     [&fraction](const GasOpticalSpec &spec) {
                                         return fraction.id == QLatin1String(spec.id);
                                     });
        if (it == kGasOpticalSpecs.end()) {
            continue;
        }

        const double partialPressureAtm = pressureAtm_ * fraction.share;
        if (partialPressureAtm <= 0.0) {
            continue;
        }

        // Сатурирующая зависимость: при росте концентрации эффективность снижается.
        // Используем степенное насыщение по частичному давлению:
        // tau_i = tau0 * (p / (p + p_sat))^alpha.
        const double saturation =
            std::pow(partialPressureAtm / (partialPressureAtm + it->saturationPressureAtm),
                     kSaturationExponent);
        const double opticalDepth = it->baseOpticalDepth * saturation * temperatureScale;
        gasOpticalDepths.insert(fraction.id, opticalDepth);
        shortwaveOpticalDepth_ += it->shortwaveShare * opticalDepth;
    }

    if (gasOpticalDepths.isEmpty()) {
        shortwaveOpticalDepth_ = 0.0;
        return;
    }

    double transmissionProduct = 1.0;
    for (auto it = gasOpticalDepths.constBegin(); it != gasOpticalDepths.constEnd(); ++it) {
        const double clampedDepth = qBound(0.0, it.value(), 0.999);
        transmissionProduct *= (1.0 - clampedDepth);
    }

    // Базовое комбинирование полос поглощения: 1 - Π(1 - tau_i).
    const double combinedDepth = 1.0 - transmissionProduct;

    // Учет перекрытия полос: уменьшаем суммарный эффект для пар газов.
    double overlapPenalty = 1.0;
    for (const auto &spec : kOverlapSpecs) {
        const double first = gasOpticalDepths.value(QLatin1String(spec.firstId), 0.0);
        const double second = gasOpticalDepths.value(QLatin1String(spec.secondId), 0.0);
        if (first <= 0.0 || second <= 0.0) {
            continue;
        }
        const double overlap = spec.overlapCoefficient * std::sqrt(first * second);
        overlapPenalty *= (1.0 - qBound(0.0, overlap, 0.95));
    }

    overlapPenalty = qBound(0.05, overlapPenalty, 1.0);
    effectiveOpticalDepth_ = combinedDepth * overlapPenalty;
    shortwaveOpticalDepth_ *= overlapPenalty;
}

double AtmosphericRadiationModel::effectiveOpticalDepth() const {
    return effectiveOpticalDepth_;
}

double AtmosphericRadiationModel::incomingTransmission() const {
    // Закон Бугера-Ламберта: I = I0 * exp(-tau).
    return std::exp(-shortwaveOpticalDepth_);
}

double AtmosphericRadiationModel::outgoingTransmission() const {
    // Эффективная прозрачность в длинноволновом диапазоне.
    return std::exp(-effectiveOpticalDepth_);
}

double AtmosphericRadiationModel::applyIncomingFlux(double flux) const {
    return flux * incomingTransmission();
}

double AtmosphericRadiationModel::applyOutgoingFlux(double flux) const {
    return flux * outgoingTransmission();
}
