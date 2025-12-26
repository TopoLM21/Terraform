#include "atmospheric_lapse_rate_model.h"

#include "surface_temperature_calculator.h"

#include <QtCore/QtMath>

#include <algorithm>
#include <cmath>

namespace {
constexpr double kKelvinOffset = 273.15;
constexpr double kPascalPerAtm = 101325.0;
constexpr double kUniversalGasConstant = 8.314462618;
constexpr double kReferenceCpDry = 1004.0;
constexpr double kReferenceCpWet = 1850.0;
constexpr double kWaterMolarMassKgPerMol = 0.018015;
constexpr double kDryAirMolarMassKgPerMol = 0.02897;
constexpr double kLatentHeatVaporization = 2.5e6;
constexpr double kMinLapseRate = 0.0015;
constexpr double kMaxLapseRate = 0.012;
constexpr double kMinAdjustmentHeightMeters = 300.0;
constexpr double kMaxAdjustmentHeightMeters = 9000.0;
}  // namespace

AtmosphericLapseRateModel::AtmosphericLapseRateModel(double atmospherePressureAtm,
                                                     const AtmosphereComposition &atmosphere,
                                                     double surfaceGravity)
    : atmospherePressureAtm_(atmospherePressureAtm),
      atmosphere_(atmosphere),
      surfaceGravity_(surfaceGravity) {}

TemperatureRangePoint AtmosphericLapseRateModel::applyLapseRate(
    const TemperatureRangePoint &point) const {
    if (atmospherePressureAtm_ <= 0.0 || surfaceGravity_ <= 0.0) {
        return point;
    }

    const double meanTemperature = qMax(1.0, point.meanDailyKelvin);
    const double lapseRate = moistAdiabaticLapseRate(meanTemperature);
    const double adjustmentHeight = surfaceAdjustmentHeightMeters(meanTemperature);
    // Температуры из SurfaceTemperatureCalculator::radiativeBalanceByLatitudeForSegment
    // относятся к поверхности (слой layers[0]), поэтому вертикальную поправку не применяем.
    // Если бы это была эффективная температура излучения на высоте, к поверхности надо
    // прибавлять delta = Γ * Δz (положительный знак, потому что температура растет вниз).
    const double delta = lapseRate * adjustmentHeight;
    const bool applyVerticalCorrection = false;

    if (!applyVerticalCorrection) {
        return point;
    }

    TemperatureRangePoint adjusted = point;
    adjusted.meanDailyKelvin = qMax(1.0, adjusted.meanDailyKelvin + delta);
    adjusted.meanDayKelvin = qMax(1.0, adjusted.meanDayKelvin + delta);
    adjusted.meanNightKelvin = qMax(1.0, adjusted.meanNightKelvin + delta);
    adjusted.minimumKelvin = qMax(1.0, adjusted.minimumKelvin + delta);
    adjusted.maximumKelvin = qMax(1.0, adjusted.maximumKelvin + delta);
    adjusted.minimumCelsius = adjusted.minimumKelvin - kKelvinOffset;
    adjusted.maximumCelsius = adjusted.maximumKelvin - kKelvinOffset;
    adjusted.meanDailyCelsius = adjusted.meanDailyKelvin - kKelvinOffset;
    adjusted.meanDayCelsius = adjusted.meanDayKelvin - kKelvinOffset;
    adjusted.meanNightCelsius = adjusted.meanNightKelvin - kKelvinOffset;
    return adjusted;
}

double AtmosphericLapseRateModel::meanMolarMassKgPerMol() const {
    const auto fractions = atmosphere_.fractions();
    if (fractions.isEmpty()) {
        return kDryAirMolarMassKgPerMol;
    }

    double totalShare = 0.0;
    double meanMolarMass = 0.0;
    const auto gases = availableGases();
    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }
        const auto it = std::find_if(gases.begin(), gases.end(),
                                     [&fraction](const GasSpec &spec) {
                                         return spec.id == fraction.id;
                                     });
        if (it == gases.end()) {
            continue;
        }
        totalShare += fraction.share;
        meanMolarMass += fraction.share * (it->molarMass / 1000.0);
    }

    if (totalShare <= 0.0) {
        return kDryAirMolarMassKgPerMol;
    }
    return meanMolarMass / totalShare;
}

double AtmosphericLapseRateModel::specificGasConstant() const {
    // R_specific = R_universal / M.
    const double molarMass = qMax(1e-6, meanMolarMassKgPerMol());
    return kUniversalGasConstant / molarMass;
}

double AtmosphericLapseRateModel::specificHeatCp() const {
    const auto fractions = atmosphere_.fractions();
    if (fractions.isEmpty()) {
        return kReferenceCpDry;
    }

    double greenhouseShare = 0.0;
    double totalShare = 0.0;
    const auto gases = availableGases();
    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }
        const auto it = std::find_if(gases.begin(), gases.end(),
                                     [&fraction](const GasSpec &spec) {
                                         return spec.id == fraction.id;
                                     });
        if (it == gases.end()) {
            continue;
        }
        totalShare += fraction.share;
        if (it->isGreenhouse) {
            greenhouseShare += fraction.share;
        }
    }

    if (totalShare <= 0.0) {
        return kReferenceCpDry;
    }

    // Cp оцениваем как смесь сухого воздуха и легких паров: влажные/парниковые газы
    // увеличивают теплоёмкость, поэтому сдвигаем Cp к верхней границе.
    const double wetShare = qBound(0.0, greenhouseShare / totalShare, 1.0);
    return kReferenceCpDry + wetShare * (kReferenceCpWet - kReferenceCpDry);
}

double AtmosphericLapseRateModel::relativeHumidityEstimate() const {
    const auto fractions = atmosphere_.fractions();
    double h2oShare = 0.0;
    for (const auto &fraction : fractions) {
        if (fraction.id == QLatin1String("h2o")) {
            h2oShare = fraction.share;
            break;
        }
    }
    // Грубая оценка относительной влажности из доли водяного пара.
    return qBound(0.0, h2oShare * 4.0, 1.0);
}

double AtmosphericLapseRateModel::saturationVaporPressureAtm(double temperatureKelvin) const {
    // Формула Тетенса для насыщенного давления водяного пара:
    // e_s(T) = 6.112 * exp(17.67 * (T - 273.15) / (T - 29.65)) [гПа].
    const double temperatureCelsius = temperatureKelvin - kKelvinOffset;
    const double exponent = 17.67 * temperatureCelsius / (temperatureCelsius + 243.5);
    const double saturationHpa = 6.112 * std::exp(exponent);
    const double saturationPascal = saturationHpa * 100.0;
    return saturationPascal / kPascalPerAtm;
}

double AtmosphericLapseRateModel::moistAdiabaticLapseRate(double temperatureKelvin) const {
    const double cp = specificHeatCp();
    const double rSpecific = specificGasConstant();
    const double gravity = surfaceGravity_;

    const double dryLapse = dryAdiabaticLapseRate();

    // Влажный градиент через латентную теплоту:
    // Γ_m = Γ_d * (1 + L_v * q_s / (R * T)) / (1 + L_v^2 * q_s / (c_p * R * T^2)).
    // q_s — массовая доля насыщенного водяного пара.
    const double saturationPressureAtm = saturationVaporPressureAtm(temperatureKelvin);
    const double saturationPressure = saturationPressureAtm * kPascalPerAtm;
    const double totalPressure = qMax(1.0, atmospherePressureAtm_ * kPascalPerAtm);
    const double mixingRatio =
        0.622 * saturationPressure / qMax(1.0, totalPressure - saturationPressure);
    const double relativeHumidity = relativeHumidityEstimate();
    const double saturationMixingRatio = mixingRatio * relativeHumidity;

    const double numerator = 1.0 + (kLatentHeatVaporization * saturationMixingRatio) /
                                       (rSpecific * temperatureKelvin);
    const double denominator =
        1.0 + (kLatentHeatVaporization * kLatentHeatVaporization *
               saturationMixingRatio) /
                  (cp * rSpecific * temperatureKelvin * temperatureKelvin);
    const double moistLapse = dryLapse * (numerator / qMax(1e-6, denominator));
    const double maxLapse = qMin(kMaxLapseRate, dryLapse);
    return qBound(kMinLapseRate, moistLapse, maxLapse);
}

double AtmosphericLapseRateModel::dryAdiabaticLapseRate() const {
    const double cp = specificHeatCp();
    // Γ_d = g / c_p.
    return surfaceGravity_ / qMax(1.0, cp);
}

double AtmosphericLapseRateModel::scaleHeightMeters(double temperatureKelvin) const {
    // Масштабная высота: H = R_specific * T / g.
    // В гидростатике плотность убывает как ρ(z) = ρ0 * exp(-z / H),
    // поэтому H задаёт характерную толщину слоя атмосферы.
    return (specificGasConstant() * temperatureKelvin) / qMax(1.0, surfaceGravity_);
}

double AtmosphericLapseRateModel::surfaceAdjustmentHeightMeters(double temperatureKelvin) const {
    // Для вертикальной поправки используем долю масштабной высоты,
    // чтобы учесть связь с давлением и составом: чем ниже давление, тем тоньше слой.
    const double pressureFactor = std::sqrt(qMax(0.01, atmospherePressureAtm_));
    const double height = scaleHeightMeters(temperatureKelvin) * (0.35 + 0.35 * pressureFactor);
    return qBound(kMinAdjustmentHeightMeters, height, kMaxAdjustmentHeightMeters);
}
