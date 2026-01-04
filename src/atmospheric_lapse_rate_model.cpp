#include "atmospheric_lapse_rate_model.h"

#include "surface_temperature_calculator.h"
#include "atmospheric_thermodynamics.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kKelvinOffset = 273.15;
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
    return AtmosphericThermodynamics::meanMolarMassKgPerMol(atmosphere_);
}

double AtmosphericLapseRateModel::specificGasConstant() const {
    return AtmosphericThermodynamics::specificGasConstant(atmosphere_);
}

double AtmosphericLapseRateModel::specificHeatCp() const {
    return AtmosphericThermodynamics::specificHeatCp(atmosphere_);
}

double AtmosphericLapseRateModel::relativeHumidityEstimate() const {
    return AtmosphericThermodynamics::relativeHumidityEstimate(atmosphere_);
}

double AtmosphericLapseRateModel::saturationVaporPressureAtm(double temperatureKelvin) const {
    return AtmosphericThermodynamics::saturationVaporPressureAtm(temperatureKelvin);
}

double AtmosphericLapseRateModel::moistAdiabaticLapseRate(double temperatureKelvin) const {
    return AtmosphericThermodynamics::moistAdiabaticLapseRate(temperatureKelvin,
                                                              atmosphere_,
                                                              atmospherePressureAtm_,
                                                              surfaceGravity_);
}

double AtmosphericLapseRateModel::dryAdiabaticLapseRate() const {
    return AtmosphericThermodynamics::dryAdiabaticLapseRate(atmosphere_, surfaceGravity_);
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
