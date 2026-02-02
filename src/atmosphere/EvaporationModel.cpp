#include "atmosphere/EvaporationModel.h"

#include <cmath>
#include <QtMath>

namespace {
constexpr double kWaterVaporGasConstant = 461.5; // Дж/(кг·К).

// Линейная интерполяция для плавной релаксации.
double lerp(double a, double b, double t) {
    return a + (b - a) * t;
}

double iceFractionFromTemperature(double temperatureKelvin,
                                  double freezingTemperatureKelvin,
                                  double freezingRangeKelvin) {
    if (freezingRangeKelvin <= 0.0) {
        return (temperatureKelvin < freezingTemperatureKelvin) ? 1.0 : 0.0;
    }
    // Допущение: смешанная фаза линейно распределяется в диапазоне
    // [T_freeze - dT, T_freeze], чтобы избежать резких скачков конденсата.
    const double fraction =
        (freezingTemperatureKelvin - temperatureKelvin) / freezingRangeKelvin;
    return qBound(0.0, fraction, 1.0);
}
} // namespace

EvaporationModel::EvaporationModel(const Settings &settings)
    : settings_(settings) {}

const EvaporationModel::Settings &EvaporationModel::settings() const {
    return settings_;
}

void EvaporationModel::setSettings(const Settings &settings) {
    settings_ = settings;
}

void EvaporationModel::initializeLayer(AtmosphericLayerState &layer,
                                       double relativeHumidity) const {
    const double clampedHumidity = qBound(0.0, relativeHumidity, 1.0);
    const double maxVaporKgPerM2 =
        saturationWaterVaporKgPerM2(layer.temperatureKelvin(), layer.thicknessMeters());
    const double vaporKgPerM2 = maxVaporKgPerM2 * clampedHumidity;

    layer.setWaterVaporKgPerM2(vaporKgPerM2);
    layer.setLiquidWaterKgPerM2(0.0);
    layer.setIceWaterKgPerM2(0.0);
    layer.setRelativeHumidity(clampedHumidity);
}

void EvaporationModel::updateColumn(AtmosphericColumn &column,
                                    double timeStepSeconds) const {
    updateColumnWithPrecipitation(column, timeStepSeconds);
}

double EvaporationModel::updateColumnWithPrecipitation(AtmosphericColumn &column,
                                                       double timeStepSeconds) const {
    auto &layers = column.layers();
    if (layers.isEmpty()) {
        return 0.0;
    }

    const double relaxation = (settings_.relaxationTimeSeconds > 0.0 && timeStepSeconds > 0.0)
        ? qBound(0.0, timeStepSeconds / settings_.relaxationTimeSeconds, 1.0)
        : 1.0;
    const double precipitationRelaxation =
        (settings_.precipitationTimeSeconds > 0.0 && timeStepSeconds > 0.0)
            ? qBound(0.0, timeStepSeconds / settings_.precipitationTimeSeconds, 1.0)
            : 1.0;

    double precipitationKgPerM2 = 0.0;

    for (auto &layer : layers) {
        const double maxVaporKgPerM2 =
            saturationWaterVaporKgPerM2(layer.temperatureKelvin(), layer.thicknessMeters());
        const double vaporKgPerM2 = qMax(0.0, layer.waterVaporKgPerM2());
        const double liquidKgPerM2 = qMax(0.0, layer.liquidWaterKgPerM2());
        const double iceKgPerM2 = qMax(0.0, layer.iceWaterKgPerM2());
        const double totalWaterKgPerM2 = vaporKgPerM2 + liquidKgPerM2 + iceKgPerM2;

        if (maxVaporKgPerM2 <= 0.0 || totalWaterKgPerM2 <= 0.0) {
            layer.setWaterVaporKgPerM2(0.0);
            layer.setLiquidWaterKgPerM2(0.0);
            layer.setIceWaterKgPerM2(0.0);
            layer.setRelativeHumidity(0.0);
            continue;
        }

        // Целевое равновесие: пар ограничен насыщением, остальное — конденсат.
        // Масса сконденсированной влаги = max(0, W_total - W_sat).
        const double targetVaporKgPerM2 = qMin(totalWaterKgPerM2, maxVaporKgPerM2);
        const double updatedVaporKgPerM2 =
            lerp(vaporKgPerM2, targetVaporKgPerM2, relaxation);
        const double updatedCondensateKgPerM2 =
            qMax(0.0, totalWaterKgPerM2 - updatedVaporKgPerM2);
        const double iceFraction =
            iceFractionFromTemperature(layer.temperatureKelvin(),
                                       settings_.freezingTemperatureKelvin,
                                       settings_.freezingRangeKelvin);
        double updatedIceKgPerM2 = updatedCondensateKgPerM2 * iceFraction;
        double updatedLiquidKgPerM2 = updatedCondensateKgPerM2 - updatedIceKgPerM2;

        if (updatedCondensateKgPerM2 > 0.0 && settings_.precipitationEfficiency > 0.0) {
            // Осадки как релаксация: P = C * eff * dt / tau_precip,
            // где C — суммарная масса конденсата (жидкость+лёд) в слое (кг/м²).
            const double precipitationShare =
                qBound(0.0, settings_.precipitationEfficiency, 1.0) *
                precipitationRelaxation;
            const double layerPrecipitation =
                qBound(0.0,
                       updatedCondensateKgPerM2 * precipitationShare,
                       updatedCondensateKgPerM2);
            const double remainingCondensate =
                qMax(0.0, updatedCondensateKgPerM2 - layerPrecipitation);
            const double remainingFraction = (updatedCondensateKgPerM2 > 0.0)
                ? (remainingCondensate / updatedCondensateKgPerM2)
                : 0.0;
            updatedLiquidKgPerM2 *= remainingFraction;
            updatedIceKgPerM2 *= remainingFraction;
            precipitationKgPerM2 += layerPrecipitation;
        }

        layer.setWaterVaporKgPerM2(updatedVaporKgPerM2);
        layer.setLiquidWaterKgPerM2(updatedLiquidKgPerM2);
        layer.setIceWaterKgPerM2(updatedIceKgPerM2);
        layer.setRelativeHumidity(
            (maxVaporKgPerM2 > 0.0)
                ? qBound(0.0, updatedVaporKgPerM2 / maxVaporKgPerM2, 1.0)
                : 0.0);
    }

    return precipitationKgPerM2;
}

double EvaporationModel::cloudAlbedoFromCondensation(const AtmosphericColumn &column) const {
    if (settings_.maxCloudAlbedo <= 0.0 || settings_.cloudAlbedoScaleKgPerM2 <= 0.0) {
        return 0.0;
    }

    double totalCondensedWaterPath = 0.0;
    const auto &layers = column.layers();
    for (const auto &layer : layers) {
        totalCondensedWaterPath +=
            qMax(0.0, layer.liquidWaterKgPerM2()) + qMax(0.0, layer.iceWaterKgPerM2());
    }

    if (totalCondensedWaterPath <= 0.0) {
        return 0.0;
    }

    // Альбедо облаков растёт по экспоненте от массы конденсата (LWP+IWP),
    // игнорируя различия в микрофизике — это упрощение для визуализации.
    const double albedo = settings_.maxCloudAlbedo *
        (1.0 - std::exp(-totalCondensedWaterPath / settings_.cloudAlbedoScaleKgPerM2));
    return qBound(0.0, albedo, settings_.maxCloudAlbedo);
}

double EvaporationModel::saturationVaporPressurePa(double temperatureKelvin) {
    if (temperatureKelvin <= 0.0) {
        return 0.0;
    }

    // Аппроксимация Тетенса для воды: e_s = 610.94 * exp(17.625 * T / (T + 243.04)).
    // Формула даёт давление насыщения в зависимости от температуры и описывает,
    // сколько пара может удерживать воздух без конденсации.
    // T в °C, диапазон подходит для тропосферы и визуализации облачности.
    const double temperatureC = temperatureKelvin - 273.15;
    const double exponent = 17.625 * temperatureC / (temperatureC + 243.04);
    return 610.94 * std::exp(exponent);
}

double EvaporationModel::saturationVaporDensityKgPerM3(double temperatureKelvin) {
    if (temperatureKelvin <= 0.0) {
        return 0.0;
    }
    const double pressurePa = saturationVaporPressurePa(temperatureKelvin);
    // Идеальный газ: rho = p / (R_v * T).
    return (pressurePa > 0.0) ? pressurePa / (kWaterVaporGasConstant * temperatureKelvin) : 0.0;
}

double EvaporationModel::saturationWaterVaporKgPerM2(double temperatureKelvin,
                                                     double thicknessMeters) {
    if (thicknessMeters <= 0.0) {
        return 0.0;
    }
    return saturationVaporDensityKgPerM3(temperatureKelvin) * thicknessMeters;
}
