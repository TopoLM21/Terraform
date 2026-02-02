#include "atmosphere/EvaporationModel.h"

#include <cmath>
#include <QtMath>

namespace {
constexpr double kWaterVaporGasConstant = 461.5; // Дж/(кг·К).

// Линейная интерполяция для плавной релаксации.
double lerp(double a, double b, double t) {
    return a + (b - a) * t;
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
    layer.setRelativeHumidity(clampedHumidity);
}

void EvaporationModel::updateColumn(AtmosphericColumn &column,
                                    double timeStepSeconds) const {
    auto &layers = column.layers();
    if (layers.isEmpty()) {
        return;
    }

    const double relaxation = (settings_.relaxationTimeSeconds > 0.0 && timeStepSeconds > 0.0)
        ? qBound(0.0, timeStepSeconds / settings_.relaxationTimeSeconds, 1.0)
        : 1.0;

    for (auto &layer : layers) {
        const double maxVaporKgPerM2 =
            saturationWaterVaporKgPerM2(layer.temperatureKelvin(), layer.thicknessMeters());
        const double vaporKgPerM2 = qMax(0.0, layer.waterVaporKgPerM2());
        const double liquidKgPerM2 = qMax(0.0, layer.liquidWaterKgPerM2());
        const double totalWaterKgPerM2 = vaporKgPerM2 + liquidKgPerM2;

        if (maxVaporKgPerM2 <= 0.0 || totalWaterKgPerM2 <= 0.0) {
            layer.setWaterVaporKgPerM2(0.0);
            layer.setLiquidWaterKgPerM2(0.0);
            layer.setRelativeHumidity(0.0);
            continue;
        }

        // Целевое равновесие: пар ограничен насыщением, остальное — конденсат.
        const double targetVaporKgPerM2 = qMin(totalWaterKgPerM2, maxVaporKgPerM2);
        const double updatedVaporKgPerM2 =
            lerp(vaporKgPerM2, targetVaporKgPerM2, relaxation);
        const double updatedLiquidKgPerM2 =
            qMax(0.0, totalWaterKgPerM2 - updatedVaporKgPerM2);

        layer.setWaterVaporKgPerM2(updatedVaporKgPerM2);
        layer.setLiquidWaterKgPerM2(updatedLiquidKgPerM2);
        layer.setRelativeHumidity(
            (maxVaporKgPerM2 > 0.0)
                ? qBound(0.0, updatedVaporKgPerM2 / maxVaporKgPerM2, 1.0)
                : 0.0);
    }
}

double EvaporationModel::cloudAlbedoFromCondensation(const AtmosphericColumn &column) const {
    if (settings_.maxCloudAlbedo <= 0.0 || settings_.cloudAlbedoScaleKgPerM2 <= 0.0) {
        return 0.0;
    }

    double totalLiquidWaterPath = 0.0;
    const auto &layers = column.layers();
    for (const auto &layer : layers) {
        totalLiquidWaterPath += qMax(0.0, layer.liquidWaterKgPerM2());
    }

    if (totalLiquidWaterPath <= 0.0) {
        return 0.0;
    }

    // Альбедо облаков растёт по экспоненте от жидкой водяной массы (LWP).
    const double albedo = settings_.maxCloudAlbedo *
        (1.0 - std::exp(-totalLiquidWaterPath / settings_.cloudAlbedoScaleKgPerM2));
    return qBound(0.0, albedo, settings_.maxCloudAlbedo);
}

double EvaporationModel::saturationVaporPressurePa(double temperatureKelvin) {
    if (temperatureKelvin <= 0.0) {
        return 0.0;
    }

    // Аппроксимация Тетенса для воды: e_s = 610.94 * exp(17.625 * T / (T + 243.04)).
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
