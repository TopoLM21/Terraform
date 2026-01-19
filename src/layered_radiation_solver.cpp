#include "layered_radiation_solver.h"

#include <QtGlobal>

#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kTwoStreamEddingtonFactor = 0.75;

// Эмиссивность слоя в LW диапазоне для двухпоточной модели Eddington.
// T_lw ≈ 1 / (1 + 3/4 * τ), ε = 1 - T_lw.
double layerEmissivity(double opticalDepthLongwave) {
    if (opticalDepthLongwave <= 0.0) {
        return 0.0;
    }
    const double transmission = 1.0 / (1.0 + kTwoStreamEddingtonFactor * opticalDepthLongwave);
    return qBound(0.0, 1.0 - transmission, 1.0);
}
} // namespace

LayeredRadiationSolver::LayeredRadiationSolver(double timeStepSeconds)
    : timeStepSeconds_(qMax(0.0, timeStepSeconds)) {}

QVector<double> LayeredRadiationSolver::solve(const AtmosphericColumn &column,
                                              double insolationWPerM2,
                                              double albedo,
                                              double cloudShortwaveTransmission) const {
    const auto &layers = column.layers();
    QVector<double> temperatureDeltas;
    temperatureDeltas.fill(0.0, layers.size());

    if (layers.isEmpty() || timeStepSeconds_ <= 0.0) {
        return temperatureDeltas;
    }

    const double surfaceAlbedo = qBound(0.0, albedo, 1.0);
    const double cloudTransmission = qBound(0.0, cloudShortwaveTransmission, 1.0);

    // Коротковолновое излучение идёт сверху вниз: S_{i+1} = S_i * exp(-τ_sw).
    // Поглощение слоя: ΔS = S_i - S_{i+1}.
    QVector<double> shortwaveAbsorbed;
    shortwaveAbsorbed.fill(0.0, layers.size());

    double shortwaveFlux = qMax(0.0, insolationWPerM2) * cloudTransmission;
    for (int i = layers.size() - 1; i >= 0; --i) {
        const double tauSw = qMax(0.0, layers.at(i).opticalDepthShortwave());
        const double transmitted = shortwaveFlux * std::exp(-tauSw);
        shortwaveAbsorbed[i] = shortwaveFlux - transmitted;
        shortwaveFlux = transmitted;
    }

    // Поток у поверхности после атмосферы (можно использовать для связки с поверхностью).
    const double surfaceShortwaveFlux = shortwaveFlux * (1.0 - surfaceAlbedo);
    Q_UNUSED(surfaceShortwaveFlux);

    // Длинноволновое излучение: двухпоточная модель Eddington по слоям.
    const double surfaceTemperature = qMax(0.0, layers.first().temperatureKelvin());
    const double surfaceEmission =
        kStefanBoltzmannConstant * std::pow(surfaceTemperature, 4.0);

    const int layerCount = layers.size();
    QVector<double> upwardFlux;
    QVector<double> downwardFlux;
    upwardFlux.fill(0.0, layerCount + 1);
    downwardFlux.fill(0.0, layerCount + 1);

    upwardFlux[0] = surfaceEmission;
    for (int i = 0; i < layerCount; ++i) {
        const double emissivity = layerEmissivity(layers.at(i).opticalDepthLongwave());
        const double emission =
            emissivity * kStefanBoltzmannConstant *
            std::pow(qMax(0.0, layers.at(i).temperatureKelvin()), 4.0);
        upwardFlux[i + 1] = upwardFlux[i] * (1.0 - emissivity) + emission;
    }

    downwardFlux[layerCount] = 0.0;
    for (int i = layerCount - 1; i >= 0; --i) {
        const double emissivity = layerEmissivity(layers.at(i).opticalDepthLongwave());
        const double emission =
            emissivity * kStefanBoltzmannConstant *
            std::pow(qMax(0.0, layers.at(i).temperatureKelvin()), 4.0);
        downwardFlux[i] = downwardFlux[i + 1] * (1.0 - emissivity) + emission;
    }

    for (int i = 0; i < layerCount; ++i) {
        const double heatCapacity = layers.at(i).heatCapacityJPerM2K();
        if (heatCapacity <= 0.0) {
            continue;
        }
        // Баланс слоя: (F_up, F_down) на границах слоя -> изменение энергии слоя.
        const double netLongwaveFlux =
            (upwardFlux[i] + downwardFlux[i + 1]) - (upwardFlux[i + 1] + downwardFlux[i]);
        const double netFlux = shortwaveAbsorbed[i] + netLongwaveFlux;
        temperatureDeltas[i] = netFlux * timeStepSeconds_ / heatCapacity;
    }

    return temperatureDeltas;
}
