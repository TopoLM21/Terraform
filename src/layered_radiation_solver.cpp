#include "layered_radiation_solver.h"

#include <QtGlobal>

#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kTwoStreamEddingtonFactor = 0.75;
constexpr double kMinHeatCapacityJPerM2K = 1.0e-3;
constexpr double kMaxDeltaKPerStep = 20.0;

// Доли Планка для каждой ИК-полосы при ~280 K.
// Far-IR (>17мкм) 30%, Window (8-13мкм) 38%, CO₂ (13-17мкм) 18%, H₂O+CH₄ (5-8мкм) 14%.
constexpr double kPlanckFraction[kRadiationBandCount] = {0.30, 0.38, 0.18, 0.14};

// Эмиссивность слоя в LW диапазоне для двухпоточной модели Eddington.
// T_lw ≈ 1 / (1 + 3/4 * τ), ε = 1 - T_lw.
double layerEmissivity(double opticalDepthLongwave) {
    if (opticalDepthLongwave <= 0.0) {
        return 0.0;
    }
    const double transmission = 1.0 / (1.0 + kTwoStreamEddingtonFactor * opticalDepthLongwave);
    return qBound(0.0, 1.0 - transmission, 1.0);
}

// Безопасное излучение чёрного тела (σT⁴).
double safeBlackbodyEmission(double temperatureKelvin) {
    if (!std::isfinite(temperatureKelvin) || temperatureKelvin <= 0.0) return 0.0;
    const double e = kStefanBoltzmannConstant * std::pow(temperatureKelvin, 4.0);
    return std::isfinite(e) ? e : 0.0;
}
} // namespace

LayeredRadiationSolver::LayeredRadiationSolver(double timeStepSeconds)
    : timeStepSeconds_(qMax(0.0, timeStepSeconds)) {}

QVector<double> LayeredRadiationSolver::solve(const AtmosphericColumn &column,
                                              double insolationWPerM2,
                                              double albedo,
                                              double cloudShortwaveTransmission,
                                              double surfaceTemperatureKelvin) const {
    SurfaceRadiativeFluxes unused;
    return solve(column, insolationWPerM2, albedo, cloudShortwaveTransmission,
                 surfaceTemperatureKelvin, unused);
}

QVector<double> LayeredRadiationSolver::solve(const AtmosphericColumn &column,
                                              double insolationWPerM2,
                                              double albedo,
                                              double cloudShortwaveTransmission,
                                              double surfaceTemperatureKelvin,
                                              SurfaceRadiativeFluxes &surfaceFluxes) const {
    const auto &layers = column.layers();
    QVector<double> temperatureDeltas;
    temperatureDeltas.fill(0.0, layers.size());
    surfaceFluxes = {};

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

    // Поток у поверхности после прохождения атмосферы: SW_surface = SW_bottom * (1 - albedo).
    surfaceFluxes.shortwaveAbsorbedWPerM2 = shortwaveFlux * (1.0 - surfaceAlbedo);

    // Длинноволновое излучение: двухпоточная модель Eddington по слоям.
    double surfaceTemperature = surfaceTemperatureKelvin;
    if (!std::isfinite(surfaceTemperature)) {
        // Защита от переполнений из-за степенной зависимости T^4.
        surfaceTemperature = 0.0;
    }
    surfaceTemperature = qMax(0.0, surfaceTemperature);
    double surfaceEmission = kStefanBoltzmannConstant * std::pow(surfaceTemperature, 4.0);
    if (!std::isfinite(surfaceEmission)) {
        // Защита от переполнений из-за степенной зависимости T^4.
        surfaceEmission = 0.0;
    }

    const int layerCount = layers.size();
    QVector<double> upwardFlux;
    QVector<double> downwardFlux;
    upwardFlux.fill(0.0, layerCount + 1);
    downwardFlux.fill(0.0, layerCount + 1);

    upwardFlux[0] = surfaceEmission;
    for (int i = 0; i < layerCount; ++i) {
        const double emissivity = layerEmissivity(layers.at(i).opticalDepthLongwave());
        double layerTemperature = layers.at(i).temperatureKelvin();
        if (!std::isfinite(layerTemperature)) {
            // Защита от переполнений из-за степенной зависимости T^4.
            layerTemperature = 0.0;
        }
        layerTemperature = qMax(0.0, layerTemperature);
        double emission = emissivity * kStefanBoltzmannConstant * std::pow(layerTemperature, 4.0);
        if (!std::isfinite(emission)) {
            // Защита от переполнений из-за степенной зависимости T^4.
            emission = 0.0;
        }
        upwardFlux[i + 1] = upwardFlux[i] * (1.0 - emissivity) + emission;
        if (!std::isfinite(upwardFlux[i + 1])) {
            // Защита от переполнений из-за степенной зависимости T^4.
            upwardFlux[i + 1] = 0.0;
        }
    }

    downwardFlux[layerCount] = 0.0;
    for (int i = layerCount - 1; i >= 0; --i) {
        const double emissivity = layerEmissivity(layers.at(i).opticalDepthLongwave());
        double layerTemperature = layers.at(i).temperatureKelvin();
        if (!std::isfinite(layerTemperature)) {
            // Защита от переполнений из-за степенной зависимости T^4.
            layerTemperature = 0.0;
        }
        layerTemperature = qMax(0.0, layerTemperature);
        double emission = emissivity * kStefanBoltzmannConstant * std::pow(layerTemperature, 4.0);
        if (!std::isfinite(emission)) {
            // Защита от переполнений из-за степенной зависимости T^4.
            emission = 0.0;
        }
        downwardFlux[i] = downwardFlux[i + 1] * (1.0 - emissivity) + emission;
        if (!std::isfinite(downwardFlux[i])) {
            // Защита от переполнений из-за степенной зависимости T^4.
            downwardFlux[i] = 0.0;
        }
    }

    // Нисходящий длинноволновый поток от атмосферы к поверхности (atmospheric backradiation).
    // Это ключевой компонент парникового эффекта: атмосфера переизлучает ИК вниз.
    surfaceFluxes.longwaveDownWPerM2 = downwardFlux[0];

    for (int i = 0; i < layerCount; ++i) {
        const double heatCapacity = layers.at(i).heatCapacityJPerM2K();
        if (heatCapacity <= 0.0) {
            continue;
        }
        // Баланс слоя: (F_up, F_down) на границах слоя -> изменение энергии слоя.
        double netLongwaveFlux =
            (upwardFlux[i] + downwardFlux[i + 1]) - (upwardFlux[i + 1] + downwardFlux[i]);
        if (!std::isfinite(netLongwaveFlux)) {
            // Защита от переполнений из-за степенной зависимости T^4.
            netLongwaveFlux = 0.0;
        }
        const double netFlux = shortwaveAbsorbed[i] + netLongwaveFlux;
        const double rawDelta = netFlux * timeStepSeconds_ / heatCapacity;
        // Ограничиваем шаг температуры для численной устойчивости в разрежённых слоях
        // (особенно у верхней границы атмосферы с малой теплоёмкостью), при этом
        // тонкие слои всё равно должны откликаться на радиационные потоки.
        temperatureDeltas[i] = qBound(-kMaxDeltaKPerStep, rawDelta, kMaxDeltaKPerStep);
    }

    return temperatureDeltas;
}

QVector<double> LayeredRadiationSolver::solveMultiband(
    const AtmosphericColumn &column,
    const QVector<BandOpticalDepths> &bandTaus,
    double insolationWPerM2,
    double albedo,
    double cloudShortwaveTransmission,
    double surfaceTemperatureKelvin,
    SurfaceRadiativeFluxes &surfaceFluxes) const {
    const auto &layers = column.layers();
    const int layerCount = layers.size();
    QVector<double> temperatureDeltas;
    temperatureDeltas.fill(0.0, layerCount);
    surfaceFluxes = {};

    if (layerCount == 0 || timeStepSeconds_ <= 0.0 ||
        bandTaus.size() != layerCount) {
        return temperatureDeltas;
    }

    const double surfaceAlbedo = qBound(0.0, albedo, 1.0);
    const double cloudTransmission = qBound(0.0, cloudShortwaveTransmission, 1.0);

    // ── Коротковолновое: один проход (солнечный спектр, без разбивки по ИК-полосам) ──
    QVector<double> shortwaveAbsorbed;
    shortwaveAbsorbed.fill(0.0, layerCount);
    double shortwaveFlux = qMax(0.0, insolationWPerM2) * cloudTransmission;
    for (int i = layerCount - 1; i >= 0; --i) {
        const double tauSw = qMax(0.0, layers.at(i).opticalDepthShortwave());
        const double transmitted = shortwaveFlux * std::exp(-tauSw);
        shortwaveAbsorbed[i] = shortwaveFlux - transmitted;
        shortwaveFlux = transmitted;
    }
    surfaceFluxes.shortwaveAbsorbedWPerM2 = shortwaveFlux * (1.0 - surfaceAlbedo);

    // ── Длинноволновое: независимый двухпоточный перенос для каждой ИК-полосы ──
    const double surfaceEmission = safeBlackbodyEmission(surfaceTemperatureKelvin);

    // Предвычисляем σT⁴ каждого слоя один раз.
    QVector<double> layerEmissions(layerCount);
    for (int i = 0; i < layerCount; ++i) {
        layerEmissions[i] = safeBlackbodyEmission(layers.at(i).temperatureKelvin());
    }

    // Суммарные потоки на границах (layerCount+1 граница).
    QVector<double> totalUpward(layerCount + 1, 0.0);
    QVector<double> totalDownward(layerCount + 1, 0.0);

    for (int b = 0; b < kRadiationBandCount; ++b) {
        const double frac = kPlanckFraction[b];

        // Восходящий проход: от поверхности (граница 0) вверх.
        QVector<double> bandUp(layerCount + 1);
        bandUp[0] = frac * surfaceEmission;
        for (int i = 0; i < layerCount; ++i) {
            const double eps = layerEmissivity(bandTaus.at(i).tauLw[b]);
            const double emission = frac * eps * layerEmissions.at(i);
            double up = bandUp.at(i) * (1.0 - eps) + emission;
            if (!std::isfinite(up)) up = 0.0;
            bandUp[i + 1] = up;
        }

        // Нисходящий проход: от верхней границы вниз.
        QVector<double> bandDown(layerCount + 1);
        bandDown[layerCount] = 0.0;
        for (int i = layerCount - 1; i >= 0; --i) {
            const double eps = layerEmissivity(bandTaus.at(i).tauLw[b]);
            const double emission = frac * eps * layerEmissions.at(i);
            double down = bandDown.at(i + 1) * (1.0 - eps) + emission;
            if (!std::isfinite(down)) down = 0.0;
            bandDown[i] = down;
        }

        // Суммируем по полосам.
        for (int j = 0; j <= layerCount; ++j) {
            totalUpward[j] += bandUp.at(j);
            totalDownward[j] += bandDown.at(j);
        }
    }

    // Нисходящий LW к поверхности = backradiation (сумма по всем полосам).
    surfaceFluxes.longwaveDownWPerM2 = totalDownward[0];

    // ── Температурные дельты ─────────────────────────────────────────────
    for (int i = 0; i < layerCount; ++i) {
        const double heatCapacity = layers.at(i).heatCapacityJPerM2K();
        if (heatCapacity <= 0.0) continue;
        double netLongwaveFlux =
            (totalUpward[i] + totalDownward[i + 1]) -
            (totalUpward[i + 1] + totalDownward[i]);
        if (!std::isfinite(netLongwaveFlux)) netLongwaveFlux = 0.0;
        const double netFlux = shortwaveAbsorbed[i] + netLongwaveFlux;
        const double rawDelta = netFlux * timeStepSeconds_ / heatCapacity;
        temperatureDeltas[i] = qBound(-kMaxDeltaKPerStep, rawDelta, kMaxDeltaKPerStep);
    }

    return temperatureDeltas;
}
