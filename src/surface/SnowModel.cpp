#include "surface/SnowModel.h"

#include <QtCore/QString>
#include <QtCore/QtMath>

namespace {
constexpr double kSecondsPerDay = 86400.0;
}

SnowModel::SnowModel(const Settings &settings)
    : settings_(settings) {}

const SnowModel::Settings &SnowModel::settings() const {
    return settings_;
}

void SnowModel::setSettings(const Settings &settings) {
    settings_ = settings;
}

void SnowModel::updatePoint(SurfacePoint &point, double dtSeconds) const {
    if (dtSeconds <= 0.0) {
        return;
    }
    // Пропускаем жидкий океан — снег не ложится на воду.
    // Замёрзший океан (морской лёд) обрабатываем как сушу.
    if (point.materialId == QLatin1String("ocean") &&
        point.waterPhase == PhaseModel::Phase::Liquid) {
        return;
    }
    if (point.snowKgPerM2 <= 0.0) {
        return;
    }

    const double surfaceTemperatureK =
        (point.temperatureK > 0.0) ? point.temperatureK : point.airTemperatureK;
    if (surfaceTemperatureK <= settings_.freezingTemperatureK) {
        return;
    }

    const double temperatureExcess = surfaceTemperatureK - settings_.freezingTemperatureK;
    const double temperatureFactor = (settings_.meltTemperatureRangeK > 0.0)
        ? qBound(0.0, temperatureExcess / settings_.meltTemperatureRangeK, 1.0)
        : 1.0;
    // Упрощение: таяние линейно по температуре и ограничено суточной нормой.
    const double meltRateKgPerM2S =
        (settings_.meltRateKgPerM2PerDay > 0.0)
            ? settings_.meltRateKgPerM2PerDay / kSecondsPerDay
            : 0.0;
    const double potentialMeltKgPerM2 = meltRateKgPerM2S * temperatureFactor * dtSeconds;
    const double meltedKgPerM2 = qMin(point.snowKgPerM2, potentialMeltKgPerM2);
    if (meltedKgPerM2 <= 0.0) {
        return;
    }

    point.snowKgPerM2 = qMax(0.0, point.snowKgPerM2 - meltedKgPerM2);
    // Излишек, который не помещается в почвенную влагу, считаем поверхностным стоком.
    point.surfaceMoisture.addPrecipitation(meltedKgPerM2);
}
