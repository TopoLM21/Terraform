#pragma once

#include "planet_surface_point.h"

class SnowModel {
public:
    struct Settings {
        // Температурный порог разделения дождя/снега и начала таяния.
        double freezingTemperatureK = 273.15;
        // Диапазон прогрева (К), в котором скорость таяния растёт от 0 до максимума.
        double meltTemperatureRangeK = 5.0;
        // Максимальная скорость таяния (кг/м² за сутки) при T >= T0 + range.
        double meltRateKgPerM2PerDay = 20.0;
    };

    SnowModel() = default;
    explicit SnowModel(const Settings &settings);

    const Settings &settings() const;
    void setSettings(const Settings &settings);

    // Обновляет снежный покров точки за шаг dtSeconds.
    void updatePoint(SurfacePoint &point, double dtSeconds) const;

private:
    Settings settings_;
};
