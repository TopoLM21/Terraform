#pragma once

#include "atmospheric_column.h"

class EvaporationModel {
public:
    struct Settings {
        // Базовая относительная влажность внизу атмосферы (0..1).
        double baseRelativeHumidity = 0.65;
        // Характерное время (с) достижения фазового равновесия.
        double relaxationTimeSeconds = 600.0;
        // Максимальный вклад облаков в альбедо (0..1).
        double maxCloudAlbedo = 0.55;
        // Масштаб роста облачного альбедо по жидкой воде (кг/м²).
        double cloudAlbedoScaleKgPerM2 = 0.12;
        // Характерное время выпадения осадков (с).
        double precipitationTimeSeconds = 1800.0;
        // Доля конденсата (жидкость+лёд), уходящая в осадки за характерное время (0..1).
        double precipitationEfficiency = 0.6;
        // Температура замерзания воды (K), по которой делим конденсат на лёд и жидкость.
        double freezingTemperatureKelvin = 273.15;
        // Диапазон сглаживания (K) для смешанной фазы вокруг точки замерзания.
        double freezingRangeKelvin = 8.0;
        // Дополнительное "адиабатическое" охлаждение (K) для расчёта насыщения.
        // Используется как упрощённый сдвиг температуры перед фазовым балансом,
        // чтобы имитировать подъём влажного воздуха без явного расчёта энергии.
        // Не влияет на радиацию и термодинамику, поэтому по умолчанию отключено.
        double adiabaticCoolingDeltaKelvin = 0.0;
    };

    EvaporationModel() = default;
    explicit EvaporationModel(const Settings &settings);

    const Settings &settings() const;
    void setSettings(const Settings &settings);

    // Инициализировать слой по заданной относительной влажности.
    void initializeLayer(AtmosphericLayerState &layer, double relativeHumidity) const;

    // Обновить фазовый баланс (пар ↔ жидкость ↔ лёд) внутри колонны.
    void updateColumn(AtmosphericColumn &column, double timeStepSeconds) const;
    // Обновить фазовый баланс и вернуть выпавшие осадки (кг/м²), ушедшие на поверхность.
    double updateColumnWithPrecipitation(AtmosphericColumn &column,
                                         double timeStepSeconds) const;

    // Оценить вклад облаков в альбедо по суммарной массе конденсата.
    double cloudAlbedoFromCondensation(const AtmosphericColumn &column) const;

    // Давление насыщенного пара воды, Па.
    static double saturationVaporPressurePa(double temperatureKelvin);
    // Плотность насыщенного водяного пара, кг/м³.
    static double saturationVaporDensityKgPerM3(double temperatureKelvin);
    // Максимальная масса водяного пара в слое на единицу площади, кг/м².
    static double saturationWaterVaporKgPerM2(double temperatureKelvin, double thicknessMeters);

private:
    Settings settings_;
};
