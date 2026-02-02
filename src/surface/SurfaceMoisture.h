#pragma once

class SurfaceMoisture {
public:
    struct Settings {
        // Максимальная доступная влага в верхнем слое поверхности, кг/м².
        double maxStorageKgPerM2 = 200.0;
    };

    SurfaceMoisture() = default;
    explicit SurfaceMoisture(const Settings &settings);

    const Settings &settings() const;
    void setSettings(const Settings &settings);

    double waterKgPerM2() const;
    void setWaterKgPerM2(double waterKgPerM2);

    // Добавить осадки (кг/м²) и вернуть объём стока (избыточная вода).
    double addPrecipitation(double precipitationKgPerM2);
    // Применить испарение (кг/м²), возвращает фактически испарённую массу.
    double applyEvaporation(double potentialEvaporationKgPerM2);
    // Обновить влагу поверхности за шаг, учитывая осадки и испарение.
    double update(double precipitationKgPerM2, double potentialEvaporationKgPerM2);

    // Доля заполнения слоя влагой (0..1), пригодна для будущих биомоделей.
    double moistureFraction() const;
    // Запросить воду для биосферы, возвращает фактически выделенную массу.
    double takeWaterKgPerM2(double requestedKgPerM2);
    // Вернуть воду обратно на поверхность (например, после транспирации).
    void returnWaterKgPerM2(double returnedKgPerM2);

private:
    double clampedWater(double waterKgPerM2) const;

    Settings settings_;
    double waterKgPerM2_ = 0.0;
};
