#pragma once

#include "planet_surface_point.h"

#include <QVector>

class VegetationModel {
public:
    struct Settings {
        // Максимальная биомасса (кг/м²) для нормировки доли покрытия.
        double maxBiomassKgPerM2 = 15.0;
        // Минимальная «семенная» биомасса, позволяющая росту стартовать (кг/м²).
        double biomassSeedKgPerM2 = 0.02;
        // Базовый темп роста (1/сутки).
        double growthRatePerDay = 0.08;
        // Базовый темп смертности (1/сутки).
        double mortalityRatePerDay = 0.02;
        // Коэффициент диффузии на соседние ячейки (1/сутки).
        double diffusionRatePerDay = 0.05;
        // Температурные пороги (K): минимум, оптимум, максимум.
        double minTempK = 260.0;
        double optimalTempK = 295.0;
        double maxTempK = 315.0;
        // Минимальная доля влажности почвы, при которой рост не нулевой.
        double minMoistureFraction = 0.15;
        // Минимальное давление (атм), ниже которого рост подавляется.
        double minPressureAtm = 0.2;
        // Минимальное парциальное давление CO₂ (атм) для фотосинтеза.
        double minCo2PartialAtm = 2.0e-4;
        // Поток света для насыщения (Вт/м²), после которого прирост почти не растёт.
        double lightSaturationWPerM2 = 400.0;
        // Усиление смертности при стрессе (когда климатический фактор близок к 0).
        double stressMortalityBoost = 1.5;
    };

    VegetationModel();
    explicit VegetationModel(const Settings &settings);

    void setSettings(const Settings &settings);
    const Settings &settings() const;

    // co2Share — доля CO₂ в атмосфере (по составу), используется вместе с
    // point.pressureAtm для расчёта парциального давления.
    // point.solarFluxWPerM2 задаёт освещённость для ограничения фотосинтеза.
    void update(QVector<SurfacePoint> &points,
                double timeStepSeconds,
                double co2Share) const;

private:
    double temperatureFactor(double temperatureK) const;
    double moistureFactor(double moistureFraction) const;
    double pressureFactor(double pressureAtm) const;
    double co2Factor(double co2PartialAtm) const;
    double lightFactor(double solarFluxWPerM2) const;

    Settings settings_;
};
