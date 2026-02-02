#pragma once

#include <cstddef>

// Упрощенная модель фазовой диаграммы воды с учетом массы и энергии.
// Компонент используется для оценки доминирующей фазы и расчета фазовых переходов.
class PhaseModel {
public:
    enum class Phase {
        Ice,
        Liquid,
        Vapor
    };

    struct PhaseState {
        double iceKg = 0.0;
        double liquidKg = 0.0;
        double vaporKg = 0.0;
        double temperatureK = 273.15;
    };

    struct AlbedoTable {
        double ice = 0.0;
        double liquid = 0.0;
        double vapor = 0.0;
    };

    PhaseModel() = default;

    // Определение доминирующей фазы по температуре и давлению.
    Phase phaseForTemperaturePressure(double temperatureK, double pressurePa) const;

    // Таблица альбедо для льда, жидкости и пара.
    static const AlbedoTable &albedoTable();
    static double albedoForPhase(Phase phase);

    // Добавляет/убирает энергию и корректирует фазы с учетом сохранения массы и энергии.
    void applyEnergyJ(PhaseState &state, double energyJ) const;

    double totalMassKg(const PhaseState &state) const;

private:
    static constexpr double kMeltingTempK = 273.15;
    static constexpr double kBoilingTempK = 373.15;
    static constexpr double kTriplePointTempK = 273.16;
    static constexpr double kTriplePointPressurePa = 611.657;
    static constexpr double kCriticalTempK = 647.096;
    static constexpr double kCriticalPressurePa = 22.064e6;

    static constexpr double kSpecificHeatIce = 2100.0;
    static constexpr double kSpecificHeatLiquid = 4182.0;
    static constexpr double kSpecificHeatVapor = 2000.0;
    static constexpr double kLatentHeatFusion = 334000.0;
    static constexpr double kLatentHeatVaporization = 2.5e6;
    static constexpr double kLatentHeatSublimation = 2.834e6;

    double boilingTemperatureK(double pressurePa) const;
    double heatCapacityJPerK(const PhaseState &state) const;
};
