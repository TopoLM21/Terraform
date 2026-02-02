#include "PhaseModel.h"

#include <algorithm>
#include <cmath>

PhaseModel::Phase PhaseModel::phaseForTemperaturePressure(double temperatureK,
                                                          double pressurePa) const {
    if (temperatureK >= kCriticalTempK && pressurePa >= kCriticalPressurePa) {
        return Phase::Vapor;
    }

    if (pressurePa <= kTriplePointPressurePa) {
        return (temperatureK < kTriplePointTempK) ? Phase::Ice : Phase::Vapor;
    }

    if (temperatureK < kMeltingTempK) {
        return Phase::Ice;
    }

    const double boilingK = boilingTemperatureK(pressurePa);
    return (temperatureK >= boilingK) ? Phase::Vapor : Phase::Liquid;
}

const PhaseModel::AlbedoTable &PhaseModel::albedoTable() {
    static const AlbedoTable table{
        0.6,  // Лед: высокое альбедо свежего льда/снега в видимом диапазоне.
        0.07, // Жидкая вода: типичные океанические значения при умеренных углах.
        0.02  // Пар: низкое альбедо чистого водяного пара без облачности.
    };
    return table;
}

double PhaseModel::albedoForPhase(PhaseModel::Phase phase) {
    const AlbedoTable &table = albedoTable();
    switch (phase) {
    case Phase::Ice:
        return table.ice;
    case Phase::Liquid:
        return table.liquid;
    case Phase::Vapor:
        return table.vapor;
    }
    return table.liquid;
}

void PhaseModel::applyEnergyJ(PhaseState &state, double energyJ) const {
    // Алгоритм двигает температуру к фазовым порогам, затем тратит энергию
    // на скрытую теплоту переходов. Масса не создается и не исчезает, а энергия
    // распределяется между нагревом и фазовыми превращениями.
    double energyRemaining = energyJ;
    if (energyRemaining == 0.0) {
        return;
    }

    if (energyRemaining > 0.0) {
        if (state.temperatureK < kMeltingTempK) {
            const double heatCapacity = heatCapacityJPerK(state);
            if (heatCapacity > 0.0) {
                const double energyToMeltTemp = (kMeltingTempK - state.temperatureK) * heatCapacity;
                if (energyRemaining < energyToMeltTemp) {
                    state.temperatureK += energyRemaining / heatCapacity;
                    return;
                }
                state.temperatureK = kMeltingTempK;
                energyRemaining -= energyToMeltTemp;
            }
        }

        if (energyRemaining > 0.0 && state.iceKg > 0.0) {
            const double meltMass = std::min(state.iceKg, energyRemaining / kLatentHeatFusion);
            state.iceKg -= meltMass;
            state.liquidKg += meltMass;
            energyRemaining -= meltMass * kLatentHeatFusion;
        }

        if (energyRemaining > 0.0 && state.temperatureK < kBoilingTempK) {
            const double heatCapacity = heatCapacityJPerK(state);
            if (heatCapacity > 0.0) {
                const double energyToBoilTemp = (kBoilingTempK - state.temperatureK) * heatCapacity;
                if (energyRemaining < energyToBoilTemp) {
                    state.temperatureK += energyRemaining / heatCapacity;
                    return;
                }
                state.temperatureK = kBoilingTempK;
                energyRemaining -= energyToBoilTemp;
            }
        }

        if (energyRemaining > 0.0 && state.liquidKg > 0.0) {
            const double vaporMass =
                std::min(state.liquidKg, energyRemaining / kLatentHeatVaporization);
            state.liquidKg -= vaporMass;
            state.vaporKg += vaporMass;
            energyRemaining -= vaporMass * kLatentHeatVaporization;
        }

        if (energyRemaining > 0.0) {
            const double heatCapacity = heatCapacityJPerK(state);
            if (heatCapacity > 0.0) {
                state.temperatureK += energyRemaining / heatCapacity;
            }
        }
        return;
    }

    energyRemaining = -energyRemaining;
    if (state.temperatureK > kBoilingTempK) {
        const double heatCapacity = heatCapacityJPerK(state);
        if (heatCapacity > 0.0) {
            const double energyToCoolTemp = (state.temperatureK - kBoilingTempK) * heatCapacity;
            if (energyRemaining < energyToCoolTemp) {
                state.temperatureK -= energyRemaining / heatCapacity;
                return;
            }
            state.temperatureK = kBoilingTempK;
            energyRemaining -= energyToCoolTemp;
        }
    }

    if (energyRemaining > 0.0 && state.vaporKg > 0.0) {
        const double condenseMass =
            std::min(state.vaporKg, energyRemaining / kLatentHeatVaporization);
        state.vaporKg -= condenseMass;
        state.liquidKg += condenseMass;
        energyRemaining -= condenseMass * kLatentHeatVaporization;
    }

    if (energyRemaining > 0.0 && state.temperatureK > kMeltingTempK) {
        const double heatCapacity = heatCapacityJPerK(state);
        if (heatCapacity > 0.0) {
            const double energyToCoolTemp = (state.temperatureK - kMeltingTempK) * heatCapacity;
            if (energyRemaining < energyToCoolTemp) {
                state.temperatureK -= energyRemaining / heatCapacity;
                return;
            }
            state.temperatureK = kMeltingTempK;
            energyRemaining -= energyToCoolTemp;
        }
    }

    if (energyRemaining > 0.0 && state.liquidKg > 0.0) {
        const double freezeMass = std::min(state.liquidKg, energyRemaining / kLatentHeatFusion);
        state.liquidKg -= freezeMass;
        state.iceKg += freezeMass;
        energyRemaining -= freezeMass * kLatentHeatFusion;
    }

    if (energyRemaining > 0.0) {
        const double heatCapacity = heatCapacityJPerK(state);
        if (heatCapacity > 0.0) {
            state.temperatureK -= energyRemaining / heatCapacity;
        }
    }
}

double PhaseModel::totalMassKg(const PhaseModel::PhaseState &state) const {
    return state.iceKg + state.liquidKg + state.vaporKg;
}

double PhaseModel::boilingTemperatureK(double pressurePa) const {
    if (pressurePa <= 0.0) {
        return kBoilingTempK;
    }

    // Упрощенная аппроксимация по уравнению Антуана для диапазона 1-10^5 Па.
    const double pressureMmHg = pressurePa / 133.322;
    const double a = 8.07131;
    const double b = 1730.63;
    const double c = 233.426;
    const double temperatureC = (b / (a - std::log10(pressureMmHg))) - c;
    return std::clamp(temperatureC + 273.15, kTriplePointTempK, kCriticalTempK);
}

double PhaseModel::heatCapacityJPerK(const PhaseModel::PhaseState &state) const {
    const double capacity =
        state.iceKg * kSpecificHeatIce + state.liquidKg * kSpecificHeatLiquid +
        state.vaporKg * kSpecificHeatVapor;
    return capacity;
}
