#include "atmospheric_cell_state.h"

AtmosphericCellState::AtmosphericCellState(double airTemperatureKelvin,
                                           double heatCapacityJPerM2K)
    : airTemperatureKelvin_(airTemperatureKelvin),
      heatCapacityJPerM2K_(heatCapacityJPerM2K) {}

double AtmosphericCellState::airTemperatureKelvin() const {
    return airTemperatureKelvin_;
}

void AtmosphericCellState::setAirTemperatureKelvin(double temperatureKelvin) {
    airTemperatureKelvin_ = temperatureKelvin;
}

double AtmosphericCellState::heatCapacityJPerM2K() const {
    return heatCapacityJPerM2K_;
}

void AtmosphericCellState::setHeatCapacityJPerM2K(double heatCapacityJPerM2K) {
    heatCapacityJPerM2K_ = heatCapacityJPerM2K;
}
