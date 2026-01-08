#pragma once

class AtmosphericCellState {
public:
    AtmosphericCellState() = default;
    AtmosphericCellState(double airTemperatureKelvin, double heatCapacityJPerM2K);

    double airTemperatureKelvin() const;
    void setAirTemperatureKelvin(double temperatureKelvin);

    double heatCapacityJPerM2K() const;
    void setHeatCapacityJPerM2K(double heatCapacityJPerM2K);

private:
    double airTemperatureKelvin_ = 0.0;
    double heatCapacityJPerM2K_ = 0.0;
};
