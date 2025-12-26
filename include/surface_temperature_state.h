#pragma once

class SurfaceTemperatureState {
public:
    SurfaceTemperatureState(double initialTemperatureKelvin,
                            double albedo,
                            double heatCapacity,
                            double greenhouseOpacity,
                            double minTemperatureKelvin);

    double temperatureKelvin() const;
    void setTemperatureKelvin(double temperatureKelvin);

    void updateTemperature(double solarIrradiance, double dtSeconds);

private:
    double temperatureKelvin_ = 3.0;
    double albedo_ = 0.0;
    double heatCapacity_ = 1.0;
    double greenhouseOpacity_ = 0.0;
    double minTemperatureKelvin_ = 3.0;
};
