#pragma once

enum class ThermalConductivityModelType {
    Constant,
    Linear
};

struct ThermalConductivityModel {
    ThermalConductivityModelType type = ThermalConductivityModelType::Constant;
    // Опорная температура для линейной модели.
    double referenceTemperatureKelvin = 250.0;
    // Линейный коэффициент dK/dT в Вт/(м·К²):
    // k(T) = k0 + coefficient * (T - T_ref).
    double temperatureCoefficient = 0.0;
};

double thermalConductivityAtTemperature(double baseConductivity,
                                        const ThermalConductivityModel &model,
                                        double temperatureKelvin);
