#pragma once

#include "atmosphere_model.h"

// Утилиты для термодинамических свойств атмосферы.
class AtmosphericThermodynamics {
public:
    static double meanMolarMassKgPerMol(const AtmosphereComposition &composition);
    static double specificGasConstant(const AtmosphereComposition &composition);
    static double specificHeatCp(const AtmosphereComposition &composition);
    static double relativeHumidityEstimate(const AtmosphereComposition &composition);
    static double saturationVaporPressureAtm(double temperatureKelvin);
    static double moistAdiabaticLapseRate(double temperatureKelvin,
                                          const AtmosphereComposition &composition,
                                          double pressureAtm,
                                          double gravity);
    static double dryAdiabaticLapseRate(const AtmosphereComposition &composition,
                                        double gravity);
    static double adiabaticGradientDlnT(const AtmosphereComposition &composition,
                                        double pressureAtm,
                                        double temperatureKelvin,
                                        double gravity);
};
