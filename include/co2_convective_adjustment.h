#pragma once

// Конвективная корректировка для CO2-атмосферы.
// Используется сухая адиабата:
//   Γ_d = g / c_p,
// а связь температуры с давлением (или оптической толщиной τ ~ p) задаёт
//   T(τ) = T_surface * (τ / τ_surface)^(R / c_p).

class Co2ConvectiveAdjustment {
public:
    static double dryAdiabaticLapseRate();
    static double adjustedEmissionTemperature(double surfaceTemperatureKelvin,
                                              double opticalDepth,
                                              double radiativeEmissionTemperatureKelvin);
};
