#pragma once

#include "atmosphere_model.h"

// Модель барометрического распределения для расчёта давления на высоте.
class AtmosphericPressureModel {
public:
    // Барометрическая формула (идеальный газ):
    // p = p0 * exp(-z / H), где H = R_specific * T / g.
    // Единицы: p, p0 [атм], z [м], T [К], g [м/с^2].
    // Дополнительно: 1 атм = 101325 Па (используется для справочных переводов).
    static double pressureAtHeightAtm(double seaLevelPressureAtm,
                                      double heightMeters,
                                      double temperatureKelvin,
                                      const AtmosphereComposition &composition,
                                      double gravity);
};
