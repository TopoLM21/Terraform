#pragma once

// Утилиты для пересчёта периодов вращения.
// Формула связывает солнечный (synodic) и сидерический периоды:
// 1 / P_solar = 1 / P_sidereal - 1 / P_orbit.
inline double solarToSiderealPeriodDays(double solarPeriodDays, double orbitalPeriodDays) {
    if (solarPeriodDays <= 0.0 || orbitalPeriodDays <= 0.0) {
        return 0.0;
    }
    return 1.0 / (1.0 / solarPeriodDays + 1.0 / orbitalPeriodDays);
}
