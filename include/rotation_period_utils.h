#pragma once

// Утилиты для пересчёта периодов вращения.
// Формула связывает солнечный (synodic) и сидерический периоды.
// Проградное вращение: 1 / P_solar = 1 / P_sidereal - 1 / P_orbit.
// Ретроградное вращение: 1 / P_solar = 1 / P_sidereal + 1 / P_orbit.
inline double solarToSiderealPeriodDays(double solarPeriodDays,
                                        double orbitalPeriodDays,
                                        bool isRetrograde) {
    if (solarPeriodDays <= 0.0 || orbitalPeriodDays <= 0.0) {
        return 0.0;
    }
    // Для ретроградного вращения период сидерического вращения длиннее,
    // поэтому знак орбитальной частоты меняется в формуле пересчёта.
    const double orbitalContribution = isRetrograde ? -1.0 : 1.0;
    const double siderealFrequency =
        1.0 / solarPeriodDays + orbitalContribution / orbitalPeriodDays;
    if (siderealFrequency <= 0.0) {
        return 0.0;
    }
    return 1.0 / siderealFrequency;
}

inline bool isRetrogradeRotation(double obliquityDegrees) {
    // Простая эвристика: наклон оси > 90° соответствует ретроградному вращению.
    return obliquityDegrees > 90.0;
}
