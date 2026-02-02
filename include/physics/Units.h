#pragma once

namespace PhysicsUnits {
// Масса: 1 метрическая тонна = 1000 кг.
constexpr double kKgPerTon = 1000.0;
// Масса: 1 гигатонна = 10^9 тонн.
constexpr double kTonsPerGigaton = 1.0e9;
// Масса: 1 терратонна = 10^12 тонн.
constexpr double kTonsPerTerraton = 1.0e12;

// Конвертации массы в килограммы.
constexpr double kKgPerGigaton = kKgPerTon * kTonsPerGigaton;
constexpr double kKgPerTerraton = kKgPerTon * kTonsPerTerraton;

// Плотность: 1 г/см^3 = 1000 кг/м^3.
constexpr double kKgPerM3PerGramPerCm3 = 1000.0;
} // namespace PhysicsUnits
