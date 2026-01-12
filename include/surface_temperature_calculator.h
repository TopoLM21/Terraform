#pragma once

#include <QtCore/QVector>
struct TemperatureRangePoint {
    double latitudeDegrees;
    bool hasInsolation = false;
    double minimumKelvin;
    double maximumKelvin;
    double meanDailyKelvin;
    double meanDayKelvin;
    double meanNightKelvin;
    double minimumCelsius;
    double maximumCelsius;
    double meanDailyCelsius;
    double meanDayCelsius;
    double meanNightCelsius;
};

struct TemperatureSummaryPoint {
    double latitudeDegrees;
    double minimumKelvin;
    double maximumKelvin;
    double meanAnnualKelvin;
    double meanAnnualDayKelvin;
    double meanAnnualNightKelvin;
    double minimumCelsius;
    double maximumCelsius;
    double meanAnnualCelsius;
    double meanAnnualDayCelsius;
    double meanAnnualNightCelsius;
};
