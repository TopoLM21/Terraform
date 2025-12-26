#include "atmospheric_circulation_model.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr int kCellsPerHemisphere = 3;
constexpr double kDegreesPerHemisphere = 90.0;
constexpr double kCellSizeDegrees =
    (kDegreesPerHemisphere / static_cast<double>(kCellsPerHemisphere));
}  // namespace

AtmosphericCirculationModel::AtmosphericCirculationModel(double dayLengthDays,
                                                         RotationMode rotationMode,
                                                         double atmospherePressureAtm)
    : dayLengthDays_(dayLengthDays),
      rotationMode_(rotationMode),
      atmospherePressureAtm_(atmospherePressureAtm) {}

void AtmosphericCirculationModel::setAtmosphereMassKg(double atmosphereMassKg) {
    atmosphereMassKg_ = atmosphereMassKg;
}

void AtmosphericCirculationModel::setMeridionalTransportSteps(int steps) {
    meridionalTransportSteps_ = qMax(1, steps);
}

QVector<TemperatureRangePoint> AtmosphericCirculationModel::applyHeatTransport(
    const QVector<TemperatureRangePoint> &points) const {
    if (points.size() < 2) {
        return points;
    }

    QVector<double> meanDaily;
    meanDaily.reserve(points.size());
    for (const auto &point : points) {
        meanDaily.push_back(point.meanDailyKelvin);
    }

    QVector<double> mixingCoefficients;
    mixingCoefficients.reserve(points.size());
    double maxMixingCoefficient = 0.0;
    for (const auto &point : points) {
        const double coefficient = meridionalMixingCoefficient(point.latitudeDegrees);
        mixingCoefficients.push_back(coefficient);
        maxMixingCoefficient = qMax(maxMixingCoefficient, coefficient);
    }
    const int transportSteps = qMax(1, meridionalTransportSteps_);
    const double cellSizeDegrees = kCellSizeDegrees;
    // Условие устойчивости явной схемы диффузии: dt <= dx^2 / (2K).
    // Физически это означает, что перенос за один шаг не должен "перепрыгивать"
    // через соседнюю ячейку широт, иначе профиль температур начнет колебаться.
    const double transportTimeStep =
        (maxMixingCoefficient > 0.0)
            ? (0.9 * cellSizeDegrees * cellSizeDegrees / (2.0 * maxMixingCoefficient))
            : 0.0;
    QVector<double> diffusionFactors;
    diffusionFactors.reserve(mixingCoefficients.size());
    for (const double coefficient : mixingCoefficients) {
        const double diffusionFactor =
            (maxMixingCoefficient > 0.0)
                ? (coefficient * transportTimeStep / (cellSizeDegrees * cellSizeDegrees))
                : 0.0;
        diffusionFactors.push_back(diffusionFactor);
    }
    QVector<double> mixed = meanDaily;
    QVector<double> next = meanDaily;
    for (int step = 0; step < transportSteps; ++step) {
        for (int i = 0; i < mixed.size(); ++i) {
            double diffusionTerm = 0.0;
            const int currentCell = cellIndexForAxis(points.at(i).latitudeDegrees);
            if (i > 0) {
                const int leftCell = cellIndexForAxis(points.at(i - 1).latitudeDegrees);
                const double coupling = cellCouplingFactor(leftCell, currentCell);
                diffusionTerm += coupling * (mixed.at(i - 1) - mixed.at(i));
            }
            if (i + 1 < mixed.size()) {
                const int rightCell = cellIndexForAxis(points.at(i + 1).latitudeDegrees);
                const double coupling = cellCouplingFactor(currentCell, rightCell);
                diffusionTerm += coupling * (mixed.at(i + 1) - mixed.at(i));
            }
            next[i] = mixed.at(i) + diffusionFactors.at(i) * diffusionTerm;
        }
        mixed.swap(next);
    }

    if (rotationMode_ == RotationMode::TidalLocked) {
        double daySum = 0.0;
        double nightSum = 0.0;
        int dayCount = 0;
        int nightCount = 0;
        for (int i = 0; i < mixed.size(); ++i) {
            if (points.at(i).hasInsolation) {
                daySum += mixed.at(i);
                ++dayCount;
            } else {
                nightSum += mixed.at(i);
                ++nightCount;
            }
        }
        const double dayMean = (dayCount > 0) ? (daySum / dayCount) : 0.0;
        const double nightMean = (nightCount > 0) ? (nightSum / nightCount) : 0.0;
        const double exchangeCoefficient = dayNightExchangeCoefficient();
        for (int i = 0; i < mixed.size(); ++i) {
            if (points.at(i).hasInsolation) {
                mixed[i] += exchangeCoefficient * (nightMean - mixed.at(i));
            } else {
                mixed[i] += exchangeCoefficient * (dayMean - mixed.at(i));
            }
        }
    }

    QVector<TemperatureRangePoint> adjusted = points;
    for (int i = 0; i < adjusted.size(); ++i) {
        const double delta = mixed.at(i) - adjusted.at(i).meanDailyKelvin;
        auto &point = adjusted[i];
        point.meanDailyKelvin = qMax(1.0, point.meanDailyKelvin + delta);
        point.meanDayKelvin = qMax(1.0, point.meanDayKelvin + delta);
        point.meanNightKelvin = qMax(1.0, point.meanNightKelvin + delta);
        point.minimumKelvin = qMax(1.0, point.minimumKelvin + delta);
        point.maximumKelvin = qMax(1.0, point.maximumKelvin + delta);
        point.minimumCelsius = point.minimumKelvin - 273.15;
        point.maximumCelsius = point.maximumKelvin - 273.15;
        point.meanDailyCelsius = point.meanDailyKelvin - 273.15;
        point.meanDayCelsius = point.meanDayKelvin - 273.15;
        point.meanNightCelsius = point.meanNightKelvin - 273.15;
    }

    return adjusted;
}

double AtmosphericCirculationModel::baseMeridionalMixingCoefficient() const {
    if (atmospherePressureAtm_ <= 0.0) {
        return 0.0;
    }
    const double clampedDayLength = qBound(0.1, dayLengthDays_, 50.0);
    // Эффективный коэффициент горизонтального переноса уменьшается с быстрой ротацией
    // (разбивка на узкие ячейки и барьерный эффект струй) и растет с ростом давления.
    // Простая параметризация: k ~ sqrt(dayLength) * sqrt(pressure).
    const double rotationFactor = std::sqrt(clampedDayLength);
    const double pressureFactor = std::sqrt(qMax(0.01, atmospherePressureAtm_));
    const double rawCoefficient = 0.06 * rotationFactor * pressureFactor;
    return qBound(0.0, rawCoefficient, 0.35);
}

double AtmosphericCirculationModel::meridionalMixingCoefficient(double latitudeDegrees) const {
    const double baseCoefficient = baseMeridionalMixingCoefficient();
    if (baseCoefficient <= 0.0) {
        return 0.0;
    }
    const double latitudeRadians = qDegreesToRadians(latitudeDegrees);
    // Параметризация циркуляционных ячеек (Хэдли/Ферреля): перенос максимален
    // в средних широтах и минимален у экватора и полюсов.
    // Профиль: K(φ) = K0 * (0.25 + 0.75 * |sin(2φ)|), максимум при φ≈45°.
    const double profileFactor =
        0.25 + 0.75 * std::abs(std::sin(2.0 * latitudeRadians));
    return baseCoefficient * profileFactor;
}

double AtmosphericCirculationModel::dayNightExchangeCoefficient() const {
    if (atmospherePressureAtm_ <= 0.0) {
        return 0.0;
    }
    // На приливно-связанной планете перенос между дневной и ночной сторонами
    // усиливается с ростом массы атмосферы: больше масса → сильнее теплоемкость и ветры.
    // Масштабируем относительно земной массы атмосферы (~5.1e18 кг) и давления.
    const double pressureFactor = std::sqrt(qMax(0.01, atmospherePressureAtm_));
    const double massFactor =
        (atmosphereMassKg_ > 0.0)
            ? std::cbrt(qMax(0.05, atmosphereMassKg_ / 5.1e18))
            : 1.0;
    const double rawCoefficient = 0.08 * pressureFactor * massFactor;
    return qBound(0.0, rawCoefficient, 0.45);
}

int AtmosphericCirculationModel::cellIndexForAxis(double axisDegrees) const {
    const double normalized =
        (rotationMode_ == RotationMode::Normal) ? (axisDegrees + kDegreesPerHemisphere)
                                                : axisDegrees;
    const int cellIndex = static_cast<int>(std::floor(normalized / kCellSizeDegrees));
    return qBound(0, cellIndex, kCellsPerHemisphere * 2 - 1);
}

double AtmosphericCirculationModel::cellCouplingFactor(int leftCell, int rightCell) const {
    if (leftCell == rightCell) {
        return 1.0;
    }
    // Между соседними ячейками перенос ослаблен: границы ячеек тормозят обмен.
    return 0.6;
}
