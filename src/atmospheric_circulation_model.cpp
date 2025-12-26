#include "atmospheric_circulation_model.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kDegreesPerHemisphere = 90.0;
constexpr int kMinCellsPerHemisphere = 1;
constexpr int kMaxCellsPerHemisphere = 6;

int calculateCellsPerHemisphere(double dayLengthDays) {
    const double clampedDayLength = qBound(0.1, dayLengthDays, 50.0);
    // Параметризация разделения циркуляции на ячейки при изменении вращения:
    // чем быстрее ротация, тем больше ячеек; при медленной — меньше.
    const double rotationScale = std::sqrt(1.0 / clampedDayLength);
    const double rawCells = 3.0 * rotationScale;
    const int cells = static_cast<int>(std::round(rawCells));
    return qBound(kMinCellsPerHemisphere, cells, kMaxCellsPerHemisphere);
}
}  // namespace

AtmosphericCirculationModel::AtmosphericCirculationModel(double dayLengthDays,
                                                         RotationMode rotationMode,
                                                         double atmospherePressureAtm)
    : dayLengthDays_(dayLengthDays),
      rotationMode_(rotationMode),
      atmospherePressureAtm_(atmospherePressureAtm),
      cellsPerHemisphere_(calculateCellsPerHemisphere(dayLengthDays_)),
      cellSizeDegrees_(kDegreesPerHemisphere /
                       static_cast<double>(cellsPerHemisphere_)) {}

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
    for (const auto &point : points) {
        mixingCoefficients.push_back(meridionalMixingCoefficient(point.latitudeDegrees));
    }
    const int transportSteps = qMax(1, meridionalTransportSteps_);
    const double cellSizeDegrees = cellSizeDegrees_;
    const double dayLengthSeconds = qMax(0.1, dayLengthDays_) * 86400.0;
    const double rotationOmega = 2.0 * M_PI / dayLengthSeconds;
    const double windScalingAlpha = 0.12;
    const double minCoriolis = 1.0e-6;
    QVector<double> totalMixingCoefficients;
    totalMixingCoefficients.reserve(mixingCoefficients.size());
    double maxMixingCoefficient = 0.0;
    for (int i = 0; i < mixingCoefficients.size(); ++i) {
        double localGradient = 0.0;
        if (i > 0 && i + 1 < meanDaily.size()) {
            localGradient =
                (meanDaily.at(i + 1) - meanDaily.at(i - 1)) / (2.0 * cellSizeDegrees);
        } else if (i + 1 < meanDaily.size()) {
            localGradient = (meanDaily.at(i + 1) - meanDaily.at(i)) / cellSizeDegrees;
        } else if (i > 0) {
            localGradient = (meanDaily.at(i) - meanDaily.at(i - 1)) / cellSizeDegrees;
        }
        const double latitudeRadians = qDegreesToRadians(points.at(i).latitudeDegrees);
        const double coriolis = 2.0 * rotationOmega * std::sin(latitudeRadians);
        // Упрощенная оценка скорости ветра:
        // U ~ alpha * |dT/dy| / f, где f = 2 * Omega * sin(phi).
        // Это не динамическая модель, а эвристика связи температурного градиента и переноса.
        const double windSpeed =
            windScalingAlpha * std::abs(localGradient) / qMax(minCoriolis, std::abs(coriolis));
        // Переводим скорость ветра в эффективный коэффициент диффузии K ~ U * L,
        // где L — характерный масштаб переноса (ширина широтной ячейки).
        const double windMixingCoefficient = windSpeed * cellSizeDegrees;
        const double coefficient = mixingCoefficients.at(i) + windMixingCoefficient;
        totalMixingCoefficients.push_back(coefficient);
        maxMixingCoefficient = qMax(maxMixingCoefficient, coefficient);
    }
    // Условие устойчивости явной схемы диффузии: dt <= dx^2 / (2K).
    // Физически это означает, что перенос за один шаг не должен "перепрыгивать"
    // через соседнюю ячейку широт, иначе профиль температур начнет колебаться.
    const double transportTimeStep =
        (maxMixingCoefficient > 0.0)
            ? (0.9 * cellSizeDegrees * cellSizeDegrees / (2.0 * maxMixingCoefficient))
            : 0.0;
    QVector<double> diffusionFactors;
    diffusionFactors.reserve(totalMixingCoefficients.size());
    for (const double coefficient : totalMixingCoefficients) {
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
    const int cellIndex = static_cast<int>(std::floor(normalized / cellSizeDegrees_));
    return qBound(0, cellIndex, cellsPerHemisphere_ * 2 - 1);
}

double AtmosphericCirculationModel::cellCouplingFactor(int leftCell, int rightCell) const {
    if (leftCell == rightCell) {
        return 1.0;
    }
    // Между соседними ячейками перенос ослаблен: границы ячеек тормозят обмен.
    const double adjustment = 0.04 * (3.0 - static_cast<double>(cellsPerHemisphere_));
    return qBound(0.4, 0.6 + adjustment, 0.8);
}
