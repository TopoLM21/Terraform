#include "surface/VegetationModel.h"

#include <QtCore/QtGlobal>
#include <QtCore/QtMath>

namespace {
constexpr double kSecondsPerDay = 86400.0;
}

VegetationModel::VegetationModel() = default;

VegetationModel::VegetationModel(const Settings &settings)
    : settings_(settings) {}

void VegetationModel::setSettings(const Settings &settings) {
    settings_ = settings;
}

const VegetationModel::Settings &VegetationModel::settings() const {
    return settings_;
}

void VegetationModel::update(QVector<SurfacePoint> &points,
                             double timeStepSeconds,
                             double co2Share) const {
    if (points.isEmpty() || timeStepSeconds <= 0.0) {
        return;
    }

    const double growthRatePerSecond = settings_.growthRatePerDay / kSecondsPerDay;
    const double mortalityRatePerSecond = settings_.mortalityRatePerDay / kSecondsPerDay;
    const double diffusionRatePerSecond = settings_.diffusionRatePerDay / kSecondsPerDay;

    QVector<double> biomassAfterGrowth(points.size(), 0.0);

    for (int i = 0; i < points.size(); ++i) {
        SurfacePoint &point = points[i];
        const bool isOcean = (point.materialId == QLatin1String("ocean"));
        if (isOcean) {
            point.vegetationFraction = 0.0;
            point.vegetationBiomass = 0.0;
            point.vegetationBlendMask = 0.0;
            biomassAfterGrowth[i] = 0.0;
            continue;
        }

        // Парциальное давление CO₂: P_CO2 = P_total * x_CO2.
        // P_total берём из SurfacePoint, x_CO2 задаётся составом атмосферы.
        const double co2PartialAtm = qMax(0.0, point.pressureAtm) * qMax(0.0, co2Share);
        // Ограничивающие коэффициенты (0..1) для CO₂ и света.
        // Эти множители дополнительно "душат" рост при дефиците углекислого газа
        // или освещённости и участвуют в расчёте смертности через климатический стресс.
        const double co2Limit = co2Factor(co2PartialAtm);
        const double lightLimit = lightFactor(point.solarFluxWPerM2);
        const double climateFactor =
            temperatureFactor(point.temperatureK) *
            moistureFactor(point.surfaceMoisture.moistureFraction()) *
            pressureFactor(point.pressureAtm) *
            co2Limit *
            lightLimit;

        // Рост модели: логистическая формула с ограничением биомассы.
        // dB/dt = r * F_climate * B * (1 - B/B_max) - m * B,
        // где r — темп роста, m — смертность, F_climate — климатический фактор (0..1).
        const double currentBiomass = qMax(point.vegetationBiomass,
                                           point.vegetationFraction * settings_.maxBiomassKgPerM2);
        const double seedBiomass = qMax(currentBiomass, settings_.biomassSeedKgPerM2);
        const double growth = growthRatePerSecond * climateFactor *
            (1.0 - currentBiomass / settings_.maxBiomassKgPerM2) * seedBiomass;
        const double stress = 1.0 - climateFactor;
        const double mortality = mortalityRatePerSecond *
            (1.0 + settings_.stressMortalityBoost * stress) * currentBiomass;

        double updatedBiomass = currentBiomass + (growth - mortality) * timeStepSeconds;
        if (updatedBiomass < 0.0) {
            updatedBiomass = 0.0;
        }

        biomassAfterGrowth[i] = updatedBiomass;
        // Маска для визуального смешивания: напрямую связана с климатическим потенциалом.
        point.vegetationBlendMask = qBound(0.0, climateFactor, 1.0);
    }

    for (int i = 0; i < points.size(); ++i) {
        SurfacePoint &point = points[i];
        if (point.materialId == QLatin1String("ocean")) {
            continue;
        }

        double updatedBiomass = biomassAfterGrowth.at(i);
        const QVector<int> &neighbors = point.neighborIndices;
        if (!neighbors.isEmpty() && diffusionRatePerSecond > 0.0) {
            double neighborSum = 0.0;
            int neighborCount = 0;
            for (int neighborIndex : neighbors) {
                if (neighborIndex < 0 || neighborIndex >= points.size()) {
                    continue;
                }
                const SurfacePoint &neighbor = points.at(neighborIndex);
                if (neighbor.materialId == QLatin1String("ocean")) {
                    continue;
                }
                neighborSum += biomassAfterGrowth.at(neighborIndex);
                ++neighborCount;
            }
            if (neighborCount > 0) {
                const double neighborAverage = neighborSum / neighborCount;
                // Диффузия: B += D * (B_neighbors - B_self), где D — коэффициент обмена.
                updatedBiomass += diffusionRatePerSecond * timeStepSeconds *
                    (neighborAverage - updatedBiomass);
            }
        }

        updatedBiomass = qBound(0.0, updatedBiomass, settings_.maxBiomassKgPerM2);
        point.vegetationBiomass = updatedBiomass;
        point.vegetationFraction =
            settings_.maxBiomassKgPerM2 > 0.0
            ? qBound(0.0, updatedBiomass / settings_.maxBiomassKgPerM2, 1.0)
            : 0.0;
    }
}

double VegetationModel::temperatureFactor(double temperatureK) const {
    if (temperatureK <= 0.0) {
        return 0.0;
    }
    const double minTemp = settings_.minTempK;
    const double maxTemp = settings_.maxTempK;
    if (temperatureK <= minTemp || temperatureK >= maxTemp) {
        return 0.0;
    }

    // Треугольная функция с максимумом в optimalTemp:
    // 0 на min/max, 1 в оптимуме, линейный спад по обе стороны.
    const double span = maxTemp - minTemp;
    if (span <= 0.0) {
        return 0.0;
    }
    const double optimal = qBound(minTemp, settings_.optimalTempK, maxTemp);
    const double leftSpan = qMax(1.0, optimal - minTemp);
    const double rightSpan = qMax(1.0, maxTemp - optimal);
    if (temperatureK <= optimal) {
        return qBound(0.0, (temperatureK - minTemp) / leftSpan, 1.0);
    }
    return qBound(0.0, (maxTemp - temperatureK) / rightSpan, 1.0);
}

double VegetationModel::moistureFactor(double moistureFraction) const {
    if (moistureFraction <= settings_.minMoistureFraction) {
        return 0.0;
    }
    const double span = 1.0 - settings_.minMoistureFraction;
    return span > 0.0
        ? qBound(0.0, (moistureFraction - settings_.minMoistureFraction) / span, 1.0)
        : 0.0;
}

double VegetationModel::pressureFactor(double pressureAtm) const {
    if (pressureAtm <= settings_.minPressureAtm) {
        return 0.0;
    }
    // Линейно доводим до 1.0 при 1 атм, чтобы отражать нехватку воздуха на высоте.
    const double span = qMax(0.1, 1.0 - settings_.minPressureAtm);
    return qBound(0.0, (pressureAtm - settings_.minPressureAtm) / span, 1.0);
}

double VegetationModel::co2Factor(double co2PartialAtm) const {
    if (co2PartialAtm <= settings_.minCo2PartialAtm) {
        return 0.0;
    }
    // Насыщение по CO₂: при 3×минимуме считаем фотосинтез почти максимальным.
    const double target = settings_.minCo2PartialAtm * 3.0;
    return target > 0.0 ? qBound(0.0, co2PartialAtm / target, 1.0) : 0.0;
}

double VegetationModel::lightFactor(double solarFluxWPerM2) const {
    if (solarFluxWPerM2 <= 0.0) {
        return 0.0;
    }
    // Гиперболическое насыщение: F = 1 - exp(-I / I_sat),
    // где I_sat задаёт поток света, при котором рост близок к максимуму.
    const double saturation = qMax(1.0, settings_.lightSaturationWPerM2);
    return 1.0 - qExp(-solarFluxWPerM2 / saturation);
}
