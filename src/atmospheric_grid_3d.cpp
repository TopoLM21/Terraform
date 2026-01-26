#include "atmospheric_grid_3d.h"

#include "atmospheric_layer_resampler.h"
#include "atmospheric_profile_initializer.h"
#include "atmospheric_thermodynamics.h"

#include <QHash>
#include <QtMath>

namespace {
constexpr double kEarthMassKg = 5.9722e24;
constexpr double kGravitationalConstant = 6.67430e-11;
constexpr int kDefaultLayerCount = 12;
constexpr double kMinBottomLayerThicknessMeters = 250.0;
constexpr double kMaxBottomLayerThicknessMeters = 500.0;
constexpr double kDefaultBottomLayerThicknessMeters = 250.0;
double greenhouseMassFraction(const AtmosphereComposition &composition) {
    const auto gases = availableGases();
    QHash<QString, bool> isGreenhouse;
    isGreenhouse.reserve(gases.size());
    for (const auto &gas : gases) {
        isGreenhouse.insert(gas.id, gas.isGreenhouse);
    }

    const double totalMass = composition.totalMassGigatons();
    if (totalMass <= 0.0) {
        return 0.0;
    }

    double greenhouseMass = 0.0;
    const auto fractions = composition.fractions();
    for (const auto &fraction : fractions) {
        if (isGreenhouse.value(fraction.id, false)) {
            greenhouseMass += fraction.massGigatons;
        }
    }

    return greenhouseMass / totalMass;
}

double surfaceGravityMps2(double planetMassEarths, double radiusKm) {
    if (planetMassEarths <= 0.0 || radiusKm <= 0.0) {
        return 0.0;
    }
    const double radiusMeters = radiusKm * 1000.0;
    const double planetMassKg = planetMassEarths * kEarthMassKg;
    // g = G * M / R^2
    return kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);
}

int resolveLayerCountForTopPressure(double surfacePressureAtm,
                                    double minTopPressureAtm,
                                    double surfaceTemperatureKelvin,
                                    double gravityMps2,
                                    double specificGasConstant) {
    if (surfacePressureAtm <= 0.0 || minTopPressureAtm <= 0.0 ||
        surfacePressureAtm <= minTopPressureAtm || surfaceTemperatureKelvin <= 0.0 ||
        gravityMps2 <= 0.0 || specificGasConstant <= 0.0) {
        return 0;
    }

    // Масштабная высота: H = R * T / g, затем используем её для оценки высоты до P_min.
    const double scaleHeight =
        (specificGasConstant * surfaceTemperatureKelvin) / gravityMps2;
    if (scaleHeight <= 0.0) {
        return 0;
    }

    const double topHeightMeters = scaleHeight * qLn(surfacePressureAtm / minTopPressureAtm);
    if (topHeightMeters <= 0.0) {
        return 0;
    }

    const double defaultTopHeight = scaleHeight * 6.0;
    const double targetLayerThickness =
        (defaultTopHeight > 0.0)
            ? (defaultTopHeight / static_cast<double>(kDefaultLayerCount))
            : 0.0;
    if (targetLayerThickness <= 0.0) {
        return 0;
    }

    return qMax(1, static_cast<int>(qCeil(topHeightMeters / targetLayerThickness)));
}
} // namespace

void AtmosphericGrid3D::resizeColumns(int columnCount, int layerCount) {
    columns_.resize(columnCount);
    if (layerCount > 0) {
        layerCount_ = layerCount;
        for (auto &column : columns_) {
            column.resize(layerCount);
        }
    }
}

int AtmosphericGrid3D::columnCount() const {
    return columns_.size();
}

int AtmosphericGrid3D::layerCount() const {
    return layerCount_;
}

double AtmosphericGrid3D::minTopPressureAtm() const {
    return minTopPressureAtm_;
}

QVector<AtmosphericColumn> &AtmosphericGrid3D::columns() {
    return columns_;
}

const QVector<AtmosphericColumn> &AtmosphericGrid3D::columns() const {
    return columns_;
}

void AtmosphericGrid3D::initialize(const AtmosphereComposition &composition,
                                  double planetMassEarths,
                                  double radiusKm,
                                  double baseTemperatureKelvin,
                                  int columnCount,
                                  int layerCount,
                                  double minBottomLayerThicknessMeters,
                                  double minTopPressureAtm) {
    minTopPressureAtm_ = minTopPressureAtm;
    // Сохраняем фактический минимум толщины нижнего слоя, чтобы повторять
    // «уплотнённую» геометрию при пересборке сетки.
    minBottomLayerThicknessMeters_ = qBound(
        kMinBottomLayerThicknessMeters,
        (minBottomLayerThicknessMeters > 0.0)
            ? minBottomLayerThicknessMeters
            : kDefaultBottomLayerThicknessMeters,
        kMaxBottomLayerThicknessMeters);
    const double surfacePressureAtm =
        composition.totalPressureAtm(planetMassEarths, radiusKm);
    const double gravity = surfaceGravityMps2(planetMassEarths, radiusKm);
    const double specificHeat = AtmosphericThermodynamics::specificHeatCp(composition);
    const double rSpecific = AtmosphericThermodynamics::specificGasConstant(composition);

    // Простейшее допущение: заданная температура на поверхности.
    const double resolvedBaseTemperatureKelvin =
        (surfacePressureAtm > 0.0) ? qMax(0.0, baseTemperatureKelvin) : 0.0;

    int resolvedLayerCount = (layerCount > 0) ? layerCount : kDefaultLayerCount;
    if (layerCount <= 0 && minTopPressureAtm > 0.0) {
        // Для сверхплотных атмосфер выбираем верх по давлению, чтобы слои покрывали
        // физически значимую высоту (например, 0.1–0.5 атм для Венеры).
        const int pressureDrivenLayerCount = resolveLayerCountForTopPressure(
            surfacePressureAtm,
            minTopPressureAtm,
            resolvedBaseTemperatureKelvin,
            gravity,
            rSpecific);
        if (pressureDrivenLayerCount > 0) {
            resolvedLayerCount = pressureDrivenLayerCount;
        }
    }

    resizeColumns(columnCount, resolvedLayerCount);
    if (columns_.isEmpty()) {
        return;
    }

    AtmosphericProfileInitializer initializer(composition);
    AtmosphericProfileInitializer::Settings profileSettings;
    profileSettings.surfaceTemperatureKelvin = resolvedBaseTemperatureKelvin;
    profileSettings.surfacePressureAtm = surfacePressureAtm;
    profileSettings.gravityMps2 = gravity;
    profileSettings.layerCount = resolvedLayerCount;
    profileSettings.useDryAdiabatic = true;
    profileSettings.minTopPressureAtm = minTopPressureAtm;
    // Жёстко ограничиваем "обитаемый слой" 250–500 м как компромисс между
    // стабильной визуализацией и физической точностью.
    profileSettings.minBottomLayerThicknessMeters = minBottomLayerThicknessMeters_;

    const QVector<AtmosphericProfileLayer> profileLayers =
        initializer.buildProfile(profileSettings);

    // Оптическая толщина: увеличивается с долей парниковых газов и распределяется по слоям.
    const double greenhouseShare = greenhouseMassFraction(composition);
    const double baseTauSw = 0.02 + 0.3 * greenhouseShare;
    const double baseTauLw = 0.05 + 1.2 * greenhouseShare;
    const double dryLapseRateKPerM = AtmosphericThermodynamics::dryAdiabaticLapseRate(
        composition, gravity);
    double topHeightMeters = 0.0;
    for (const auto &profileLayer : profileLayers) {
        topHeightMeters += profileLayer.thicknessMeters;
    }

    for (auto &column : columns_) {
        column.resize(resolvedLayerCount);
        auto &layers = column.layers();
        for (int i = 0; i < layers.size(); ++i) {
            const AtmosphericProfileLayer profileLayer =
                (i < profileLayers.size()) ? profileLayers[i] : AtmosphericProfileLayer{};
            // Теплоёмкость слоя на единицу площади: C = rho * Cp * dz.
            const double heatCapacityJPerM2K =
                (specificHeat > 0.0)
                    ? profileLayer.densityKgPerM3 * specificHeat * profileLayer.thicknessMeters
                    : 0.0;

            const bool convectionEnabled =
                profileLayer.temperatureKelvin > 0.0 && profileLayer.thicknessMeters > 0.0;
            const double convectionMixingCoefficient = convectionEnabled
                ? qBound(0.0,
                         (dryLapseRateKPerM * profileLayer.thicknessMeters) /
                             qMax(1.0, profileLayer.temperatureKelvin),
                         1.0)
                : 0.0;
            const double tauScale = (topHeightMeters > 0.0)
                ? (profileLayer.thicknessMeters / topHeightMeters)
                : 0.0;

            layers[i].setTemperatureKelvin(profileLayer.temperatureKelvin);
            layers[i].setPressureAtm(profileLayer.pressureAtm);
            layers[i].setDensityKgPerM3(profileLayer.densityKgPerM3);
            layers[i].setHeightMeters(profileLayer.heightMeters);
            layers[i].setThicknessMeters(profileLayer.thicknessMeters);
            layers[i].setHeatCapacityJPerM2K(heatCapacityJPerM2K);
            layers[i].setWindUMps(0.0);
            layers[i].setWindVMps(0.0);
            layers[i].setOpticalDepthShortwave(baseTauSw * tauScale);
            layers[i].setOpticalDepthLongwave(baseTauLw * tauScale);
            layers[i].setConvectionEnabled(convectionEnabled);
            layers[i].setConvectionMixingCoefficient(convectionMixingCoefficient);
        }
    }
}

bool AtmosphericGrid3D::updateLayerCountForTopPressure(double surfacePressureAtm,
                                                       double minTopPressureAtm,
                                                       double surfaceTemperatureKelvin,
                                                       double gravityMps2,
                                                       double specificGasConstant,
                                                       double changeThresholdFraction) {
    if (surfacePressureAtm <= 0.0 || minTopPressureAtm <= 0.0 ||
        surfacePressureAtm <= minTopPressureAtm || surfaceTemperatureKelvin <= 0.0 ||
        gravityMps2 <= 0.0 || specificGasConstant <= 0.0) {
        return false;
    }

    minTopPressureAtm_ = minTopPressureAtm;
    // Масштабная высота: H = R * T / g (см. resolveLayerCountForTopPressure).
    const double scaleHeightMeters =
        (specificGasConstant * surfaceTemperatureKelvin) / gravityMps2;
    if (scaleHeightMeters <= 0.0) {
        return false;
    }

    const double topHeightMeters =
        scaleHeightMeters * qLn(surfacePressureAtm / minTopPressureAtm);
    if (topHeightMeters <= 0.0) {
        return false;
    }

    const int newLayerCount = resolveLayerCountForTopPressure(surfacePressureAtm,
                                                              minTopPressureAtm,
                                                              surfaceTemperatureKelvin,
                                                              gravityMps2,
                                                              specificGasConstant);
    if (newLayerCount <= 0) {
        return false;
    }

    auto relativeChange = [](double value, double previous) -> double {
        return (previous > 0.0) ? qAbs(value - previous) / previous : 0.0;
    };

    const double scaleHeightDelta =
        relativeChange(scaleHeightMeters, lastScaleHeightMeters_);
    const double surfacePressureDelta =
        relativeChange(surfacePressureAtm, lastSurfacePressureAtm_);
    const double topHeightDelta =
        relativeChange(topHeightMeters, lastTopHeightMeters_);

    const bool hasBaseline =
        lastScaleHeightMeters_ > 0.0 && lastSurfacePressureAtm_ > 0.0 &&
        lastTopHeightMeters_ > 0.0;
    const bool thresholdExceeded = (changeThresholdFraction <= 0.0)
        ? true
        : (!hasBaseline ||
           scaleHeightDelta > changeThresholdFraction ||
           surfacePressureDelta > changeThresholdFraction ||
           topHeightDelta > changeThresholdFraction);

    lastScaleHeightMeters_ = scaleHeightMeters;
    lastSurfacePressureAtm_ = surfacePressureAtm;
    lastTopHeightMeters_ = topHeightMeters;

    if (!thresholdExceeded || columns_.isEmpty()) {
        return false;
    }

    rebuildLayersWithInterpolation(newLayerCount, topHeightMeters);
    return true;
}

void AtmosphericGrid3D::rebuildLayersWithInterpolation(int newLayerCount,
                                                       double topHeightMeters) {
    if (newLayerCount <= 0 || topHeightMeters <= 0.0) {
        return;
    }

    for (auto &column : columns_) {
        const QVector<AtmosphericLayerState> oldLayers = column.layers();
        const int remainingLayers = newLayerCount - 1;
        const double requestedBottomThickness = qBound(
            kMinBottomLayerThicknessMeters,
            (minBottomLayerThicknessMeters_ > 0.0)
                ? minBottomLayerThicknessMeters_
                : kDefaultBottomLayerThicknessMeters,
            kMaxBottomLayerThicknessMeters);
        // Если остаётся один слой, покрываем всю высоту, сохраняя нижнюю «плотность»
        // только для многослойной сетки.
        const double bottomThickness = (remainingLayers == 0)
            ? topHeightMeters
            : qMin(requestedBottomThickness, topHeightMeters);
        const double remainingHeight =
            (remainingLayers > 0) ? qMax(0.0, topHeightMeters - bottomThickness) : 0.0;
        const double uniformThickness =
            (remainingLayers > 0) ? (remainingHeight / remainingLayers) : 0.0;

        QVector<AtmosphericLayerResampler::TargetLayer> targets;
        targets.reserve(newLayerCount);
        double bottomEdgeMeters = 0.0;
        targets.push_back(AtmosphericLayerResampler::TargetLayer{
            bottomEdgeMeters + bottomThickness * 0.5, bottomThickness});
        bottomEdgeMeters += bottomThickness;
        for (int i = 0; i < remainingLayers; ++i) {
            const double heightMeters = bottomEdgeMeters + uniformThickness * 0.5;
            targets.push_back(AtmosphericLayerResampler::TargetLayer{
                heightMeters, uniformThickness});
            bottomEdgeMeters += uniformThickness;
        }

        const QVector<AtmosphericLayerState> newLayers =
            AtmosphericLayerResampler::resample(oldLayers, targets);
        if (newLayers.isEmpty()) {
            return;
        }
        column.resize(newLayerCount);
        auto &layers = column.layers();
        for (int i = 0; i < layers.size(); ++i) {
            layers[i] = newLayers.at(i);
        }
    }

    layerCount_ = newLayerCount;
}
