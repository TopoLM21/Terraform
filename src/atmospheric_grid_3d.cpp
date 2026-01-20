#include "atmospheric_grid_3d.h"

#include "atmospheric_profile_initializer.h"
#include "atmospheric_thermodynamics.h"

#include <QHash>
#include <QtMath>

namespace {
constexpr double kEarthMassKg = 5.9722e24;
constexpr double kGravitationalConstant = 6.67430e-11;
constexpr int kDefaultLayerCount = 12;
constexpr double kDefaultBottomLayerThicknessMeters = 100.0;

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
                                  double minBottomLayerThicknessMeters) {
    const int resolvedLayerCount = (layerCount > 0) ? layerCount : kDefaultLayerCount;
    resizeColumns(columnCount, resolvedLayerCount);
    if (columns_.isEmpty()) {
        return;
    }

    const double surfacePressureAtm =
        composition.totalPressureAtm(planetMassEarths, radiusKm);
    const double gravity = surfaceGravityMps2(planetMassEarths, radiusKm);
    const double specificHeat = AtmosphericThermodynamics::specificHeatCp(composition);

    // Простейшее допущение: заданная температура на поверхности.
    const double resolvedBaseTemperatureKelvin =
        (surfacePressureAtm > 0.0) ? qMax(0.0, baseTemperatureKelvin) : 0.0;

    AtmosphericProfileInitializer initializer(composition);
    AtmosphericProfileInitializer::Settings profileSettings;
    profileSettings.surfaceTemperatureKelvin = resolvedBaseTemperatureKelvin;
    profileSettings.surfacePressureAtm = surfacePressureAtm;
    profileSettings.gravityMps2 = gravity;
    profileSettings.layerCount = resolvedLayerCount;
    profileSettings.useDryAdiabatic = true;
    profileSettings.minBottomLayerThicknessMeters =
        (minBottomLayerThicknessMeters > 0.0)
            ? minBottomLayerThicknessMeters
            : kDefaultBottomLayerThicknessMeters;

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
