#include "atmospheric_grid_3d.h"

#include "atmospheric_thermodynamics.h"

#include <QHash>
#include <QtMath>

namespace {
constexpr double kEarthMassKg = 5.9722e24;
constexpr double kGravitationalConstant = 6.67430e-11;
constexpr double kStandardPressurePa = 101325.0;

constexpr int kDefaultLayerCount = 12;

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
                                  int columnCount,
                                  int layerCount) {
    const int resolvedLayerCount = (layerCount > 0) ? layerCount : kDefaultLayerCount;
    resizeColumns(columnCount, resolvedLayerCount);
    if (columns_.isEmpty()) {
        return;
    }

    const double surfacePressureAtm = composition.totalPressureAtm(planetMassEarths, radiusKm);
    const double gravity = surfaceGravityMps2(planetMassEarths, radiusKm);
    const double rSpecific = AtmosphericThermodynamics::specificGasConstant(composition);
    const double specificHeat = AtmosphericThermodynamics::specificHeatCp(composition);

    // Простейшее допущение: изотермическая атмосфера у поверхности.
    const double baseTemperatureKelvin = (surfacePressureAtm > 0.0) ? 288.0 : 0.0;

    // Масштабная высота H = R * T / g.
    const double scaleHeightMeters =
        (rSpecific > 0.0 && gravity > 0.0 && baseTemperatureKelvin > 0.0)
            ? (rSpecific * baseTemperatureKelvin) / gravity
            : 0.0;

    // Высоту атмосферы берем как ~6 масштабных высот, где давление падает до ~0.2%.
    const double topHeightMeters = (scaleHeightMeters > 0.0)
        ? qMax(1000.0, scaleHeightMeters * 6.0)
        : 0.0;
    const double layerThicknessMeters = (resolvedLayerCount > 0)
        ? topHeightMeters / static_cast<double>(resolvedLayerCount)
        : 0.0;

    // Сухой адиабатический градиент: dT/dz = g / Cp.
    const double lapseRateKPerM = AtmosphericThermodynamics::dryAdiabaticLapseRate(
        composition, gravity);

    // Оптическая толщина: увеличивается с долей парниковых газов и распределяется по слоям.
    const double greenhouseShare = greenhouseMassFraction(composition);
    const double baseTauSw = 0.02 + 0.3 * greenhouseShare;
    const double baseTauLw = 0.05 + 1.2 * greenhouseShare;
    const double tauScale = (topHeightMeters > 0.0) ? (layerThicknessMeters / topHeightMeters) : 0.0;

    for (auto &column : columns_) {
        column.resize(resolvedLayerCount);
        auto &layers = column.layers();
        for (int i = 0; i < layers.size(); ++i) {
            const double heightMidMeters = (static_cast<double>(i) + 0.5) * layerThicknessMeters;
            const double temperatureKelvin = (baseTemperatureKelvin > 0.0)
                ? qMax(1.0, baseTemperatureKelvin - lapseRateKPerM * heightMidMeters)
                : 0.0;
            const double pressureAtm =
                (surfacePressureAtm > 0.0 && scaleHeightMeters > 0.0)
                    ? surfacePressureAtm * qExp(-heightMidMeters / scaleHeightMeters)
                    : 0.0;
            const double pressurePa = pressureAtm * kStandardPressurePa;
            // Плотность по уравнению идеального газа: rho = P / (R * T).
            const double densityKgPerM3 =
                (rSpecific > 0.0 && temperatureKelvin > 0.0)
                    ? (pressurePa / (rSpecific * temperatureKelvin))
                    : 0.0;
            // Теплоёмкость слоя на единицу площади: C = rho * Cp * dz.
            const double heatCapacityJPerM2K =
                (specificHeat > 0.0)
                    ? densityKgPerM3 * specificHeat * layerThicknessMeters
                    : 0.0;

            const bool convectionEnabled = lapseRateKPerM > 0.0 && temperatureKelvin > 0.0;
            const double convectionMixingCoefficient = convectionEnabled
                ? qBound(0.0,
                         (lapseRateKPerM * layerThicknessMeters) /
                             qMax(1.0, baseTemperatureKelvin),
                         1.0)
                : 0.0;

            layers[i].setTemperatureKelvin(temperatureKelvin);
            layers[i].setPressureAtm(pressureAtm);
            layers[i].setDensityKgPerM3(densityKgPerM3);
            layers[i].setHeightMeters(heightMidMeters);
            layers[i].setThicknessMeters(layerThicknessMeters);
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
