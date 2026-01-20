#include "atmospheric_profile_initializer.h"

#include "atmospheric_thermodynamics.h"

#include <QtMath>

namespace {
constexpr double kStandardPressurePa = 101325.0;
}

AtmosphericProfileInitializer::AtmosphericProfileInitializer(
    const AtmosphereComposition &composition)
    : composition_(composition) {}

QVector<AtmosphericProfileLayer> AtmosphericProfileInitializer::buildProfile(
    const Settings &settings) const {
    QVector<AtmosphericProfileLayer> layers;
    if (settings.layerCount <= 0) {
        return layers;
    }

    layers.resize(settings.layerCount);

    const double rSpecific = AtmosphericThermodynamics::specificGasConstant(composition_);
    const double lapseRateKPerM = resolveLapseRateKPerM(settings);
    const double scaleHeight =
        scaleHeightMeters(settings.surfaceTemperatureKelvin, settings.gravityMps2);

    const double resolvedTopHeight = (settings.topHeightMeters > 0.0)
        ? settings.topHeightMeters
        : ((scaleHeight > 0.0) ? qMax(1000.0, scaleHeight * 6.0) : 0.0);

    const double uniformLayerThickness = (resolvedTopHeight > 0.0)
        ? resolvedTopHeight / static_cast<double>(settings.layerCount)
        : 0.0;
    double bottomLayerThickness = uniformLayerThickness;
    double remainingLayerThickness = uniformLayerThickness;
    bool useNonUniformLayers = settings.layerCount > 1 && resolvedTopHeight > 0.0 &&
        settings.minBottomLayerThicknessMeters > 0.0;
    if (useNonUniformLayers) {
        bottomLayerThickness =
            qMin(settings.minBottomLayerThicknessMeters, resolvedTopHeight);
        const double remainingHeight = resolvedTopHeight - bottomLayerThickness;
        if (remainingHeight > 0.0) {
            remainingLayerThickness =
                remainingHeight / static_cast<double>(settings.layerCount - 1);
        } else {
            useNonUniformLayers = false;
            bottomLayerThickness = uniformLayerThickness;
            remainingLayerThickness = uniformLayerThickness;
        }
    }

    // Разбиваем профиль так, чтобы первый слой имел заданную толщину, а остальные
    // делили оставшуюся высоту. Так сетка плотнее у поверхности — там важнее
    // детализация для визуализации и обитаемого объёма.
    double bottomEdgeMeters = 0.0;
    for (int i = 0; i < layers.size(); ++i) {
        const double layerThickness = useNonUniformLayers
            ? ((i == 0) ? bottomLayerThickness : remainingLayerThickness)
            : uniformLayerThickness;
        const double heightMidMeters = bottomEdgeMeters + layerThickness * 0.5;
        bottomEdgeMeters += layerThickness;
        // Температура по линейному градиенту: T(z) = T0 - Γ * z.
        const double temperatureKelvin = (settings.surfaceTemperatureKelvin > 0.0)
            ? qMax(1.0, settings.surfaceTemperatureKelvin - lapseRateKPerM * heightMidMeters)
            : 0.0;
        // Барометрическая формула: P(z) = P0 * exp(-z / H).
        const double pressureAtm = pressureAtmAtHeight(settings.surfacePressureAtm,
                                                       heightMidMeters,
                                                       scaleHeight);
        const double pressurePa = pressureAtm * kStandardPressurePa;
        // Уравнение состояния идеального газа: rho = P / (R * T).
        const double densityKgPerM3 =
            (rSpecific > 0.0 && temperatureKelvin > 0.0)
                ? (pressurePa / (rSpecific * temperatureKelvin))
                : 0.0;

        layers[i].heightMeters = heightMidMeters;
        layers[i].thicknessMeters = layerThickness;
        layers[i].temperatureKelvin = temperatureKelvin;
        layers[i].pressureAtm = pressureAtm;
        layers[i].densityKgPerM3 = densityKgPerM3;
    }

    return layers;
}

double AtmosphericProfileInitializer::resolveLapseRateKPerM(const Settings &settings) const {
    if (settings.lapseRateKPerM > 0.0) {
        return settings.lapseRateKPerM;
    }
    if (settings.useDryAdiabatic && settings.gravityMps2 > 0.0) {
        // Сухой адиабатический градиент: Γ_d = g / Cp.
        return AtmosphericThermodynamics::dryAdiabaticLapseRate(composition_,
                                                                settings.gravityMps2);
    }
    return 0.0;
}

double AtmosphericProfileInitializer::scaleHeightMeters(double temperatureKelvin,
                                                        double gravityMps2) const {
    const double rSpecific = AtmosphericThermodynamics::specificGasConstant(composition_);
    if (rSpecific <= 0.0 || gravityMps2 <= 0.0 || temperatureKelvin <= 0.0) {
        return 0.0;
    }
    // Масштабная высота: H = R * T / g.
    return (rSpecific * temperatureKelvin) / gravityMps2;
}

double AtmosphericProfileInitializer::pressureAtmAtHeight(double surfacePressureAtm,
                                                          double heightMeters,
                                                          double scaleHeightMeters) const {
    if (surfacePressureAtm <= 0.0 || scaleHeightMeters <= 0.0) {
        return 0.0;
    }
    return surfacePressureAtm * qExp(-heightMeters / scaleHeightMeters);
}
