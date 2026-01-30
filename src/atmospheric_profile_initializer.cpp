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

    const double rSpecific = AtmosphericThermodynamics::specificGasConstant(composition_);
    const double lapseRateKPerM = resolveLapseRateKPerM(settings);
    const double scaleHeight =
        scaleHeightMeters(settings.surfaceTemperatureKelvin, settings.gravityMps2);
    double pressureStopAtm = 0.0;
    if (settings.minTopPressureAtm > 0.0) {
        pressureStopAtm = settings.minTopPressureAtm;
    }
    if (settings.minTopPressureFraction > 0.0 && settings.surfacePressureAtm > 0.0) {
        const double fractionPressureAtm =
            settings.surfacePressureAtm * settings.minTopPressureFraction;
        pressureStopAtm = (pressureStopAtm > 0.0)
            ? qMax(pressureStopAtm, fractionPressureAtm)
            : fractionPressureAtm;
    }

    double resolvedTopHeight = 0.0;
    if (pressureStopAtm > 0.0 &&
        settings.surfacePressureAtm > pressureStopAtm &&
        scaleHeight > 0.0) {
        // Верхняя граница задаётся давлением: прекращаем слои, когда P падает ниже порога,
        // чтобы моделировать только атмосферу с заметной массой и оптическим вкладом.
        resolvedTopHeight =
            scaleHeight * qLn(settings.surfacePressureAtm / pressureStopAtm);
    } else if (settings.topHeightMeters > 0.0) {
        resolvedTopHeight = settings.topHeightMeters;
    } else {
        resolvedTopHeight = (scaleHeight > 0.0) ? qMax(1000.0, scaleHeight * 6.0) : 0.0;
    }

    layers.reserve(settings.layerCount);

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
    double pressureAtmAtBottom = settings.surfacePressureAtm;
    for (int i = 0; i < settings.layerCount; ++i) {
        const double layerThickness = useNonUniformLayers
            ? ((i == 0) ? bottomLayerThickness : remainingLayerThickness)
            : uniformLayerThickness;
        const double heightMidMeters = bottomEdgeMeters + layerThickness * 0.5;
        bottomEdgeMeters += layerThickness;
        // Температура по линейному градиенту: T(z) = T0 - Γ * z.
        const double temperatureKelvin = (settings.surfaceTemperatureKelvin > 0.0)
            ? qMax(1.0, settings.surfaceTemperatureKelvin - lapseRateKPerM * heightMidMeters)
            : 0.0;
        const double scaleHeightLayer = (rSpecific > 0.0 &&
                                         temperatureKelvin > 0.0 &&
                                         settings.gravityMps2 > 0.0)
            ? (rSpecific * temperatureKelvin) / settings.gravityMps2
            : 0.0;
        double pressureAtm = 0.0;
        if (pressureAtmAtBottom > 0.0 && scaleHeightLayer > 0.0) {
            // Гидростатическое равновесие с переменной температурой.
            pressureAtm = pressureAtmAtBottom *
                qExp(-0.5 * layerThickness / scaleHeightLayer);
            pressureAtmAtBottom =
                pressureAtmAtBottom * qExp(-layerThickness / scaleHeightLayer);
        }
        if (pressureStopAtm > 0.0 && pressureAtm > 0.0 && pressureAtm <= pressureStopAtm) {
            // Слой удаляется, если давление упало ниже порога — чтобы не плодить
            // «мертвые» ячейки с нулевой теплоёмкостью и фиксированной температурой.
            break;
        }
        if (pressureStopAtm > 0.0 && pressureAtm > 0.0) {
            // Минимальный физический порог для защиты от численного ухода давления в ноль,
            // иначе теплоёмкость слоя станет нулевой.
            pressureAtm = qMax(pressureAtm, pressureStopAtm);
        }
        const double pressurePa = pressureAtm * kStandardPressurePa;
        // Уравнение состояния идеального газа: rho = P / (R * T).
        const double densityKgPerM3 =
            (rSpecific > 0.0 && temperatureKelvin > 0.0)
                ? (pressurePa / (rSpecific * temperatureKelvin))
                : 0.0;

        AtmosphericProfileLayer layer;
        layer.heightMeters = heightMidMeters;
        layer.thicknessMeters = layerThickness;
        layer.temperatureKelvin = temperatureKelvin;
        layer.pressureAtm = pressureAtm;
        layer.densityKgPerM3 = densityKgPerM3;
        layers.append(layer);
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
