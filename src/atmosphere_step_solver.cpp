#include "atmosphere_step_solver.h"

#include "atmospheric_thermodynamics.h"
#include "surface_atmosphere_coupler.h"

#include <QtCore/QLoggingCategory>
#include <QtCore/QtMath>

Q_LOGGING_CATEGORY(atmosphereProfileLog, "solar.atmosphere.profile")

AtmosphereStepSolver::AtmosphereStepSolver(const AtmosphereComposition &composition,
                                           double gravityMps2,
                                           double timeStepSeconds,
                                           double dayLengthSeconds,
                                           bool isRetrograde)
    : radiationSolver_(timeStepSeconds)
    , convectiveSolver_(composition, gravityMps2)
    , gravityMps2_(gravityMps2)
    , rSpecific_(AtmosphericThermodynamics::specificGasConstant(composition))
    , specificHeatCp_(AtmosphericThermodynamics::specificHeatCp(composition))
    , timeStepSeconds_(timeStepSeconds)
    , dayLengthSeconds_(dayLengthSeconds)
    , isRetrograde_(isRetrograde) {}

const SurfaceMaterial &AtmosphereStepSolver::materialForPoint(
    const QHash<QString, SurfaceMaterial> &materialsById,
    const SurfaceMaterial &defaultMaterial,
    const QString &materialId) const {
    const auto it = materialsById.constFind(materialId);
    return it != materialsById.cend() ? *it : defaultMaterial;
}

void AtmosphereStepSolver::updateLayerPressures(double surfacePressureAtm,
                                                QVector<AtmosphericLayerState> &layers) const {
    if (surfacePressureAtm <= 0.0 || rSpecific_ <= 0.0 || gravityMps2_ <= 0.0) {
        return;
    }

    double pressureAtmAtBottom = surfacePressureAtm;
    for (int layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
        const double layerTemperatureKelvin = layers.at(layerIndex).temperatureKelvin();
        const double layerThicknessMeters = layers.at(layerIndex).thicknessMeters();
        double layerPressureAtm = 0.0;
        if (pressureAtmAtBottom > 0.0 && layerTemperatureKelvin > 0.0 &&
            layerThicknessMeters > 0.0) {
            // Масштабная высота: H = R_specific * T / g, где
            // R_specific [Дж/(кг·К)], T [К], g [м/с²], H [м].
            const double scaleHeightMeters =
                (rSpecific_ * layerTemperatureKelvin) / gravityMps2_;
            if (scaleHeightMeters > 0.0) {
                // Гидростатика: dP/dz = -P/H ⇒ P(z+dz) = P(z) * exp(-dz/H).
                // Давление в середине слоя: P_mid = P_bottom * exp(-0.5 * dz / H) [атм].
                layerPressureAtm =
                    pressureAtmAtBottom * qExp(-0.5 * layerThicknessMeters / scaleHeightMeters);
                pressureAtmAtBottom =
                    pressureAtmAtBottom * qExp(-layerThicknessMeters / scaleHeightMeters);
            }
        }
        layers[layerIndex].setPressureAtm(layerPressureAtm);

        // Согласуем P–T–rho–C (давление–температура–плотность–теплоёмкость),
        // чтобы термодинамика оставалась устойчивой при изменении профиля.
        double densityKgPerM3 = 0.0;
        if (layerPressureAtm > 0.0 && layerTemperatureKelvin > 0.0 && rSpecific_ > 0.0) {
            const double pressurePa = layerPressureAtm * 101325.0;
            // Уравнение состояния идеального газа: rho = P / (R_specific * T).
            densityKgPerM3 = pressurePa / (rSpecific_ * layerTemperatureKelvin);
        }
        layers[layerIndex].setDensityKgPerM3(densityKgPerM3);

        const double heatCapacityJPerM2K =
            (densityKgPerM3 > 0.0 && specificHeatCp_ > 0.0 && layerThicknessMeters > 0.0)
                ? densityKgPerM3 * specificHeatCp_ * layerThicknessMeters
                : 0.0;
        layers[layerIndex].setHeatCapacityJPerM2K(heatCapacityJPerM2K);
    }
}

void AtmosphereStepSolver::runLayeredStep(const LayeredStepInput &input) {
    const int pointCount = input.surfaceGrid.points().size();
    const int columnCount = input.atmosphereGrid.columnCount();
    const int processedCount = qMin(pointCount, columnCount);
    if (processedCount <= 0) {
        return;
    }
    verticalWindMixingSolver_.setMixingCoefficient(input.verticalWindMixingCoefficientKz);
    int logPointIndex = input.logPointIndex;
    if (logPointIndex < 0 || logPointIndex >= processedCount) {
        // Если индекс не задан, пишем лог для первой ячейки: так проще сравнивать шаги.
        logPointIndex = 0;
    }

    const double minTopPressureAtm = input.atmosphereGrid.minTopPressureAtm();

    // Последовательность шага атмосферы:
    // (1) радиация по слоям, (2) конвективная коррекция,
    // (3) обмен с поверхностью, (4) обновление ветра, (5) горизонтальная адвекция.
    for (int i = 0; i < processedCount; ++i) {
        auto &point = input.surfaceGrid.points()[i];
        AtmosphericColumn &column = input.atmosphereGrid.columns()[i];
        const double localInsolation =
            (i < input.localInsolations.size()) ? input.localInsolations.at(i) : 0.0;
        const SurfaceMaterial &material =
            materialForPoint(input.materialsById, input.defaultMaterial, point.materialId);
        const double surfaceAlbedo = qBound(0.0, material.albedo, 1.0);

        // Поверхность и нижний слой — разные сущности: передаём температуру поверхности явно,
        // чтобы не подменять её температурой слоя и не получать самоподогрев атмосферы.
        const QVector<double> layerDeltas =
            radiationSolver_.solve(column,
                                   localInsolation,
                                   surfaceAlbedo,
                                   input.cloudShortwaveTransmission,
                                   point.state.temperatureKelvin());
        auto &layers = column.layers();
        const int layerCount = qMin(layers.size(), layerDeltas.size());
        for (int layerIndex = 0; layerIndex < layerCount; ++layerIndex) {
            const double updatedTemperature =
                qMax(0.0, layers.at(layerIndex).temperatureKelvin() + layerDeltas.at(layerIndex));
            layers[layerIndex].setTemperatureKelvin(updatedTemperature);
        }

        convectiveSolver_.adjust(column);

        const double initialAirTemperature =
            (point.airTemperatureK > 0.0) ? point.airTemperatureK : point.temperatureK;
        if (layers.isEmpty()) {
            point.airTemperatureK = initialAirTemperature;
            continue;
        }

        point.airTemperatureK = layers.first().temperatureKelvin();

        SurfaceAtmosphereCoupler coupler(input.heatTransferCoefficientWPerM2K);
        double surfaceAirFluxWPerM2 = 0.0;
        coupler.exchangeHeat(point.state,
                             layers[0],
                             material.roughnessLengthMeters,
                             timeStepSeconds_,
                             &surfaceAirFluxWPerM2);
        point.airTemperatureK = layers[0].temperatureKelvin();
        point.temperatureK = point.state.temperatureKelvin();
        point.surfaceAirFluxWPerM2 = surfaceAirFluxWPerM2;
        // Поток в грунт учитывает радиацию и обмен с воздухом.
        point.subsurfaceFluxWPerM2 += -surfaceAirFluxWPerM2;

        updateLayerPressures(point.pressureAtm, layers);

        if (minTopPressureAtm > 0.0 && !layers.isEmpty()) {
            const AtmosphericLayerState &topLayer = layers.last();
            const double fixedLayerThicknessMeters = topLayer.thicknessMeters();
            if (fixedLayerThicknessMeters > 0.0 &&
                topLayer.temperatureKelvin() > 0.0 &&
                topLayer.pressureAtm() > 0.0 &&
                gravityMps2_ > 0.0 &&
                rSpecific_ > 0.0) {
                // Масштабная высота: H = R_specific * T / g, тогда
                // P(z + dz) = P(z) * exp(-dz / H). Для нового слоя берём dz = толщину слоя.
                // Толщины фиксируем, чтобы новые слои добавлялись без перерасчёта геометрии.
                const double scaleHeightMeters =
                    (rSpecific_ * topLayer.temperatureKelvin()) / gravityMps2_;
                if (scaleHeightMeters > 0.0) {
                    const double nextLayerPressureAtm =
                        topLayer.pressureAtm() *
                        qExp(-fixedLayerThicknessMeters / scaleHeightMeters);
                    if (nextLayerPressureAtm > minTopPressureAtm) {
                        input.atmosphereGrid.updateColumnLayerCountFixedThickness(
                            i,
                            layers.size() + 1,
                            fixedLayerThicknessMeters);
                    } else if (topLayer.pressureAtm() < minTopPressureAtm &&
                               layers.size() > 1) {
                        input.atmosphereGrid.updateColumnLayerCountFixedThickness(
                            i,
                            layers.size() - 1,
                            fixedLayerThicknessMeters);
                    }
                }
            }
        }

        if (i == logPointIndex) {
            qCInfo(atmosphereProfileLog) << "Atmosphere profile (layered step)"
                                         << "index=" << i
                                         << "layerCount=" << layers.size();
            for (int layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
                const AtmosphericLayerState &layer = layers.at(layerIndex);
                qCInfo(atmosphereProfileLog) << "Layer"
                                             << layerIndex
                                             << "heightKm=" << layer.heightMeters() / 1000.0
                                             << "temperatureK=" << layer.temperatureKelvin()
                                             << "pressureAtm=" << layer.pressureAtm();
            }
        }
    }

    dynamicsSolver_.updateLayerWinds(input.surfaceGrid,
                                    input.atmosphereGrid,
                                    dayLengthSeconds_,
                                    isRetrograde_,
                                    timeStepSeconds_,
                                    1);
    for (int i = 0; i < processedCount; ++i) {
        verticalWindMixingSolver_.mix(input.atmosphereGrid.columns()[i], timeStepSeconds_);
    }
    advectionSolver_.advectLayerWinds(input.surfaceGrid,
                                      input.atmosphereGrid,
                                      dayLengthSeconds_,
                                      timeStepSeconds_,
                                      1);
}
