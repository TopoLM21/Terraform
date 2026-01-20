#include "atmosphere_step_solver.h"

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

void AtmosphereStepSolver::runLayeredStep(const LayeredStepInput &input) {
    const int pointCount = input.surfaceGrid.points().size();
    const int columnCount = input.atmosphereGrid.columnCount();
    const int processedCount = qMin(pointCount, columnCount);
    if (processedCount <= 0) {
        return;
    }
    int logPointIndex = input.logPointIndex;
    if (logPointIndex < 0 || logPointIndex >= processedCount) {
        // Если индекс не задан, пишем лог для первой ячейки: так проще сравнивать шаги.
        logPointIndex = 0;
    }

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

        const QVector<double> layerDeltas =
            radiationSolver_.solve(column,
                                   localInsolation,
                                   surfaceAlbedo,
                                   input.cloudShortwaveTransmission);
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
    advectionSolver_.advectLayerWinds(input.surfaceGrid,
                                      input.atmosphereGrid,
                                      dayLengthSeconds_,
                                      timeStepSeconds_,
                                      1);
}
