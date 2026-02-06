#pragma once

#include "surface_point_state.h"
#include "surface/SurfaceMoisture.h"
#include "fluids/PhaseModel.h"

#include <QVector>
#include <QString>

struct SurfacePoint {
    double latitudeDeg = 0.0;
    double longitudeDeg = 0.0;
    double latitudeRadians = 0.0;
    double longitudeRadians = 0.0;
    double sinLatitude = 0.0;
    double cosLatitude = 1.0;
    double radiusKm = 0.0;
    double temperatureK = 0.0;
    double airTemperatureK = 0.0;
    double solarFluxWPerM2 = 0.0;
    // Потоки в W/м²: shortwaveSurface — SW у поверхности,
    // longwaveUp/Down — LW вверх/вниз на границе поверхность↔атмосфера.
    // surfaceAirFlux положителен при переносе энергии от поверхности в воздух.
    // subsurfaceFlux положителен при уходе энергии в грунт (вниз).
    double shortwaveSurfaceWPerM2 = 0.0;
    double longwaveUpWPerM2 = 0.0;
    double longwaveDownWPerM2 = 0.0;
    double surfaceAirFluxWPerM2 = 0.0;
    double subsurfaceFluxWPerM2 = 0.0;
    // Осадки за шаг моделирования, кг/м² за текущий интервал времени.
    double precipitationRateKgPerM2 = 0.0;
    // Сглаженная (EMA) интенсивность осадков, кг/м² за интервал; это не интеграл за всё время.
    double precipitationKgPerM2 = 0.0;
    // Запас снега/льда на суше, кг/м² (водный эквивалент).
    double snowKgPerM2 = 0.0;
    SurfaceMoisture surfaceMoisture;
    // Состояние хранится отдельно для каждой точки, чтобы учитывать локальную
    // тепловую инерцию и парниковую поправку без глобального пересчёта.
    SurfacePointState state;
    double heightKm = 0.0;
    double pressureAtm = 0.0;
    double windEastMps = 0.0;
    double windNorthMps = 0.0;
    double windSpeedMps = 0.0;
    QVector<double> layerHeightsKm;
    QString materialId;
    // Фаза воды для океанических точек (используется для льда/альбедо/отрисовки).
    PhaseModel::Phase waterPhase = PhaseModel::Phase::Liquid;
    // Визуальная непрозрачность облаков над точкой (0..1), используется в отрисовке.
    double cloudOpacity = 0.0;
    // Доля площади, покрытой растительностью (0..1), обновляется моделью биома.
    double vegetationFraction = 0.0;
    // Биомасса растительности (кг/м²), учитывается в росте/диффузии.
    double vegetationBiomass = 0.0;
    // Маска/вес для визуального смешивания растительности (0..1).
    double vegetationBlendMask = 0.0;
    QVector<int> neighborIndices;
};
