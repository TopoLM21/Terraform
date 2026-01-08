#pragma once

#include "atmosphere_model.h"

#include <memory>

// Режим радиационной модели атмосферы.
// Fast: однопотоковая аппроксимация с экспоненциальным законом Бугера-Ламберта.
// Layered: многослойная серая атмосфера (двухпотоковая аппроксимация/Eddington).
enum class RadiationModelType {
    Fast,
    Layered
};

class RadiationModel {
public:
    virtual ~RadiationModel() = default;

    // Коэффициент пропускания входящего (коротковолнового) излучения.
    virtual double incomingTransmission() const = 0;
    // Коэффициент пропускания исходящего (длинноволнового) излучения.
    virtual double outgoingTransmission() const = 0;
    // Эффективная оптическая толщина атмосферы в длинноволновом диапазоне.
    virtual double effectiveOpticalDepth() const = 0;
};

std::unique_ptr<RadiationModel> makeRadiationModel(const AtmosphereComposition &composition,
                                                   double pressureAtm,
                                                   double baseTemperatureKelvin,
                                                   double effectiveTemperatureKelvin,
                                                   double surfaceGravity,
                                                   RadiationModelType type);
