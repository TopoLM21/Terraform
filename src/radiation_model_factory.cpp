#include "radiation_model.h"

#include "atmospheric_radiation_model.h"
#include "layered_radiation_model.h"

std::unique_ptr<RadiationModel> makeRadiationModel(const AtmosphereComposition &composition,
                                                   double pressureAtm,
                                                   double baseTemperatureKelvin,
                                                   double effectiveTemperatureKelvin,
                                                   double surfaceGravity,
                                                   RadiationModelType type) {
    // Fast оставляем однопоточной параметризацией, Layered — отдельной
    // многослойной реализацией для радиационно-конвективного профиля.
    if (type == RadiationModelType::Layered) {
        return std::make_unique<LayeredRadiationModel>(composition,
                                                       pressureAtm,
                                                       baseTemperatureKelvin,
                                                       effectiveTemperatureKelvin,
                                                       surfaceGravity);
    }
    return std::make_unique<AtmosphericRadiationModel>(composition,
                                                       pressureAtm,
                                                       baseTemperatureKelvin,
                                                       effectiveTemperatureKelvin,
                                                       surfaceGravity,
                                                       RadiationModelType::Fast);
}
