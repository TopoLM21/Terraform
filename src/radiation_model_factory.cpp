#include "radiation_model.h"

#include "atmospheric_radiation_model.h"

std::unique_ptr<RadiationModel> makeRadiationModel(const AtmosphereComposition &composition,
                                                   double pressureAtm,
                                                   double baseTemperatureKelvin,
                                                   RadiationModelType type) {
    // Пока поддерживаем единственную модель, но выбор через фабрику
    // позволяет расширять набор реализаций без изменения кода расчёта.
    return std::make_unique<AtmosphericRadiationModel>(composition,
                                                       pressureAtm,
                                                       baseTemperatureKelvin,
                                                       type);
}
