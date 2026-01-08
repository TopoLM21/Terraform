#include "thermal_conductivity_model.h"

#include <QtCore/QtMath>

namespace {
constexpr double kMinConductivity = 1e-6;
} // namespace

double thermalConductivityAtTemperature(double baseConductivity,
                                        const ThermalConductivityModel &model,
                                        double temperatureKelvin) {
    double conductivity = baseConductivity;
    switch (model.type) {
    case ThermalConductivityModelType::Constant:
        break;
    case ThermalConductivityModelType::Linear:
        conductivity = baseConductivity +
            model.temperatureCoefficient * (temperatureKelvin - model.referenceTemperatureKelvin);
        break;
    }
    return qMax(kMinConductivity, conductivity);
}
