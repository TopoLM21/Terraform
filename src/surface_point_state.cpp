#include "surface_point_state.h"

#include "emission_layer_model.h"
#include "radiation_model_utils.h"
#include "planet_presets.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kWaterFreezeTemperatureKelvin = 273.15;
// Параметры льда: λ ≈ 2.2 Вт/(м·К), ρ ≈ 920 кг/м³, c ≈ 2100 Дж/(кг·К).
constexpr double kIceThermalConductivity = 2.2;
constexpr double kIceDensity = 920.0;
constexpr double kIceSpecificHeat = 2100.0;
// Параметры воды: λ ≈ 0.6 Вт/(м·К), ρ ≈ 1000 кг/м³, c ≈ 4180 Дж/(кг·К).
constexpr double kWaterThermalConductivity = 0.6;
constexpr double kWaterDensity = 1000.0;
constexpr double kWaterSpecificHeat = 4180.0;
}

SurfacePointState::SurfacePointState(double initialTemperatureKelvin,
                                     double albedo,
                                     double greenhouseOpacity,
                                     RadiationModelType radiationModelType,
                                     double minTemperatureKelvin,
                                     const SurfaceMaterial &material,
                                     const SubsurfaceModelSettings &subsurfaceSettings)
    : albedo_(qBound(0.0, albedo, 1.0)),
      greenhouseOpacity_(qBound(0.0, greenhouseOpacity, 0.999)),
      emissivity_(qBound(0.0, material.emissivity, 1.0)),
      minTemperatureKelvin_(qMax(0.0, minTemperatureKelvin)),
      radiationModelType_(radiationModelType) {
    const SubsurfaceGrid grid(subsurfaceSettings.layerCount,
                              subsurfaceSettings.topLayerThicknessMeters,
                              subsurfaceSettings.bottomDepthMeters);
    solver_.reset(grid,
                  material.thermalConductivity,
                  material.thermalConductivityModel,
                  material.density,
                  material.specificHeat,
                  subsurfaceSettings.bottomBoundary,
                  initialTemperatureKelvin,
                  subsurfaceSettings.geothermalFluxWPerM2);
    if (material.id == QStringLiteral("ocean")) {
        const QVector<double> layerDepths = grid.layerDepthsMeters();
        const double phaseBoundaryDepth =
            qMax(subsurfaceSettings.topLayerThicknessMeters,
                 subsurfaceSettings.oceanPhaseBoundaryDepthMeters);
        // Верхние слои океана (по глубине центров слоёв) могут замерзать при T < 273.15 K.
        // Граница задаётся как характерная глубина сезонного промерзания: ниже неё
        // считаем воду всегда жидкой, чтобы не фиксировать весь столбик как лёд.
        solver_.setLayerPropertiesProvider(
            [layerDepths, phaseBoundaryDepth](int layerIndex, double temperatureKelvin) {
                const bool isUpperLayer =
                    (layerIndex >= 0 && layerIndex < layerDepths.size())
                        ? (layerDepths.at(layerIndex) <= phaseBoundaryDepth)
                        : true;
                const bool isIce = isUpperLayer && temperatureKelvin < kWaterFreezeTemperatureKelvin;
                if (isIce) {
                    return SubsurfaceLayerProperties{
                        kIceThermalConductivity,
                        {ThermalConductivityModelType::Constant, kWaterFreezeTemperatureKelvin, 0.0},
                        kIceDensity,
                        kIceSpecificHeat};
                }
                return SubsurfaceLayerProperties{
                    kWaterThermalConductivity,
                    {ThermalConductivityModelType::Constant, kWaterFreezeTemperatureKelvin, 0.0},
                    kWaterDensity,
                    kWaterSpecificHeat};
            });
    }
    solver_.setInitialTemperature(qMax(minTemperatureKelvin_, initialTemperatureKelvin));
}

double SurfacePointState::temperatureKelvin() const {
    return solver_.surfaceTemperatureKelvin();
}

double SurfacePointState::greenhouseOpacity() const {
    return greenhouseOpacity_;
}

double SurfacePointState::emissivity() const {
    return emissivity_;
}

void SurfacePointState::setTemperatureKelvin(double temperatureKelvin) {
    solver_.setSurfaceLayerTemperatureKelvin(qMax(minTemperatureKelvin_, temperatureKelvin));
}

void SurfacePointState::setProfileTemperatureKelvin(double temperatureKelvin) {
    solver_.setInitialTemperature(qMax(minTemperatureKelvin_, temperatureKelvin));
}

void SurfacePointState::setSurfaceLayerTemperatureKelvin(double temperatureKelvin) {
    solver_.setSurfaceLayerTemperatureKelvin(qMax(minTemperatureKelvin_, temperatureKelvin));
}

void SurfacePointState::setGreenhouseOpacity(double greenhouseOpacity) {
    greenhouseOpacity_ = qBound(0.0, greenhouseOpacity, 0.999);
}

void SurfacePointState::setAlbedo(double albedo) {
    albedo_ = qBound(0.0, albedo, 1.0);
}

const SubsurfaceTemperatureSolver &SurfacePointState::solver() const {
    return solver_;
}

double SurfacePointState::absorbedFlux(double solarIrradiance) const {
    return solarIrradiance * (1.0 - albedo_);
}

double SurfacePointState::surfaceEmittedFlux() const {
    const double surfaceTemperature = temperatureKelvin();
    // Локальный баланс поверхности: поток излучения задаётся T_surface.
    return emissivity_ * kStefanBoltzmannConstant * std::pow(surfaceTemperature, 4.0);
}

double SurfacePointState::toaEmittedFlux() const {
    // emissionLayerTemperature применяется только для оценки потока на ТОА,
    // а не для локального баланса поверхности.
    const double tauSurface =
        opticalDepthFromGreenhouseOpacity(greenhouseOpacity_, radiationModelType_);
    const double emissionTemperature =
        emissionLayerTemperature(temperatureKelvin(), tauSurface);
    return emissivity_ * kStefanBoltzmannConstant * std::pow(emissionTemperature, 4.0);
}

double SurfacePointState::topLayerHeatCapacityJPerM2K() const {
    return solver_.topLayerHeatCapacity();
}

double SurfacePointState::stabilizedRadiativeFlux(double absorbedFlux,
                                                  double dtSeconds) const {
    const double surfaceTemperature = temperatureKelvin();
    const double emittedNow =
        emissivity_ * kStefanBoltzmannConstant * std::pow(surfaceTemperature, 4.0);
    // Линеаризованный радиационный поток:
    // F_rad(T) ≈ F_rad(T_n) + (dF/dT)|_{T_n} (T - T_n),
    // dF/dT = 4 * sigma * T_surface^3.
    // Тогда неявная стабилизация даёт
    // F_eff = (F_in - F_rad(T_n)) / (1 + alpha), alpha = (dF/dT) * dt / C,
    // где C — теплоёмкость верхнего слоя в J/(м²·K).
    const double heatCapacity = solver_.topLayerHeatCapacity();
    const double radiativeDerivative =
        4.0 * emissivity_ * kStefanBoltzmannConstant * std::pow(surfaceTemperature, 3.0);
    const double alpha = (heatCapacity > 0.0)
        ? (radiativeDerivative * dtSeconds / heatCapacity)
        : 0.0;
    return (absorbedFlux - emittedNow) / (1.0 + qMax(0.0, alpha));
}

void SurfacePointState::updateTemperature(double absorbedFlux,
                                          double emittedFlux,
                                          double dtSeconds) {
    Q_UNUSED(emittedFlux)
    // Обновляем температуру верхнего слоя с учётом стабилизированного
    // радиационного потока (W/м²).
    const double netFlux = stabilizedRadiativeFlux(absorbedFlux, dtSeconds);
    solver_.stepImplicit(netFlux, dtSeconds);
    clampProfile();
}

void SurfacePointState::applySurfaceFlux(double netFlux, double dtSeconds) {
    solver_.stepImplicit(netFlux, dtSeconds);
    clampProfile();
}

void SurfacePointState::clampProfile() {
    auto profile = solver_.temperatures();
    for (double &value : profile) {
        value = qMax(minTemperatureKelvin_, value);
    }
    solver_.setTemperatures(profile);
}
