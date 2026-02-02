#include "atmospheric_layer_state.h"

AtmosphericLayerState::AtmosphericLayerState(double temperatureKelvin,
                                             double pressureAtm,
                                             double densityKgPerM3,
                                             double heightMeters,
                                             double thicknessMeters,
                                             double heatCapacityJPerM2K,
                                             double windUMps,
                                             double windVMps,
                                             double opticalDepthShortwave,
                                             double opticalDepthLongwave,
                                             bool convectionEnabled,
                                             double convectionMixingCoefficient,
                                             double waterVaporKgPerM2,
                                             double liquidWaterKgPerM2,
                                             double iceWaterKgPerM2,
                                             double relativeHumidity)
    : temperatureKelvin_(temperatureKelvin)
    , pressureAtm_(pressureAtm)
    , densityKgPerM3_(densityKgPerM3)
    , heightMeters_(heightMeters)
    , thicknessMeters_(thicknessMeters)
    , heatCapacityJPerM2K_(heatCapacityJPerM2K)
    , windUMps_(windUMps)
    , windVMps_(windVMps)
    , opticalDepthShortwave_(opticalDepthShortwave)
    , opticalDepthLongwave_(opticalDepthLongwave)
    , convectionEnabled_(convectionEnabled)
    , convectionMixingCoefficient_(convectionMixingCoefficient)
    , waterVaporKgPerM2_(waterVaporKgPerM2)
    , liquidWaterKgPerM2_(liquidWaterKgPerM2)
    , iceWaterKgPerM2_(iceWaterKgPerM2)
    , relativeHumidity_(relativeHumidity) {}

double AtmosphericLayerState::temperatureKelvin() const {
    return temperatureKelvin_;
}

void AtmosphericLayerState::setTemperatureKelvin(double temperatureKelvin) {
    temperatureKelvin_ = temperatureKelvin;
}

double AtmosphericLayerState::pressureAtm() const {
    return pressureAtm_;
}

void AtmosphericLayerState::setPressureAtm(double pressureAtm) {
    pressureAtm_ = pressureAtm;
}

double AtmosphericLayerState::densityKgPerM3() const {
    return densityKgPerM3_;
}

void AtmosphericLayerState::setDensityKgPerM3(double densityKgPerM3) {
    densityKgPerM3_ = densityKgPerM3;
}

double AtmosphericLayerState::heightMeters() const {
    return heightMeters_;
}

void AtmosphericLayerState::setHeightMeters(double heightMeters) {
    heightMeters_ = heightMeters;
}

double AtmosphericLayerState::thicknessMeters() const {
    return thicknessMeters_;
}

void AtmosphericLayerState::setThicknessMeters(double thicknessMeters) {
    thicknessMeters_ = thicknessMeters;
}

double AtmosphericLayerState::heatCapacityJPerM2K() const {
    return heatCapacityJPerM2K_;
}

void AtmosphericLayerState::setHeatCapacityJPerM2K(double heatCapacityJPerM2K) {
    heatCapacityJPerM2K_ = heatCapacityJPerM2K;
}

double AtmosphericLayerState::windUMps() const {
    return windUMps_;
}

void AtmosphericLayerState::setWindUMps(double windUMps) {
    windUMps_ = windUMps;
}

double AtmosphericLayerState::windVMps() const {
    return windVMps_;
}

void AtmosphericLayerState::setWindVMps(double windVMps) {
    windVMps_ = windVMps;
}

double AtmosphericLayerState::opticalDepthShortwave() const {
    return opticalDepthShortwave_;
}

void AtmosphericLayerState::setOpticalDepthShortwave(double opticalDepthShortwave) {
    opticalDepthShortwave_ = opticalDepthShortwave;
}

double AtmosphericLayerState::opticalDepthLongwave() const {
    return opticalDepthLongwave_;
}

void AtmosphericLayerState::setOpticalDepthLongwave(double opticalDepthLongwave) {
    opticalDepthLongwave_ = opticalDepthLongwave;
}

bool AtmosphericLayerState::convectionEnabled() const {
    return convectionEnabled_;
}

void AtmosphericLayerState::setConvectionEnabled(bool convectionEnabled) {
    convectionEnabled_ = convectionEnabled;
}

double AtmosphericLayerState::convectionMixingCoefficient() const {
    return convectionMixingCoefficient_;
}

void AtmosphericLayerState::setConvectionMixingCoefficient(double convectionMixingCoefficient) {
    convectionMixingCoefficient_ = convectionMixingCoefficient;
}

double AtmosphericLayerState::waterVaporKgPerM2() const {
    return waterVaporKgPerM2_;
}

void AtmosphericLayerState::setWaterVaporKgPerM2(double waterVaporKgPerM2) {
    waterVaporKgPerM2_ = waterVaporKgPerM2;
}

double AtmosphericLayerState::liquidWaterKgPerM2() const {
    return liquidWaterKgPerM2_;
}

void AtmosphericLayerState::setLiquidWaterKgPerM2(double liquidWaterKgPerM2) {
    liquidWaterKgPerM2_ = liquidWaterKgPerM2;
}

double AtmosphericLayerState::iceWaterKgPerM2() const {
    return iceWaterKgPerM2_;
}

void AtmosphericLayerState::setIceWaterKgPerM2(double iceWaterKgPerM2) {
    iceWaterKgPerM2_ = iceWaterKgPerM2;
}

double AtmosphericLayerState::relativeHumidity() const {
    return relativeHumidity_;
}

void AtmosphericLayerState::setRelativeHumidity(double relativeHumidity) {
    relativeHumidity_ = relativeHumidity;
}
