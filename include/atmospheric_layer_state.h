#pragma once

class AtmosphericLayerState {
public:
    AtmosphericLayerState() = default;
    AtmosphericLayerState(double temperatureKelvin,
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
                          double relativeHumidity);

    double temperatureKelvin() const;
    void setTemperatureKelvin(double temperatureKelvin);

    double pressureAtm() const;
    void setPressureAtm(double pressureAtm);

    double densityKgPerM3() const;
    void setDensityKgPerM3(double densityKgPerM3);

    double heightMeters() const;
    void setHeightMeters(double heightMeters);

    double thicknessMeters() const;
    void setThicknessMeters(double thicknessMeters);

    double heatCapacityJPerM2K() const;
    void setHeatCapacityJPerM2K(double heatCapacityJPerM2K);

    double windUMps() const;
    void setWindUMps(double windUMps);

    double windVMps() const;
    void setWindVMps(double windVMps);

    double opticalDepthShortwave() const;
    void setOpticalDepthShortwave(double opticalDepthShortwave);

    double opticalDepthLongwave() const;
    void setOpticalDepthLongwave(double opticalDepthLongwave);

    bool convectionEnabled() const;
    void setConvectionEnabled(bool convectionEnabled);

    double convectionMixingCoefficient() const;
    void setConvectionMixingCoefficient(double convectionMixingCoefficient);

    double waterVaporKgPerM2() const;
    void setWaterVaporKgPerM2(double waterVaporKgPerM2);

    double liquidWaterKgPerM2() const;
    void setLiquidWaterKgPerM2(double liquidWaterKgPerM2);

    double iceWaterKgPerM2() const;
    void setIceWaterKgPerM2(double iceWaterKgPerM2);

    double relativeHumidity() const;
    void setRelativeHumidity(double relativeHumidity);

private:
    double temperatureKelvin_ = 0.0;
    double pressureAtm_ = 0.0;
    double densityKgPerM3_ = 0.0;
    double heightMeters_ = 0.0;
    double thicknessMeters_ = 0.0;
    double heatCapacityJPerM2K_ = 0.0;
    double windUMps_ = 0.0;
    double windVMps_ = 0.0;
    // Оптическая толщина слоя (τ) для коротковолнового излучения.
    double opticalDepthShortwave_ = 0.0;
    // Оптическая толщина слоя (τ) для длинноволнового излучения.
    double opticalDepthLongwave_ = 0.0;
    bool convectionEnabled_ = false;
    double convectionMixingCoefficient_ = 0.0;
    // Масса водяного пара в слое, кг/м².
    double waterVaporKgPerM2_ = 0.0;
    // Масса жидкой воды (конденсата) в слое, кг/м².
    double liquidWaterKgPerM2_ = 0.0;
    // Масса льда (конденсата) в слое, кг/м².
    double iceWaterKgPerM2_ = 0.0;
    // Относительная влажность слоя (0..1), вычисляется по насыщению.
    double relativeHumidity_ = 0.0;
};
