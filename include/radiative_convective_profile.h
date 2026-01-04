#pragma once

#include "atmosphere_model.h"

#include <QVector>

struct RadiativeConvectiveLayer {
    double pressureAtm = 0.0;
    double temperatureKelvin = 0.0;
    double densityKgPerM3 = 0.0;
    double thicknessMeters = 0.0;
    double opticalDepthLongwave = 0.0;
    double opticalDepthShortwave = 0.0;
};

class RadiativeConvectiveProfile {
public:
    RadiativeConvectiveProfile(const AtmosphereComposition &composition,
                               double surfacePressureAtm,
                               double surfaceTemperatureKelvin,
                               double surfaceGravity,
                               int layerCount = 32);

    const QVector<RadiativeConvectiveLayer> &layers() const;
    double totalOpticalDepthLongwave() const;
    double totalOpticalDepthShortwave() const;

private:
    void buildProfile();
    double kappaLongwave(double temperatureKelvin, double pressurePa) const;
    double kappaShortwave(double temperatureKelvin, double pressurePa) const;
    double radiativeGradientDlnT(double temperatureKelvin, double pressurePa) const;
    double adiabaticGradientDlnT(double temperatureKelvin, double pressurePa) const;
    double mixedOpacity(double temperatureKelvin,
                        double pressurePa,
                        bool shortwave) const;

    AtmosphereComposition composition_;
    double surfacePressureAtm_ = 0.0;
    double surfaceTemperatureKelvin_ = 0.0;
    double surfaceGravity_ = 0.0;
    int layerCount_ = 0;
    QVector<RadiativeConvectiveLayer> layers_;
    double totalOpticalDepthLongwave_ = 0.0;
    double totalOpticalDepthShortwave_ = 0.0;
};
