#pragma once

#include <QtCore/QtGlobal>

struct GeothermalBodyParameters {
    double massKg = 0.0;
    double radiusM = 0.0;
    double ageGyr = 0.0;
    // Масштаб радиогенного тепловыделения относительно земного эталона.
    double radiogenicScale = 1.0;
    // Масштаб остаточного тепла относительно земного эталона.
    double residualScale = 1.0;
    // Характерное время распада радиогенного тепла (в млрд лет).
    double radiogenicDecayGyr = 8.0;
    // Характерное время остывания остаточного тепла (в млрд лет).
    double residualDecayGyr = 4.0;
    // Параметры приливного разогрева.
    double parentMassKg = 0.0;
    double semiMajorAxisM = 0.0;
    double eccentricity = 0.0;
    // k2/Q: число Лава, деленное на добротность диссипации.
    double k2OverQ = 0.0;
};

class GeothermalFluxModel {
public:
    double surfaceFluxWPerM2(const GeothermalBodyParameters &parameters) const;

private:
    double radiogenicPowerW(const GeothermalBodyParameters &parameters) const;
    double residualPowerW(const GeothermalBodyParameters &parameters) const;
    double tidalPowerW(const GeothermalBodyParameters &parameters) const;
};
