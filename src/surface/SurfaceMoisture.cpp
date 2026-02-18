#include "surface/SurfaceMoisture.h"

#include <QtCore/QtMath>

SurfaceMoisture::SurfaceMoisture(const Settings &settings)
    : settings_(settings) {}

const SurfaceMoisture::Settings &SurfaceMoisture::settings() const {
    return settings_;
}

void SurfaceMoisture::setSettings(const Settings &settings) {
    settings_ = settings;
    waterKgPerM2_ = clampedWater(waterKgPerM2_);
}

double SurfaceMoisture::waterKgPerM2() const {
    return waterKgPerM2_;
}

void SurfaceMoisture::setWaterKgPerM2(double waterKgPerM2) {
    waterKgPerM2_ = clampedWater(waterKgPerM2);
}

double SurfaceMoisture::addPrecipitation(double precipitationKgPerM2) {
    if (precipitationKgPerM2 <= 0.0) {
        return 0.0;
    }
    const double before = waterKgPerM2_;
    waterKgPerM2_ = clampedWater(waterKgPerM2_ + precipitationKgPerM2);
    // Разница между приходом и накоплением — это сток или инфильтрация.
    const double stored = waterKgPerM2_ - before;
    return qMax(0.0, precipitationKgPerM2 - stored);
}

double SurfaceMoisture::applyEvaporation(double potentialEvaporationKgPerM2) {
    if (potentialEvaporationKgPerM2 <= 0.0 || waterKgPerM2_ <= 0.0) {
        return 0.0;
    }
    const double evaporation = qMin(waterKgPerM2_, potentialEvaporationKgPerM2);
    waterKgPerM2_ = clampedWater(waterKgPerM2_ - evaporation);
    return evaporation;
}

double SurfaceMoisture::update(double precipitationKgPerM2,
                               double potentialEvaporationKgPerM2) {
    addPrecipitation(precipitationKgPerM2);
    return applyEvaporation(potentialEvaporationKgPerM2);
}

double SurfaceMoisture::moistureFraction() const {
    if (settings_.maxStorageKgPerM2 <= 0.0) {
        return 0.0;
    }
    return qBound(0.0, waterKgPerM2_ / settings_.maxStorageKgPerM2, 1.0);
}

double SurfaceMoisture::takeWaterKgPerM2(double requestedKgPerM2) {
    if (requestedKgPerM2 <= 0.0) {
        return 0.0;
    }
    const double delivered = qMin(requestedKgPerM2, waterKgPerM2_);
    waterKgPerM2_ = clampedWater(waterKgPerM2_ - delivered);
    return delivered;
}

void SurfaceMoisture::returnWaterKgPerM2(double returnedKgPerM2) {
    if (returnedKgPerM2 <= 0.0) {
        return;
    }
    waterKgPerM2_ = clampedWater(waterKgPerM2_ + returnedKgPerM2);
}

double SurfaceMoisture::clampedWater(double waterKgPerM2) const {
    // Нет верхнего предела: избыток стекает через дренажную фазу.
    // maxStorageKgPerM2 используется только в moistureFraction() (полевая влагоёмкость).
    return qMax(0.0, waterKgPerM2);
}
