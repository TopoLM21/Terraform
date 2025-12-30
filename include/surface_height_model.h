#pragma once

#include "planet_presets.h"
#include "surface_heightmap.h"

#include <QVector3D>

#include <array>

class SurfaceHeightModel {
public:
    static constexpr int kContinentCount = 5;

    SurfaceHeightModel();
    SurfaceHeightModel(HeightSourceType sourceType,
                       const QString &heightmapPath,
                       double heightmapScaleKm,
                       quint32 heightSeed,
                       bool useContinentsHeight);
    double heightKmAt(double latitudeDeg, double longitudeDeg) const;

private:
    void rebuildContinentCentersCache() const;

    HeightSourceType sourceType_ = HeightSourceType::Procedural;
    SurfaceHeightmap heightmap_;
    quint32 heightSeed_ = 0;
    bool useContinentsHeight_ = false;
    mutable std::array<QVector3D, kContinentCount> continentCenters_{};
    mutable bool continentCentersValid_ = false;
    mutable quint32 continentCentersSeed_ = 0;
};
