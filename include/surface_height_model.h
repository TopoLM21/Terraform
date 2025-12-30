#pragma once

#include "planet_presets.h"
#include "surface_heightmap.h"

class SurfaceHeightModel {
public:
    SurfaceHeightModel();
    SurfaceHeightModel(HeightSourceType sourceType,
                       const QString &heightmapPath,
                       double heightmapScaleKm,
                       quint32 heightSeed,
                       bool useContinentsHeight);
    double heightKmAt(double latitudeDeg, double longitudeDeg) const;

private:
    HeightSourceType sourceType_ = HeightSourceType::Procedural;
    SurfaceHeightmap heightmap_;
    quint32 heightSeed_ = 0;
    bool useContinentsHeight_ = false;
};
