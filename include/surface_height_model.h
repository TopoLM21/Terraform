#pragma once

#include "planet_presets.h"
#include "surface_heightmap.h"

class SurfaceHeightModel {
public:
    SurfaceHeightModel();
    SurfaceHeightModel(HeightSourceType sourceType,
                       const QString &heightmapPath,
                       double heightmapScaleKm);
    double heightKmAt(double latitudeDeg, double longitudeDeg) const;

private:
    HeightSourceType sourceType_ = HeightSourceType::Procedural;
    SurfaceHeightmap heightmap_;
};
