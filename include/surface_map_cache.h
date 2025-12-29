#pragma once

#include <QSize>
#include <QtGlobal>
#include <QVector>

class MollweideProjection;
class PlanetSurfaceGrid;

struct PixelWeights {
    QVector<int> indices;
    QVector<float> weights;
};

struct SurfaceMapCacheKey {
    QSize widgetSize;
    double renderScale = 1.0;
    int neighborCount = 0;
    double power = 0.0;

    bool operator==(const SurfaceMapCacheKey &other) const;
};

class SurfaceMapCache {
public:
    void clear();

    bool isValidFor(const SurfaceMapCacheKey &key, const PlanetSurfaceGrid *grid) const;
    void rebuild(const PlanetSurfaceGrid *grid,
                 const MollweideProjection *projection,
                 const SurfaceMapCacheKey &key);

    const QVector<PixelWeights> &pixelWeights() const;
    const QVector<quint8> &insideMask() const;
    QSize scaledSize() const;

private:
    SurfaceMapCacheKey key_;
    QSize scaledSize_;
    const PlanetSurfaceGrid *grid_ = nullptr;
    QVector<PixelWeights> pixelWeights_;
    QVector<quint8> insideMask_;
};
