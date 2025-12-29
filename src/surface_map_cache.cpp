#include "surface_map_cache.h"

#include "mollweide_projection.h"
#include "planet_surface_grid.h"
#include "surface_map_interpolator.h"

#include <QtMath>

namespace {
bool fuzzyCompare(double left, double right) {
    return qFuzzyCompare(left + 1.0, right + 1.0);
}

QSize scaledSizeForKey(const SurfaceMapCacheKey &key) {
    const int width = qMax(1, qRound(key.widgetSize.width() * key.renderScale));
    const int height = qMax(1, qRound(key.widgetSize.height() * key.renderScale));
    return QSize(width, height);
}
} // namespace

bool SurfaceMapCacheKey::operator==(const SurfaceMapCacheKey &other) const {
    return widgetSize == other.widgetSize
        && neighborCount == other.neighborCount
        && fuzzyCompare(renderScale, other.renderScale)
        && fuzzyCompare(power, other.power);
}

void SurfaceMapCache::clear() {
    key_ = SurfaceMapCacheKey{};
    scaledSize_ = QSize();
    grid_ = nullptr;
    pixelWeights_.clear();
    insideMask_.clear();
}

bool SurfaceMapCache::isValidFor(const SurfaceMapCacheKey &key,
                                 const PlanetSurfaceGrid *grid) const {
    return grid_ == grid && key_ == key && !pixelWeights_.isEmpty();
}

void SurfaceMapCache::rebuild(const PlanetSurfaceGrid *grid,
                              const MollweideProjection *projection,
                              const SurfaceMapCacheKey &key) {
    clear();
    grid_ = grid;
    key_ = key;
    scaledSize_ = scaledSizeForKey(key);

    if (!grid_ || !projection || grid_->points().isEmpty() || scaledSize_.isEmpty()) {
        return;
    }

    SurfaceMapInterpolator interpolator(grid_, projection);
    pixelWeights_ = interpolator.buildPixelWeights(scaledSize_,
                                                   key_.neighborCount,
                                                   key_.power,
                                                   &insideMask_);
}

const QVector<PixelWeights> &SurfaceMapCache::pixelWeights() const {
    return pixelWeights_;
}

const QVector<quint8> &SurfaceMapCache::insideMask() const {
    return insideMask_;
}

QSize SurfaceMapCache::scaledSize() const {
    return scaledSize_;
}
