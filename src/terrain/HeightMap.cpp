#include "terrain/HeightMap.h"

#include <algorithm>

HeightMap::HeightMap() = default;

HeightMap::HeightMap(int width, int height, double cellAreaSquareMeters)
    : m_width(width)
    , m_height(height)
    , m_cellAreaSquareMeters(cellAreaSquareMeters) {
    m_heightsMeters.resize(cellCount());
}

int HeightMap::width() const {
    return m_width;
}

int HeightMap::height() const {
    return m_height;
}

int HeightMap::cellCount() const {
    return m_width * m_height;
}

double HeightMap::cellAreaSquareMeters() const {
    return m_cellAreaSquareMeters;
}

void HeightMap::resize(int width, int height) {
    m_width = width;
    m_height = height;
    m_heightsMeters.resize(cellCount());
}

void HeightMap::setCellAreaSquareMeters(double areaSquareMeters) {
    m_cellAreaSquareMeters = areaSquareMeters;
}

void HeightMap::setHeights(const QVector<double> &heightsMeters) {
    m_heightsMeters = heightsMeters;
}

const QVector<double> &HeightMap::heightsMeters() const {
    return m_heightsMeters;
}

HeightMap::FloodingResult HeightMap::computeFlooding(double liquidVolumeCubicMeters) const {
    FloodingResult result;
    const int totalCells = std::min(cellCount(), m_heightsMeters.size());
    if (totalCells <= 0 || liquidVolumeCubicMeters <= 0.0) {
        return result;
    }

    QVector<double> sortedHeights = m_heightsMeters.mid(0, totalCells);
    std::sort(sortedHeights.begin(), sortedHeights.end());

    const double cellArea = m_cellAreaSquareMeters;
    double remainingVolume = liquidVolumeCubicMeters;
    double currentLevel = sortedHeights.front();

    for (int i = 0; i < totalCells; ++i) {
        const double nextLevel = (i + 1 < totalCells) ? sortedHeights[i + 1] : sortedHeights[i];
        const double levelDelta = nextLevel - currentLevel;
        if (levelDelta > 0.0) {
            // Volume needed to raise the filled surface by levelDelta over (i+1) lowest cells:
            // V = delta_h * area_cell * N.
            const double requiredVolume = levelDelta * cellArea * (i + 1);
            if (remainingVolume < requiredVolume) {
                const double rise = remainingVolume / (cellArea * (i + 1));
                result.seaLevelMeters = currentLevel + rise;
                result.floodedAreaSquareMeters = cellArea * (i + 1);
                return result;
            }
            remainingVolume -= requiredVolume;
            currentLevel = nextLevel;
        }
    }

    // If we have more volume than needed to cover the highest cell, extend the sea level above it.
    result.seaLevelMeters = currentLevel + remainingVolume / (cellArea * totalCells);
    result.floodedAreaSquareMeters = cellArea * totalCells;
    return result;
}
