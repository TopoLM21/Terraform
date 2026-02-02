#pragma once

#include <QVector>

class HeightMap {
public:
    struct FloodingResult {
        double seaLevelMeters = 0.0;
        double floodedAreaSquareMeters = 0.0;
    };

    HeightMap();
    HeightMap(int width, int height, double cellAreaSquareMeters);

    int width() const;
    int height() const;
    int cellCount() const;
    double cellAreaSquareMeters() const;

    void resize(int width, int height);
    void setCellAreaSquareMeters(double areaSquareMeters);

    void setHeights(const QVector<double> &heightsMeters);
    const QVector<double> &heightsMeters() const;

    // Assumptions:
    // - The height map is discretized into rectangular cells of equal area.
    // - Each cell has a single representative height (flat topography inside the cell).
    // - Liquid is incompressible and fills from the lowest cells upward until the
    //   specified volume is exhausted.
    FloodingResult computeFlooding(double liquidVolumeCubicMeters) const;

private:
    int m_width = 0;
    int m_height = 0;
    double m_cellAreaSquareMeters = 1.0;
    QVector<double> m_heightsMeters;
};
