#pragma once

#include <QtCore/QVector>

class SubsurfaceGrid {
public:
    SubsurfaceGrid() = default;
    SubsurfaceGrid(int layerCount, double topLayerThicknessMeters, double bottomDepthMeters);

    void rebuild(int layerCount, double topLayerThicknessMeters, double bottomDepthMeters);

    int layerCount() const;
    double bottomDepthMeters() const;
    const QVector<double> &layerThicknessesMeters() const;
    const QVector<double> &layerDepthsMeters() const;

private:
    QVector<double> layerThicknessesMeters_;
    QVector<double> layerDepthsMeters_;
    double bottomDepthMeters_ = 0.0;
};
