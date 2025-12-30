#pragma once

#include "subsurface_grid.h"

#include <QtCore/QVector>

enum class SubsurfaceBottomBoundaryCondition {
    Insulating,
    FixedTemperature
};

struct SubsurfaceModelSettings {
    int layerCount = 24;
    double topLayerThicknessMeters = 0.02;
    double bottomDepthMeters = 2.0;
    SubsurfaceBottomBoundaryCondition bottomBoundary =
        SubsurfaceBottomBoundaryCondition::Insulating;
};

class SubsurfaceTemperatureSolver {
public:
    SubsurfaceTemperatureSolver() = default;
    SubsurfaceTemperatureSolver(const SubsurfaceGrid &grid,
                                double thermalConductivity,
                                double density,
                                double specificHeat,
                                SubsurfaceBottomBoundaryCondition bottomBoundary,
                                double bottomTemperatureKelvin);

    void reset(const SubsurfaceGrid &grid,
               double thermalConductivity,
               double density,
               double specificHeat,
               SubsurfaceBottomBoundaryCondition bottomBoundary,
               double bottomTemperatureKelvin);

    void setInitialTemperature(double temperatureKelvin);
    void setTemperatures(const QVector<double> &temperatures);

    const QVector<double> &temperatures() const;
    double surfaceTemperatureKelvin() const;

    void stepImplicit(double netSurfaceFlux, double dtSeconds);

    void setBottomBoundary(SubsurfaceBottomBoundaryCondition bottomBoundary,
                           double bottomTemperatureKelvin);

private:
    void solveTridiagonal(QVector<double> &a,
                          QVector<double> &b,
                          QVector<double> &c,
                          QVector<double> &d) const;

    SubsurfaceGrid grid_;
    QVector<double> temperatures_;
    double thermalConductivity_ = 1.0;
    double density_ = 1.0;
    double specificHeat_ = 1.0;
    double alpha_ = 1.0;
    SubsurfaceBottomBoundaryCondition bottomBoundary_ =
        SubsurfaceBottomBoundaryCondition::Insulating;
    double bottomTemperatureKelvin_ = 3.0;
};
