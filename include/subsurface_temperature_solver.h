#pragma once

#include "subsurface_grid.h"
#include "thermal_conductivity_model.h"

#include <QtCore/QVector>

enum class SubsurfaceBottomBoundaryCondition {
    Insulating,
    FixedTemperature
};

struct SubsurfaceModelSettings {
    int layerCount = 24;
    // Верхний слой ~5 см: ближе к типичной суточной глубине прогрева лунного реголита.
    double topLayerThicknessMeters = 0.05;
    double bottomDepthMeters = 2.0;
    SubsurfaceBottomBoundaryCondition bottomBoundary =
        SubsurfaceBottomBoundaryCondition::Insulating;
    double geothermalFluxWPerM2 = 0.0;
};

class SubsurfaceTemperatureSolver {
public:
    SubsurfaceTemperatureSolver() = default;
    SubsurfaceTemperatureSolver(const SubsurfaceGrid &grid,
                                double thermalConductivity,
                                const ThermalConductivityModel &thermalConductivityModel,
                                double density,
                                double specificHeat,
                                SubsurfaceBottomBoundaryCondition bottomBoundary,
                                double bottomTemperatureKelvin,
                                double geothermalFluxWPerM2);

    void reset(const SubsurfaceGrid &grid,
               double thermalConductivity,
               const ThermalConductivityModel &thermalConductivityModel,
               double density,
               double specificHeat,
               SubsurfaceBottomBoundaryCondition bottomBoundary,
               double bottomTemperatureKelvin,
               double geothermalFluxWPerM2);

    void setInitialTemperature(double temperatureKelvin);
    void setTemperatures(const QVector<double> &temperatures);

    const QVector<double> &temperatures() const;
    double surfaceTemperatureKelvin() const;
    double topLayerHeatCapacity() const;

    void stepImplicit(double netSurfaceFlux, double dtSeconds);

    void setBottomBoundary(SubsurfaceBottomBoundaryCondition bottomBoundary,
                           double bottomTemperatureKelvin);

private:
    double conductivityAtTemperature(double temperatureKelvin) const;
    void solveTridiagonal(QVector<double> &a,
                          QVector<double> &b,
                          QVector<double> &c,
                          QVector<double> &d) const;

    SubsurfaceGrid grid_;
    QVector<double> temperatures_;
    double thermalConductivity_ = 1.0;
    ThermalConductivityModel thermalConductivityModel_{};
    double density_ = 1.0;
    double specificHeat_ = 1.0;
    double alpha_ = 1.0;
    SubsurfaceBottomBoundaryCondition bottomBoundary_ =
        SubsurfaceBottomBoundaryCondition::Insulating;
    double bottomTemperatureKelvin_ = 3.0;
    double geothermalFluxWPerM2_ = 0.0;
};
