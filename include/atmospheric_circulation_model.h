#pragma once

#include "rotation_mode.h"
#include "surface_temperature_calculator.h"

#include <QtCore/QVector>

class AtmosphericCirculationModel {
public:
    AtmosphericCirculationModel(double dayLengthDays,
                                RotationMode rotationMode,
                                double atmospherePressureAtm);

    void setAtmosphereMassKg(double atmosphereMassKg);
    void setMeridionalTransportSteps(int steps);

    QVector<TemperatureRangePoint> applyHeatTransport(
        const QVector<TemperatureRangePoint> &points) const;

private:
    double dayLengthDays_ = 1.0;
    RotationMode rotationMode_ = RotationMode::Normal;
    double atmospherePressureAtm_ = 0.0;
    double atmosphereMassKg_ = 0.0;
    int meridionalTransportSteps_ = 8;
    int cellsPerHemisphere_ = 3;
    double cellSizeDegrees_ = 30.0;

    double baseMeridionalMixingCoefficient() const;
    double meridionalMixingCoefficient(double latitudeDegrees) const;
    double dayNightExchangeCoefficient() const;
    int cellIndexForAxis(double axisDegrees) const;
    double cellCouplingFactor(int leftCell, int rightCell) const;
};
