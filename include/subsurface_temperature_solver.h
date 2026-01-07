#pragma once

#include "subsurface_grid.h"

#include <QtCore/QMetaType>
#include <QtCore/QString>
#include <QtCore/QStringList>
#include <QtCore/QVector>

#include <algorithm>
#include <cmath>

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
    // Профиль теплопроводности (λ), плотности (ρ) и теплоемкости (c).
    // Если задан один элемент — он копируется на все слои; если два — интерполируются
    // от поверхности к дну; при большем числе значений интерполяция линейная
    // по нормированной глубине.
    QVector<double> thermalConductivityByLayer;
    QVector<double> densityByLayer;
    QVector<double> specificHeatByLayer;

    QVector<double> thermalConductivityProfile(int layers, double defaultValue) const;
    QVector<double> densityProfile(int layers, double defaultValue) const;
    QVector<double> specificHeatProfile(int layers, double defaultValue) const;
    QString profileSignature() const;

private:
    static QVector<double> expandProfile(const QVector<double> &values,
                                         int layers,
                                         double defaultValue);
};

Q_DECLARE_METATYPE(SubsurfaceModelSettings)

class SubsurfaceTemperatureSolver {
public:
    SubsurfaceTemperatureSolver() = default;
    SubsurfaceTemperatureSolver(const SubsurfaceGrid &grid,
                                const QVector<double> &thermalConductivityByLayer,
                                const QVector<double> &densityByLayer,
                                const QVector<double> &specificHeatByLayer,
                                SubsurfaceBottomBoundaryCondition bottomBoundary,
                                double bottomTemperatureKelvin);

    void reset(const SubsurfaceGrid &grid,
               const QVector<double> &thermalConductivityByLayer,
               const QVector<double> &densityByLayer,
               const QVector<double> &specificHeatByLayer,
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
    QVector<double> thermalConductivityByLayer_;
    QVector<double> densityByLayer_;
    QVector<double> specificHeatByLayer_;
    SubsurfaceBottomBoundaryCondition bottomBoundary_ =
        SubsurfaceBottomBoundaryCondition::Insulating;
    double bottomTemperatureKelvin_ = 3.0;
};

inline QVector<double> SubsurfaceModelSettings::thermalConductivityProfile(
    int layers,
    double defaultValue) const {
    return expandProfile(thermalConductivityByLayer, layers, defaultValue);
}

inline QVector<double> SubsurfaceModelSettings::densityProfile(int layers,
                                                               double defaultValue) const {
    return expandProfile(densityByLayer, layers, defaultValue);
}

inline QVector<double> SubsurfaceModelSettings::specificHeatProfile(
    int layers,
    double defaultValue) const {
    return expandProfile(specificHeatByLayer, layers, defaultValue);
}

inline QString SubsurfaceModelSettings::profileSignature() const {
    const auto serialize = [](const QVector<double> &values) -> QString {
        if (values.isEmpty()) {
            return QStringLiteral("none");
        }
        QStringList parts;
        parts.reserve(values.size());
        for (double value : values) {
            parts << QString::number(value, 'g', 10);
        }
        return parts.join(',');
    };
    return QStringLiteral("k=%1|rho=%2|cp=%3")
        .arg(serialize(thermalConductivityByLayer),
             serialize(densityByLayer),
             serialize(specificHeatByLayer));
}

inline QVector<double> SubsurfaceModelSettings::expandProfile(const QVector<double> &values,
                                                              int layers,
                                                              double defaultValue) {
    QVector<double> result;
    if (layers <= 0) {
        return result;
    }
    if (values.isEmpty()) {
        result.fill(defaultValue, layers);
        return result;
    }
    if (values.size() == 1) {
        result.fill(values.first(), layers);
        return result;
    }
    if (values.size() == layers) {
        return values;
    }
    const int sourceCount = values.size();
    result.resize(layers);
    if (layers == 1) {
        result[0] = values.first();
        return result;
    }
    for (int i = 0; i < layers; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(layers - 1);
        const double srcPosition = t * static_cast<double>(sourceCount - 1);
        const int leftIndex = std::max(0, std::min(sourceCount - 1,
                                                   static_cast<int>(std::floor(srcPosition))));
        const int rightIndex = std::min(sourceCount - 1, leftIndex + 1);
        const double mix = srcPosition - static_cast<double>(leftIndex);
        result[i] = values[leftIndex] * (1.0 - mix) + values[rightIndex] * mix;
    }
    return result;
}
