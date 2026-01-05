#pragma once

#include "atmosphere_model.h"

#include <QtCore/QHash>
#include <QtCore/QString>
#include <QtCore/QVector>

struct SpectralBandParameters {
    QString label;
    double centerMicron = 0.0;
    double widthMicron = 0.0;
    // Массовая непрозрачность κ0 в м^2/кг для полосы при (p0, T0).
    double referenceOpacity = 0.0;
    double referencePressurePa = 101325.0;
    double referenceTemperatureKelvin = 288.0;
    double pressureExponent = 0.0;
    double temperatureExponent = 0.0;
    bool shortwave = false;
    QString source;
};

class SpectralBandModel {
public:
    SpectralBandModel() = default;

    static SpectralBandModel defaultModel();
    void setBands(const QString &gasId, const QVector<SpectralBandParameters> &bands);

    // Суммарная непрозрачность κ(T, p) по спектральным полосам и газам
    // в зависимости от состава атмосферы.
    double kappa(const AtmosphereComposition &composition,
                 double temperatureKelvin,
                 double pressurePa,
                 bool shortwave) const;

private:
    QHash<QString, QVector<SpectralBandParameters>> gasBands_;
};
