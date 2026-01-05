#include "spectral_band_model.h"

#include <QtCore/QtMath>

#include <cmath>

namespace {
constexpr double kMinimumOpacity = 1e-6;

double bandOpacity(const SpectralBandParameters &band,
                   double temperatureKelvin,
                   double pressurePa) {
    // Формула масштабирования полосы:
    // κ_band(T, p) = κ0 * (p / p0)^a * (T / T0)^b,
    // где κ0 — табличная непрозрачность в полосе (HITRAN/MT_CKD/калибровка).
    const double pressureScale =
        std::pow(qMax(1e-12, pressurePa / band.referencePressurePa),
                 band.pressureExponent);
    const double temperatureScale =
        std::pow(qMax(1e-12, temperatureKelvin / band.referenceTemperatureKelvin),
                 band.temperatureExponent);
    return qMax(kMinimumOpacity, band.referenceOpacity * pressureScale * temperatureScale);
}
}  // namespace

SpectralBandModel SpectralBandModel::defaultModel() {
    SpectralBandModel model;

    // Параметры полос задаём так, чтобы грубо воспроизводить ИК-поглощение CO2 и окно.
    // Значения κ0 подбирались калибровкой для Венеры (P≈90 атм, 96.5% CO2 → T≈737 K),
    // с опорой на HITRAN/MT_CKD для ориентиров по центрам полос.
    const QVector<SpectralBandParameters> co2Bands = {
        // 4.3 мкм (асимметричный режим), сильная полоса.
        {QStringLiteral("CO2 4.3um"), 4.3, 0.7, 0.25, 101325.0, 288.0, 0.7, -0.2, false,
         QStringLiteral("HITRAN ориентир, κ0 калибровано")},
        // 15 мкм (основная тепловая полоса).
        {QStringLiteral("CO2 15um"), 15.0, 3.0, 0.40, 101325.0, 288.0, 0.8, -0.25, false,
         QStringLiteral("HITRAN ориентир, κ0 калибровано")},
        // Атмосферное окно 8-12 мкм (эффективное слабое поглощение).
        {QStringLiteral("CO2 window"), 10.0, 4.0, 0.05, 101325.0, 288.0, 0.6, -0.1, false,
         QStringLiteral("MT_CKD ориентир окна, κ0 калибровано")},
        // Видимый/ближний ИК вклад в коротковолновый диапазон.
        {QStringLiteral("CO2 SW"), 1.6, 0.5, 0.01, 101325.0, 288.0, 0.5, 0.0, true,
         QStringLiteral("HITRAN ориентир, κ0 калибровано")},
    };
    model.setBands(QStringLiteral("co2"), co2Bands);

    const QVector<SpectralBandParameters> h2oBands = {
        {QStringLiteral("H2O 6.3um"), 6.3, 2.0, 0.20, 101325.0, 288.0, 0.8, -0.3, false,
         QStringLiteral("HITRAN ориентир, κ0 по умолчанию")},
        {QStringLiteral("H2O window"), 12.0, 4.0, 0.04, 101325.0, 288.0, 0.6, -0.15, false,
         QStringLiteral("MT_CKD ориентир окна, κ0 по умолчанию")},
        {QStringLiteral("H2O SW"), 1.4, 0.4, 0.02, 101325.0, 288.0, 0.5, 0.0, true,
         QStringLiteral("HITRAN ориентир, κ0 по умолчанию")},
    };
    model.setBands(QStringLiteral("h2o"), h2oBands);

    const QVector<SpectralBandParameters> so2Bands = {
        {QStringLiteral("SO2 7.3um"), 7.3, 1.2, 0.10, 101325.0, 288.0, 0.7, -0.2, false,
         QStringLiteral("HITRAN ориентир, κ0 по умолчанию")},
        {QStringLiteral("SO2 SW"), 0.3, 0.1, 0.05, 101325.0, 288.0, 0.4, 0.0, true,
         QStringLiteral("HITRAN ориентир, κ0 по умолчанию")},
    };
    model.setBands(QStringLiteral("so2"), so2Bands);

    return model;
}

void SpectralBandModel::setBands(const QString &gasId,
                                 const QVector<SpectralBandParameters> &bands) {
    gasBands_[gasId.toLower()] = bands;
}

double SpectralBandModel::kappa(const AtmosphereComposition &composition,
                                double temperatureKelvin,
                                double pressurePa,
                                bool shortwave) const {
    const auto fractions = composition.fractions();
    if (fractions.isEmpty()) {
        return kMinimumOpacity;
    }

    double kappaSum = 0.0;
    double totalShare = 0.0;

    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }
        const auto it = gasBands_.constFind(fraction.id.toLower());
        if (it == gasBands_.constEnd()) {
            continue;
        }
        for (const auto &band : it.value()) {
            if (band.shortwave != shortwave) {
                continue;
            }
            // Суммарная непрозрачность смеси = сумма по газам доли * сумма полос.
            kappaSum += fraction.share * bandOpacity(band, temperatureKelvin, pressurePa);
        }
        totalShare += fraction.share;
    }

    if (totalShare <= 0.0) {
        return kMinimumOpacity;
    }

    // Нормируем по суммарной доле, чтобы смесь сохраняла ожидаемый масштаб.
    return qMax(kMinimumOpacity, kappaSum / totalShare);
}
