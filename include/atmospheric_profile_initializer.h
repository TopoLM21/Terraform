#pragma once

#include "atmosphere_model.h"

#include <QVector>

struct AtmosphericProfileLayer {
    double heightMeters = 0.0;
    double thicknessMeters = 0.0;
    double temperatureKelvin = 0.0;
    double pressureAtm = 0.0;
    double densityKgPerM3 = 0.0;
};

// Утилита для построения стартового вертикального профиля атмосферы.
class AtmosphericProfileInitializer {
public:
    struct Settings {
        double surfaceTemperatureKelvin = 0.0;
        double surfacePressureAtm = 0.0;
        double gravityMps2 = 0.0;
        int layerCount = 0;
        // Если <= 0, будет выбран сухой адиабатический градиент (при useDryAdiabatic).
        double lapseRateKPerM = 0.0;
        // Разрешает использовать сухой адиабатический градиент при отсутствии явного lapse rate.
        bool useDryAdiabatic = true;
        // Можно указать вручную; если 0, высота берётся как ~6 масштабных высот.
        double topHeightMeters = 0.0;
        // Верхняя граница по давлению (атм): если > 0, профиль строится до P <= порога.
        double minTopPressureAtm = 0.0;
        // Минимальная толщина нижнего слоя (м). Используется для уплотнения сетки у поверхности.
        double minBottomLayerThicknessMeters = 0.0;
    };

    explicit AtmosphericProfileInitializer(const AtmosphereComposition &composition);

    QVector<AtmosphericProfileLayer> buildProfile(const Settings &settings) const;

private:
    double resolveLapseRateKPerM(const Settings &settings) const;
    double scaleHeightMeters(double temperatureKelvin, double gravityMps2) const;
    double pressureAtmAtHeight(double surfacePressureAtm,
                               double heightMeters,
                               double scaleHeightMeters) const;

    AtmosphereComposition composition_;
};
