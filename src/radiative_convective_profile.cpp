#include "radiative_convective_profile.h"

#include "atmospheric_thermodynamics.h"

#include <QtCore/QtMath>

#include <algorithm>
#include <cmath>

namespace {
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kPascalPerAtm = 101325.0;
constexpr double kReferencePressurePa = 101325.0;
constexpr double kReferenceTemperatureKelvin = 288.0;
constexpr double kOpacityPressureExponent = 0.5;
constexpr double kOpacityTemperatureExponent = 0.3;
constexpr double kMinimumOpacity = 1e-5;
constexpr double kTopPressureFloorAtm = 1e-4;
constexpr double kMinPressureRatio = 1e-4;
}  // namespace

RadiativeConvectiveProfile::RadiativeConvectiveProfile(const AtmosphereComposition &composition,
                                                       double surfacePressureAtm,
                                                       double surfaceTemperatureKelvin,
                                                       double effectiveTemperatureKelvin,
                                                       double surfaceGravity,
                                                       int layerCount)
    : composition_(composition),
      surfacePressureAtm_(surfacePressureAtm),
      surfaceTemperatureKelvin_(surfaceTemperatureKelvin),
      effectiveTemperatureKelvin_(effectiveTemperatureKelvin),
      surfaceGravity_(surfaceGravity),
      layerCount_(layerCount) {
    buildProfile();
}

const QVector<RadiativeConvectiveLayer> &RadiativeConvectiveProfile::layers() const {
    return layers_;
}

double RadiativeConvectiveProfile::totalOpticalDepthLongwave() const {
    return totalOpticalDepthLongwave_;
}

double RadiativeConvectiveProfile::totalOpticalDepthShortwave() const {
    return totalOpticalDepthShortwave_;
}

double RadiativeConvectiveProfile::mixedOpacity(double temperatureKelvin,
                                                double pressurePa,
                                                bool shortwave) const {
    double kappaMix = 0.0;
    double totalShare = 0.0;
    const auto gases = availableGases();
    const auto fractions = composition_.fractions();
    for (const auto &fraction : fractions) {
        if (fraction.share <= 0.0) {
            continue;
        }
        const auto it = std::find_if(gases.begin(), gases.end(),
                                     [&fraction](const GasSpec &spec) {
                                         return spec.id == fraction.id;
                                     });
        if (it == gases.end()) {
            continue;
        }
        const double kappaBase =
            shortwave
                ? (it->isGreenhouse ? 0.03 : 0.005)
                : (it->isGreenhouse ? 0.20 : 0.01);
        kappaMix += fraction.share * kappaBase;
        totalShare += fraction.share;
    }

    if (totalShare <= 0.0) {
        kappaMix = shortwave ? 0.005 : 0.01;
    } else {
        kappaMix /= totalShare;
    }

    // Учитываем зависимость κ от давления и температуры:
    // κ(T, p) = κ0 * (p / p0)^a * (T / T0)^b, где a ~ 0.5, b ~ 0.3.
    const double pressureScale =
        std::pow(qMax(1e-6, pressurePa / kReferencePressurePa),
                 kOpacityPressureExponent);
    const double temperatureScale =
        std::pow(qMax(1e-6, temperatureKelvin / kReferenceTemperatureKelvin),
                 kOpacityTemperatureExponent);
    return qMax(kMinimumOpacity, kappaMix * pressureScale * temperatureScale);
}

double RadiativeConvectiveProfile::kappaLongwave(double temperatureKelvin,
                                                 double pressurePa) const {
    return mixedOpacity(temperatureKelvin, pressurePa, false);
}

double RadiativeConvectiveProfile::kappaShortwave(double temperatureKelvin,
                                                  double pressurePa) const {
    return mixedOpacity(temperatureKelvin, pressurePa, true);
}

double RadiativeConvectiveProfile::radiativeGradientDlnT(double temperatureKelvin,
                                                         double pressurePa) const {
    if (temperatureKelvin <= 0.0 || pressurePa <= 0.0 || surfaceGravity_ <= 0.0) {
        return 0.0;
    }

    // Радиативный градиент в терминах d ln T / d ln P для серой атмосферы:
    // ∇_rad = (3 κ P F) / (16 g σ T^4),
    // где F = σ T_eff^4 — вертикальный поток из ТОА-баланса.
    const double kappa = kappaLongwave(temperatureKelvin, pressurePa);
    // Если эффективная температура не задана, падаем обратно на температуру поверхности,
    // чтобы профиль оставался устойчивым в упрощённых режимах.
    const double effectiveTemperature =
        (effectiveTemperatureKelvin_ > 0.0) ? effectiveTemperatureKelvin_
                                            : surfaceTemperatureKelvin_;
    const double flux =
        kStefanBoltzmannConstant * std::pow(qMax(1.0, effectiveTemperature), 4.0);
    const double numerator = 3.0 * kappa * pressurePa * flux;
    const double denominator =
        16.0 * qMax(1.0, surfaceGravity_) * kStefanBoltzmannConstant *
        std::pow(temperatureKelvin, 4.0);
    return numerator / qMax(1e-6, denominator);
}

double RadiativeConvectiveProfile::adiabaticGradientDlnT(double temperatureKelvin,
                                                         double pressurePa) const {
    const double pressureAtm = pressurePa / kPascalPerAtm;
    return AtmosphericThermodynamics::adiabaticGradientDlnT(composition_,
                                                            pressureAtm,
                                                            temperatureKelvin,
                                                            surfaceGravity_);
}

void RadiativeConvectiveProfile::buildProfile() {
    layers_.clear();
    totalOpticalDepthLongwave_ = 0.0;
    totalOpticalDepthShortwave_ = 0.0;

    if (surfacePressureAtm_ <= 0.0 || surfaceTemperatureKelvin_ <= 0.0 ||
        surfaceGravity_ <= 0.0 || layerCount_ <= 0) {
        return;
    }

    const double topPressureAtm =
        qMax(kTopPressureFloorAtm, surfacePressureAtm_ * kMinPressureRatio);
    const double surfacePressurePa = surfacePressureAtm_ * kPascalPerAtm;
    const double topPressurePa = topPressureAtm * kPascalPerAtm;

    const double logSurface = std::log(surfacePressurePa);
    const double logTop = std::log(topPressurePa);

    QVector<double> edges;
    edges.reserve(layerCount_ + 1);
    for (int i = 0; i <= layerCount_; ++i) {
        const double frac = static_cast<double>(i) / layerCount_;
        const double logP = logSurface + (logTop - logSurface) * frac;
        edges.push_back(std::exp(logP));
    }

    layers_.reserve(layerCount_);
    double temperatureLower = surfaceTemperatureKelvin_;
    const double rSpecific = AtmosphericThermodynamics::specificGasConstant(composition_);

    for (int i = 0; i < layerCount_; ++i) {
        const double pLower = edges[i];
        const double pUpper = edges[i + 1];
        const double pMid = std::sqrt(pLower * pUpper);
        const double dlnP = std::log(pUpper / pLower);

        // Конвективно-радиативный профиль: выбираем меньший из радиативного и адиабатического
        // градиентов. Если ∇_rad > ∇_ad, то слой становится конвективным (сухая/влажная адиабата).
        const double gradRad = radiativeGradientDlnT(temperatureLower, pLower);
        const double gradAd = adiabaticGradientDlnT(temperatureLower, pLower);
        const double grad = qMax(0.0, qMin(gradRad, gradAd));

        const double temperatureUpper =
            temperatureLower * std::exp(grad * dlnP);
        const double temperatureMid = std::sqrt(temperatureLower * temperatureUpper);

        const double density = pMid / qMax(1e-6, rSpecific * temperatureMid);
        const double deltaP = pLower - pUpper;
        // Гидростатика: dP/dz = -ρ g ⇒ dz = -dP / (ρ g).
        const double thickness = deltaP / qMax(1e-6, density * surfaceGravity_);

        const double kappaLw = kappaLongwave(temperatureMid, pMid);
        const double kappaSw = kappaShortwave(temperatureMid, pMid);
        const double tauLw = kappaLw * density * thickness;
        const double tauSw = kappaSw * density * thickness;

        RadiativeConvectiveLayer layer;
        layer.pressureAtm = pMid / kPascalPerAtm;
        layer.temperatureKelvin = temperatureMid;
        layer.densityKgPerM3 = density;
        layer.thicknessMeters = thickness;
        layer.opticalDepthLongwave = tauLw;
        layer.opticalDepthShortwave = tauSw;
        layers_.push_back(layer);

        totalOpticalDepthLongwave_ += tauLw;
        totalOpticalDepthShortwave_ += tauSw;
        temperatureLower = temperatureUpper;
    }
}
