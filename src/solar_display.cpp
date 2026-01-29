#include "mode_illustration_widget.h"
#include "orbit_animation_model.h"
#include "planet_presets.h"
#include "segment_selector_widget.h"
#include "solar_calculator.h"
#include "solar_display.h"
#include "atmosphere_widget.h"
#include "heightmap_storage.h"
#include "surface_tile_temperature_calculator.h"
#include "surface_temperature_plot.h"
#include "surface_map_widget.h"
#include "surface_globe_widget.h"
#include "surface_point_status_dialog.h"
#include "surface_atmosphere_coupler.h"
#include "surface_temperature_scale_widget.h"
#include "surface_height_scale_widget.h"
#include "surface_heightmap.h"
#include "surface_wind_scale_widget.h"
#include "surface_pressure_scale_widget.h"
#include "surface_realistic_scale_widget.h"
#include "surface_map_mode.h"
#include "StarSystemTopView.h"
#include "surface_advection_model.h"
#include "surface_pressure_transport_model.h"
#include "atmosphere_step_solver.h"
#include "wind_field_model.h"
#include "rotation_period_utils.h"
#include "planet_surface_grid.h"
#include "subsurface_temperature_solver.h"
#include "orbit_segment_calculator.h"
#include "radiation_model.h"
#include "layered_radiation_model.h"
#include "radiation_model_utils.h"
#include "atmospheric_pressure_model.h"
#include "atmospheric_profile_initializer.h"
#include "atmospheric_cell_state.h"
#include "atmosphere_model.h"

#include <QtCore/QCommandLineOption>
#include <QtCore/QCommandLineParser>
#include <QtCore/QElapsedTimer>
#include <QtCore/QLocale>
#include <QtCore/QPointer>
#include <QtCore/QSignalBlocker>
#include <QtCore/QThread>
#include <QtCore/QThreadPool>
#include <QtCore/QTimer>
// #include <QtCore/QOverload>
#include <QtGlobal>
#include <QtCore/QtMath>
#include <QtCore/QSet>
#include <QtCore/QFutureWatcher>
#include <QtCore/QHash>
#include <QtCore/QStringList>
#include <QtCore/QLoggingCategory>
#include <QtGui/QDoubleValidator>
#include <QtGui/QImage>
#include <QtGui/QVector3D>
#include <QtConcurrent/QtConcurrent>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QProgressDialog>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QToolButton>
#include <algorithm>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QStyle>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <QtCore/QPointF>

#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>

namespace {
Q_LOGGING_CATEGORY(solarRadiationLog, "solar.radiation")

QString solarLoggingBaseRules() {
    static const QString baseRules = QString::fromLocal8Bit(qgetenv("QT_LOGGING_RULES"));
    return baseRules;
}

void setSolarDebugLoggingEnabled(bool enabled) {
    QString rules = solarLoggingBaseRules();
    if (!rules.isEmpty()) {
        rules.append('\n');
    }
    const QString flagValue = enabled ? QStringLiteral("true") : QStringLiteral("false");
    rules.append(QStringLiteral("solar.radiation.info=%1\n"
                                "solar.radiation.debug=%1\n"
                                "solar.atmosphere.profile.info=%1\n"
                                "solar.atmosphere.profile.debug=%1\n"
                                "solar.atmosphere.vertical_mixing.info=%1\n"
                                "solar.atmosphere.vertical_mixing.debug=%1")
                     .arg(flagValue));
    QLoggingCategory::setFilterRules(rules);
}

bool isSolarRadiationLoggingEnabledFromEnvironment() {
    if (!qEnvironmentVariableIsSet("SOLAR_RADIATION_LOG")) {
        return false;
    }

    const QByteArray rawValue = qgetenv("SOLAR_RADIATION_LOG").trimmed().toLower();
    if (rawValue.isEmpty()) {
        return true;
    }
    return rawValue != "0" && rawValue != "false" && rawValue != "off";
}

void enableSolarRadiationLogging() {
    setSolarDebugLoggingEnabled(true);
}

double denseAtmosphereMinTemperatureKelvin(double effectiveTemperatureKelvin,
                                           double pressureAtm,
                                           double greenhouseOpacity,
                                           double co2Share,
                                           double minDenseAtmosphereTemperatureK) {
    if (effectiveTemperatureKelvin <= 0.0) {
        return 0.0;
    }
    // t_eff отражает баланс излучения на ВГБ (TOA), но в плотных атмосферах
    // излучающий слой поднимается вверх, а нижние слои дополнительно греются
    // из-за оптической толщины, давления и уширения линий (особенно CO₂).
    // Поэтому для базовой температуры берём повышающую поправку по давлению,
    // составу и (для сверхплотных атмосфер) отдельный коэффициент.
    const double safePressureAtm = qMax(0.0, pressureAtm);
    const double co2Fraction = qBound(0.0, co2Share, 1.0);
    const double pressureBoost = 1.0 + 0.035 * std::log1p(safePressureAtm);
    const double greenhouseBoost = 1.0 + 0.2 * qBound(0.0, greenhouseOpacity, 1.0);
    const double compositionBoost = 1.0 + 0.8 * co2Fraction;
    const double superDenseBoost =
        (safePressureAtm > 10.0)
            ? (1.0 + 0.25 * std::log1p(safePressureAtm / 10.0))
            : 1.0;
    double minTemperature =
        effectiveTemperatureKelvin * pressureBoost * greenhouseBoost * compositionBoost * superDenseBoost;
    if (minDenseAtmosphereTemperatureK > 0.0) {
        minTemperature = qMax(minTemperature, minDenseAtmosphereTemperatureK);
    }
    return minTemperature;
}

constexpr int kRoleSemiMajorAxis = Qt::UserRole;
constexpr int kRoleIsCustom = Qt::UserRole + 1;
constexpr int kRolePlanetName = Qt::UserRole + 2;
constexpr int kRoleMaterialId = Qt::UserRole + 3;
constexpr int kRoleDayLength = Qt::UserRole + 4;
constexpr int kRoleYearLengthDays = Qt::UserRole + 5;
constexpr int kRoleEccentricity = Qt::UserRole + 6;
constexpr int kRoleObliquity = Qt::UserRole + 7;
constexpr int kRolePerihelionArgument = Qt::UserRole + 8;
constexpr int kRoleMassEarths = Qt::UserRole + 9;
constexpr int kRoleRadiusKm = Qt::UserRole + 10;
constexpr int kRoleRotationMode = Qt::UserRole + 11;
constexpr int kRoleAtmosphere = Qt::UserRole + 12;
constexpr int kRoleGreenhouseOpacity = Qt::UserRole + 13;
constexpr int kRoleCloudAlbedo = Qt::UserRole + 14;
constexpr int kRoleHeightSourceType = Qt::UserRole + 15;
constexpr int kRoleHeightmapPath = Qt::UserRole + 16;
constexpr int kRoleHeightmapScaleKm = Qt::UserRole + 17;
constexpr int kRoleHeightSeed = Qt::UserRole + 18;
constexpr int kRoleUseContinentsHeight = Qt::UserRole + 19;
constexpr int kRoleHasSeaLevel = Qt::UserRole + 20;
constexpr int kRoleFlatHeight = Qt::UserRole + 21;
constexpr int kRoleManualGreenhouseOnTopOfAtmosphere = Qt::UserRole + 22;
constexpr int kRoleAdvancedRadiationModel = Qt::UserRole + 23;
constexpr int kRoleGeothermalFlux = Qt::UserRole + 24;
constexpr int kRoleDiurnalCoolingBiasK = Qt::UserRole + 25;
constexpr int kRolePrimaryStarRadius = Qt::UserRole + 26;
constexpr int kRolePrimaryStarTemperature = Qt::UserRole + 27;
constexpr int kRoleSecondaryStarRadius = Qt::UserRole + 28;
constexpr int kRoleSecondaryStarTemperature = Qt::UserRole + 29;
constexpr int kRoleHasSecondaryStar = Qt::UserRole + 30;
constexpr int kRoleStarPresetId = Qt::UserRole + 31;
constexpr int kRoleRotationModeOverride = Qt::UserRole + 32;
constexpr int kRoleSpinOrbitP = Qt::UserRole + 33;
constexpr int kRoleSpinOrbitQ = Qt::UserRole + 34;
constexpr int kRoleBinarySemiMajorAxis = Qt::UserRole + 35;
constexpr int kRoleBinaryPeriodDays = Qt::UserRole + 36;
constexpr int kRoleBinaryEccentricity = Qt::UserRole + 37;
constexpr int kRoleBinaryInclinationDegrees = Qt::UserRole + 38;
constexpr int kRoleBinaryArgumentPericenterA = Qt::UserRole + 39;
constexpr int kRoleBinaryArgumentPericenterB = Qt::UserRole + 40;
constexpr int kRoleHasBinaryOrbit = Qt::UserRole + 41;
constexpr int kRoleAtmosphereDisabled = Qt::UserRole + 42;
constexpr int kRoleAtmosphereBottomLayerThickness = Qt::UserRole + 43;
constexpr int kRoleMinDenseAtmosphereTemperature = Qt::UserRole + 44;
constexpr int kRoleVerticalWindMixingCoefficient = Qt::UserRole + 45;
constexpr int kRoleMinTopPressureAtm = Qt::UserRole + 46;
constexpr int kRoleDisableWaterAndClouds = Qt::UserRole + 47;
constexpr double kKelvinOffset = 273.15;
constexpr double kEarthRadiusKm = 6371.0;
constexpr double kEarthMassKg = 5.9722e24;
constexpr double kGravitationalConstant = 6.67430e-11;
constexpr int kSurfaceOrbitSegmentsPerYear = 360;
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;
constexpr double kStandardPressurePa = 101325.0;
constexpr double kDryAirSpecificHeatJPerKgK = 1004.0;
constexpr double kDefaultHeatTransferWPerM2K = 8.0;
constexpr double kEarthAreaKm2 = 510072000.0;
constexpr double kDefaultSurfaceRoughness = 20.0;
constexpr double kDefaultBasinShape = 3.5;
constexpr double kEarthWaterGigatons = 1.4e9;
constexpr int kSurfaceTemperatureHistoryDays = kSurfaceOrbitSegmentsPerYear;
constexpr double kDefaultAtmosphereBottomLayerThicknessMeters = 100.0;
constexpr double kDefaultVerticalWindMixingCoefficientKz = 1.0;
// 0.001 атм ≈ 1 гПа, используем как базовый порог верхней границы профиля.
constexpr double kDefaultMinTopPressureAtm = 0.001;
// Один реальный секундный тик соответствует часу системного времени (1/24 дня).
constexpr double kStarSystemDaysPerSecond = 1.0 / 24.0;
constexpr int kKeplerIterations = 8;

struct LongwaveFluxes {
    double downwardFluxWPerM2 = 0.0;
    double upwardFluxWPerM2 = 0.0;
};

LongwaveFluxes computeLongwaveFluxes(const RadiationModel &radiationModel,
                                     RadiationModelType radiationModelType,
                                     double surfaceEmittedFluxWPerM2,
                                     double airTemperatureKelvin,
                                     double greenhouseOpacity,
                                     bool useAtmosphericModel,
                                     bool manualGreenhouseOnTop) {
    LongwaveFluxes fluxes;
    if (radiationModelType == RadiationModelType::Layered) {
        const auto *layeredModel = dynamic_cast<const LayeredRadiationModel *>(&radiationModel);
        if (layeredModel) {
            const auto profile = layeredModel->longwaveFluxProfile(surfaceEmittedFluxWPerM2);
            if (!profile.downwardFluxWPerM2.isEmpty()) {
                fluxes.downwardFluxWPerM2 = profile.downwardFluxWPerM2.first();
            }
            if (!profile.upwardFluxWPerM2.isEmpty()) {
                fluxes.upwardFluxWPerM2 = profile.upwardFluxWPerM2.last();
            }
            return fluxes;
        }
    }

    const double baseTransmission = qBound(0.0, radiationModel.outgoingTransmission(), 1.0);
    const double greenhouseTransmission = qBound(0.0, 1.0 - greenhouseOpacity, 1.0);
    double transmission = baseTransmission;
    if (!useAtmosphericModel) {
        transmission = greenhouseTransmission;
    } else if (manualGreenhouseOnTop) {
        // Добавляем дополнительную оптическую толщину τ_gh к атмосферной τ_atm:
        // T_lw = exp(-(τ_atm + τ_gh)) ≈ T_atm * T_gh. Это серый двухпотоковый подход,
        // где излучение атмосферы задаётся ε(τ) = 1 - T_lw.
        const double extraTransmission = (baseTransmission > 0.0)
            ? qBound(0.0, greenhouseTransmission / baseTransmission, 1.0)
            : greenhouseTransmission;
        transmission = qBound(0.0, baseTransmission * extraTransmission, 1.0);
    }
    const double emissivity = 1.0 - transmission;
    // Серый двухпотоковый баланс: атмосфера излучает вверх и вниз как σ T_air^4 * ε(τ).
    const double atmosphericEmission =
        emissivity * kStefanBoltzmannConstant * std::pow(qMax(0.0, airTemperatureKelvin), 4.0);
    fluxes.downwardFluxWPerM2 = atmosphericEmission;
    fluxes.upwardFluxWPerM2 = atmosphericEmission + surfaceEmittedFluxWPerM2 * transmission;
    return fluxes;
}

double resolveSurfaceLongwaveUpFlux(const RadiationModel &radiationModel,
                                    RadiationModelType radiationModelType,
                                    double surfaceEmittedFluxWPerM2) {
    if (radiationModelType == RadiationModelType::Layered) {
        const auto *layeredModel = dynamic_cast<const LayeredRadiationModel *>(&radiationModel);
        if (layeredModel) {
            const auto profile = layeredModel->longwaveFluxProfile(surfaceEmittedFluxWPerM2);
            if (!profile.upwardFluxWPerM2.isEmpty()) {
                return profile.upwardFluxWPerM2.first();
            }
        }
    }
    // В однопотоковой модели поток вверх у поверхности равен собственному излучению.
    return surfaceEmittedFluxWPerM2;
}

enum class TemperaturePlotSource {
    SurfaceGrid = 0,
    LatitudeModel1D = 1,
};

double normalizeLongitudeRadians(double radians) {
    double wrapped = std::fmod(radians + M_PI, 2.0 * M_PI);
    if (wrapped < 0.0) {
        wrapped += 2.0 * M_PI;
    }
    return wrapped - M_PI;
}

double longitudeDistanceRadians(double a, double b) {
    return std::abs(normalizeLongitudeRadians(a - b));
}

double normalizeOrbitAngleRadians(double angle) {
    double wrapped = std::fmod(angle, 2.0 * M_PI);
    if (wrapped < 0.0) {
        wrapped += 2.0 * M_PI;
    }
    return wrapped;
}

double solveEccentricAnomalyRadians(double meanAnomaly, double eccentricity) {
    // Решаем уравнение Кеплера: M = E - e * sin(E).
    const double normalizedM = normalizeOrbitAngleRadians(meanAnomaly);
    double eccentricAnomaly = (eccentricity < 0.8) ? normalizedM : M_PI;
    for (int i = 0; i < kKeplerIterations; ++i) {
        const double f = eccentricAnomaly - eccentricity * qSin(eccentricAnomaly) - normalizedM;
        const double fPrime = 1.0 - eccentricity * qCos(eccentricAnomaly);
        if (qFuzzyIsNull(fPrime)) {
            break;
        }
        eccentricAnomaly -= f / fPrime;
    }
    return eccentricAnomaly;
}

double trueAnomalyFromMeanAnomalyRadians(double meanAnomaly, double eccentricity) {
    const double eccentricAnomaly = solveEccentricAnomalyRadians(meanAnomaly, eccentricity);
    // tan(ν/2) = sqrt((1+e)/(1-e)) * tan(E/2).
    const double factor = qSqrt((1.0 + eccentricity) / (1.0 - eccentricity));
    return 2.0 * std::atan2(factor * qSin(eccentricAnomaly / 2.0),
                            qCos(eccentricAnomaly / 2.0));
}

QPointF orbitPointAuInclined(double semiMajorAxis,
                             double eccentricity,
                             double perihelionArgumentRadians,
                             double inclinationRadians,
                             double trueAnomalyRadians) {
    // Полярное уравнение эллипса относительно фокуса:
    // r = a * (1 - e^2) / (1 + e * cos(ν)).
    const double radius =
        semiMajorAxis * (1.0 - eccentricity * eccentricity) /
        (1.0 + eccentricity * qCos(trueAnomalyRadians));
    const double xOrbital = radius * qCos(trueAnomalyRadians);
    const double yOrbital = radius * qSin(trueAnomalyRadians);

    // Поворот орбиты на аргумент перицентра (ω).
    const double cosArg = qCos(perihelionArgumentRadians);
    const double sinArg = qSin(perihelionArgumentRadians);
    const double xRotated = xOrbital * cosArg - yOrbital * sinArg;
    const double yRotated = xOrbital * sinArg + yOrbital * cosArg;

    // Наклонение i учитываем как поворот вокруг оси X (проекция на плоскость).
    return QPointF(xRotated, yRotated * qCos(inclinationRadians));
}

QPointF orbitPointFromPhase(double distanceAu, double orbitalPhaseRadians) {
    return QPointF(distanceAu * std::cos(orbitalPhaseRadians),
                   distanceAu * std::sin(orbitalPhaseRadians));
}

RotationMode resolveRotationModeFromPeriods(double dayLengthDays,
                                            double yearLengthDays,
                                            bool hasOverride,
                                            RotationMode manualMode) {
    if (hasOverride) {
        return manualMode;
    }
    const bool isTidallyLocked = isTidallyLockedByPeriods(dayLengthDays, yearLengthDays);
    return isTidallyLocked ? RotationMode::TidalLocked : RotationMode::Normal;
}

double resolveSiderealDayLengthDays(double solarDayLengthDays,
                                    double yearLengthDays,
                                    bool isRetrograde) {
    // Для расчёта часового угла нужна сидерическая длительность суток.
    const double siderealDays =
        solarToSiderealPeriodDays(solarDayLengthDays, yearLengthDays, isRetrograde);
    return (siderealDays > 0.0) ? siderealDays : solarDayLengthDays;
}

double resolveSubstellarLongitudeRadians(double elapsedTimeDays,
                                         double rotationDayLengthDays,
                                         RotationMode rotationMode,
                                         double orbitalPhaseRadians,
                                         int spinOrbitP,
                                         int spinOrbitQ,
                                         bool isRetrograde) {
    if (rotationMode == RotationMode::TidalLocked) {
        return 0.0;
    }
    if (rotationDayLengthDays <= 0.0) {
        return 0.0;
    }
    const int safeSpinOrbitP = qMax(1, spinOrbitP);
    const int safeSpinOrbitQ = qMax(1, spinOrbitQ);
    const double resonanceRatio = static_cast<double>(safeSpinOrbitP) /
        static_cast<double>(safeSpinOrbitQ);
    const double phase = 2.0 * M_PI * (elapsedTimeDays / rotationDayLengthDays);
    // Фаза вращения планеты (λ_spin) выражается через непрерывное время:
    // t берём в земных днях (elapsedDays * solarDayLength), а P_spin — сидерические сутки,
    // чтобы фазу не определяла дискретизация шагов солнечного дня.
    const double spinDirection = isRetrograde ? -1.0 : 1.0;
    // Ретроградное вращение инвертирует направление изменения часового угла.
    const double spinPhase = spinDirection * phase * resonanceRatio;
    const double hourAngle = spinPhase - M_PI;
    // Субзвёздная долгота λ* задаётся разностью орбитальной и спиновой фаз:
    // λ* = λ_orbit - λ_spin, где λ_orbit — истинная долгота звезды (фаза орбиты),
    // а λ_spin = hourAngle — часовой угол на нулевом меридиане. Знак выбран так,
    // чтобы положительная спиновая фаза сдвигала подсолнечную точку к западу.
    return normalizeLongitudeRadians(orbitalPhaseRadians - hourAngle);
}

QVector3D directionFromLatLonRadians(double latitudeRad, double longitudeRad) {
    const double cosLat = std::cos(latitudeRad);
    return QVector3D(static_cast<float>(cosLat * std::sin(longitudeRad)),
                     static_cast<float>(std::sin(latitudeRad)),
                     static_cast<float>(cosLat * std::cos(longitudeRad)));
}

QVector3D resolveStarDirection(double declinationDegrees,
                               double elapsedTimeDays,
                               double rotationDayLengthDays,
                               RotationMode rotationMode,
                               double orbitalPhaseRadians,
                               int spinOrbitP,
                               int spinOrbitQ,
                               bool isRetrograde) {
    const double declinationRadians = qDegreesToRadians(declinationDegrees);
    const double substellarLongitude =
        resolveSubstellarLongitudeRadians(elapsedTimeDays,
                                          rotationDayLengthDays,
                                          rotationMode,
                                          orbitalPhaseRadians,
                                          spinOrbitP,
                                          spinOrbitQ,
                                          isRetrograde);
    // Направление на звезду совпадает с нормалью подсолнечной точки.
    return directionFromLatLonRadians(declinationRadians, substellarLongitude);
}

double estimateSurfaceWaterGigatons(const SurfaceMaterial &material) {
    // В текущем UI нет явного управления гидросферой, поэтому применяем мягкую эвристику:
    // океаническую поверхность считаем водной, остальные материалы — сухими.
    if (material.id == QLatin1String("ocean")) {
        return kEarthWaterGigatons;
    }
    return 0.0;
}

double estimateAirTemperatureKelvin(const SurfacePointState &surfaceState,
                                    double pressureAtm,
                                    double gravity,
                                    double initialAirTemperature,
                                    double blendedInsolation,
                                    double cloudShortwaveTransmission,
                                    RadiationModelType radiationModelType,
                                    const RadiationModel &radiationModel,
                                    bool useAtmosphericModel,
                                    bool manualGreenhouseOnTop,
                                    double timeStepSeconds) {
    if (pressureAtm <= 0.0 || gravity <= 0.0) {
        return surfaceState.temperatureKelvin();
    }
    const double columnMassKgPerM2 =
        (pressureAtm * kStandardPressurePa) / gravity;
    const double airHeatCapacity = columnMassKgPerM2 * kDryAirSpecificHeatJPerKgK;
    if (airHeatCapacity <= 0.0) {
        return surfaceState.temperatureKelvin();
    }

    AtmosphericCellState airState(qMax(0.0, initialAirTemperature), airHeatCapacity);
    // Масштабируем теплообмен по давлению, чтобы плотные атмосферы (например, Венеры)
    // давали более интенсивный обмен, чем разреженные.
    const double couplingScale = qMax(0.0, pressureAtm);
    const double heatTransfer = kDefaultHeatTransferWPerM2K * couplingScale;
    if (timeStepSeconds <= 0.0) {
        return airState.airTemperatureKelvin();
    }

    // Радиационный нагрев воздуха:
    // Q_sw_air = S_blend * T_cloud * (1 - T_atm), где T_atm = incomingTransmission().
    // Длинноволновый баланс берём в двухпотоковой серой аппроксимации:
    // атмосферный слой излучает вверх/вниз как σ T_air^4 * ε(τ).
    const double emittedFlux = surfaceState.surfaceEmittedFlux();
    const double shortwaveAbsorbedByAir =
        blendedInsolation * cloudShortwaveTransmission *
        (1.0 - radiationModel.incomingTransmission());
    const LongwaveFluxes longwaveFluxes =
        computeLongwaveFluxes(radiationModel,
                              radiationModelType,
                              emittedFlux,
                              airState.airTemperatureKelvin(),
                              surfaceState.greenhouseOpacity(),
                              useAtmosphericModel,
                              manualGreenhouseOnTop);
    const double longwaveNetAir =
        emittedFlux - (longwaveFluxes.upwardFluxWPerM2 + longwaveFluxes.downwardFluxWPerM2);
    const double airRadiativeHeatingFlux = shortwaveAbsorbedByAir + longwaveNetAir;
    const double airRadiativeDelta =
        airRadiativeHeatingFlux * timeStepSeconds / airHeatCapacity;
    airState.setAirTemperatureKelvin(airState.airTemperatureKelvin() + airRadiativeDelta);

    if (heatTransfer <= 0.0) {
        return airState.airTemperatureKelvin();
    }

    // Оцениваем изменение температуры воздуха без изменения поверхности,
    // чтобы не модифицировать вычисленное состояние карты.
    const double surfaceTemperature = surfaceState.temperatureKelvin();
    const double airTemperature = airState.airTemperatureKelvin();
    const double fluxWPerM2 = heatTransfer * (surfaceTemperature - airTemperature);
    const double maxStableDt = 0.5 * airHeatCapacity / heatTransfer;
    const double stableDt = qMin(timeStepSeconds, maxStableDt);
    const double airDelta = fluxWPerM2 * stableDt / airHeatCapacity;
    airState.setAirTemperatureKelvin(airTemperature + airDelta);
    return airState.airTemperatureKelvin();
}

double resolveAirTemperatureKelvin(const SurfacePointState &surfaceState,
                                   double pressureAtm,
                                   double gravity,
                                   double initialAirTemperature,
                                   double blendedInsolation,
                                   double cloudShortwaveTransmission,
                                   RadiationModelType radiationModelType,
                                   const RadiationModel &radiationModel,
                                   bool useAtmosphericModel,
                                   bool manualGreenhouseOnTop,
                                   double timeStepSeconds) {
    if (radiationModelType == RadiationModelType::Layered) {
        const auto *layeredModel = dynamic_cast<const LayeredRadiationModel *>(&radiationModel);
        if (layeredModel) {
            return layeredModel->bottomLayerTemperatureKelvin();
        }
    }
    return estimateAirTemperatureKelvin(surfaceState,
                                        pressureAtm,
                                        gravity,
                                        initialAirTemperature,
                                        blendedInsolation,
                                        cloudShortwaveTransmission,
                                        radiationModelType,
                                        radiationModel,
                                        useAtmosphericModel,
                                        manualGreenhouseOnTop,
                                        timeStepSeconds);
}

double computeLocalGreenhouseOpacity(const AtmosphereComposition &atmosphere,
                                     const SurfaceMaterial &material,
                                     double cloudAlbedo,
                                     bool disableWaterAndClouds,
                                     double pressureAtm,
                                     double planetRadiusKm,
                                     double surfaceGravity,
                                     double blendedInsolation,
                                     double meanToaFlux,
                                     double baseTemperatureKelvin,
                                     double manualGreenhouseOpacity,
                                     double minDenseAtmosphereTemperatureK,
                                     bool useAtmosphericModel,
                                     RadiationModelType radiationModelType,
                                     bool manualGreenhouseOnTopOfAtmosphere,
                                     bool logDetails) {
    const double safeRadiusKm = qMax(0.1, planetRadiusKm);
    const double areaScale = std::pow(safeRadiusKm / kEarthRadiusKm, 2.0);

    const double waterGigatons = estimateSurfaceWaterGigatons(material);
    const double planetAreaKm2 = kEarthAreaKm2 * areaScale;
    const double avgDepth = (planetAreaKm2 > 0.0) ? waterGigatons / planetAreaKm2 : 0.0;
    double potentialCoverage = 0.0;
    if (waterGigatons > 0.0) {
        const double fillFactor = (avgDepth * kDefaultBasinShape) / kDefaultSurfaceRoughness;
        potentialCoverage = 1.0 - std::exp(-fillFactor * 3.0);
    }
    potentialCoverage = qBound(0.0, potentialCoverage, 1.0);
    if (disableWaterAndClouds) {
        potentialCoverage = 0.0;
    }

    // Давление также усиливает облачность: давление выше порога повышает отражение.
    const double pressureClouds =
        pressureAtm > 0.05 ? 0.25 * (1.0 - std::exp(-pressureAtm)) : 0.0;
    const double albedo = qBound(0.0, material.albedo, 1.0);
    const double surfAlbedoPre =
        (1.0 - potentialCoverage) * albedo + potentialCoverage * 0.06;
    // Для плотных атмосфер профиль нельзя привязывать к локальной инсоляции:
    // ночью она обнуляется, а радиативный профиль должен оставаться заданным
    // глобальным средним ТОА-потоком.
    // Планетарное альбедо учитываем ровно один раз: meanToaFlux уже скорректирован
    // на облака, поэтому повторно не умножаем на (1 - max(альбедо поверхности, облаков)).
    const double absorbedMeanFlux = meanToaFlux;
    const double tEffPre = std::pow(qMax(0.0, absorbedMeanFlux) / kStefanBoltzmannConstant, 0.25);
    const double atmosphereMassGt = atmosphere.totalMassGigatons();
    const double co2MassGt = atmosphere.massGigatons(QStringLiteral("co2"));
    const double co2Share = atmosphereMassGt > 0.0 ? co2MassGt / atmosphereMassGt : 0.0;
    const double denseAtmosphereThresholdAtm = 0.1;
    double profileBaseTemperatureKelvin = baseTemperatureKelvin;
    if (pressureAtm >= denseAtmosphereThresholdAtm) {
        // Для плотных атмосфер (например, Венера) t_eff не отражает сильный парник:
        // эффективная температура задаётся балансом излучения на верхней границе,
        // но при больших давлениях и непрозрачности нижние слои остаются тёплыми.
        // Поэтому минимальная база зависит от давления и состава.
        const double tMinDense =
            denseAtmosphereMinTemperatureKelvin(tEffPre,
                                                pressureAtm,
                                                manualGreenhouseOpacity,
                                                co2Share,
                                                minDenseAtmosphereTemperatureK);
        profileBaseTemperatureKelvin = qMax(baseTemperatureKelvin, tMinDense);
    }
    // В разреженных атмосферах (Марс/Луна) оставляем прежнюю базу, чтобы не
    // искажать локальные условия.
    const auto preRadiationModel = makeRadiationModel(atmosphere,
                                                      pressureAtm,
                                                      profileBaseTemperatureKelvin,
                                                      tEffPre,
                                                      surfaceGravity,
                                                      radiationModelType);
    const double baseLongwaveTransmission =
        qMax(1e-6, preRadiationModel->outgoingTransmission());
    const double tBasePre =
        tEffPre * std::pow(1.0 / baseLongwaveTransmission, 0.25);

    double evaporation = 0.0;
    if (!disableWaterAndClouds && potentialCoverage > 0.0 && tBasePre > 263.0) {
        evaporation = potentialCoverage * std::exp((tBasePre - 280.0) / 15.0);
    }
    const double waterTau = disableWaterAndClouds ? 0.0 : qMin(8.0, evaporation * 1.5);
    double extraTau = waterTau;
    // Облака повышают альбедо, но также усиливают ИК-поглощение.
    // Держим вклад умеренным и связываем его с облачным альбедо и давлением.
    const double clampedCloudAlbedo = qBound(0.0, cloudAlbedo, 1.0);
    const double cloudLongwaveFactor =
        qBound(0.0, clampedCloudAlbedo + pressureClouds, 1.0);
    constexpr double kCloudLongwaveTau = 0.3;
    const double tauCloudLW =
        disableWaterAndClouds ? 0.0 : kCloudLongwaveTau * cloudLongwaveFactor;
    extraTau += tauCloudLW;
    const double co2PartialPressureAtm = pressureAtm * co2Share;
    if (co2PartialPressureAtm > 0.0) {
        // Серый парник от CO₂: растущее поглощение ограничиваем логарифмом,
        // чтобы имитировать насыщение линий при больших давлениях.
        // k и P0 подобраны эмпирически для мягкого вклада в общий τ.
        constexpr double kTauCo2 = 0.35;
        constexpr double kReferencePressureAtm = 0.1;
        const double tauCo2 =
            kTauCo2 * std::log(1.0 + co2PartialPressureAtm / kReferencePressureAtm);
        extraTau += tauCo2;
    }
    if (manualGreenhouseOpacity > 0.0 &&
        (!useAtmosphericModel || manualGreenhouseOnTopOfAtmosphere)) {
        // Дополнительная непрозрачность: либо без атмосферной модели,
        // либо поверх неё по явному переключателю.
        extraTau += opticalDepthFromGreenhouseOpacity(manualGreenhouseOpacity,
                                                      radiationModelType);
    }

    const double extraLongwaveTransmission =
        longwaveTransmissionForOpticalDepth(extraTau, radiationModelType);
    const double totalLongwaveTransmission =
        qMax(1e-6, baseLongwaveTransmission * extraLongwaveTransmission);

    // Приводим пропускание к коэффициенту парникового эффекта для SurfacePointState.
    const double greenhouseOpacity = 1.0 - totalLongwaveTransmission;
    if (logDetails) {
        qCInfo(solarRadiationLog) << "Local greenhouse opacity"
                                 << "pressureAtm=" << pressureAtm
                                 << "blendedInsolation=" << blendedInsolation
                                 << "surfAlbedoPre=" << surfAlbedoPre
                                 << "pressureClouds=" << pressureClouds
                                 << "tEffPre=" << tEffPre
                                 << "baseLongwaveTransmission=" << baseLongwaveTransmission
                                 << "tBasePre=" << tBasePre
                                 << "evaporation=" << evaporation
                                 << "waterTau=" << waterTau
                                 << "tauCloudLW=" << tauCloudLW
                                 << "extraTau=" << extraTau
                                 << "totalLongwaveTransmission=" << totalLongwaveTransmission
                                 << "greenhouseOpacity=" << greenhouseOpacity;
    }
    return qBound(0.0, greenhouseOpacity, 0.999);
}

struct StellarCacheKey {
    double primaryRadius = 0.0;
    double primaryTemperature = 0.0;
    bool hasSecondary = false;
    double secondaryRadius = 0.0;
    double secondaryTemperature = 0.0;

    bool operator==(const StellarCacheKey &other) const {
        return primaryRadius == other.primaryRadius &&
               primaryTemperature == other.primaryTemperature &&
               hasSecondary == other.hasSecondary &&
               secondaryRadius == other.secondaryRadius &&
               secondaryTemperature == other.secondaryTemperature;
    }

    bool operator!=(const StellarCacheKey &other) const {
        return !(*this == other);
    }
};

struct SurfacePointStateDefaults {
    double albedo = 0.0;
    double greenhouseOpacity = 0.0;
    double manualGreenhouseOpacity = 0.0;
    double minTemperatureKelvin = 3.0;
    double diurnalCoolingBiasK = 0.0;
    RadiationModelType radiationModelType = RadiationModelType::Fast;
    SurfaceMaterial material;
    SubsurfaceModelSettings subsurfaceSettings;
};

struct SurfaceGridComputationInput {
    PlanetSurfaceGrid grid;
    AtmosphereComposition atmosphere;
    SurfacePointStateDefaults stateDefaults;
    QHash<QString, SurfaceMaterial> materialsById;
    QString baseMaterialId;
    double massEarths = 0.0;
    double radiusKm = 0.0;
    double cloudAlbedo = 0.0;
    bool disableWaterAndClouds = false;
    double dayLengthDays = 0.0;
    double yearLengthDays = 0.0;
    double elapsedDays = 0.0;
    int stepsPerDay = 1;
    double manualGreenhouseOpacity = 0.0;
    bool manualGreenhouseOnTop = false;
    RadiationModelType radiationModelType = RadiationModelType::Fast;
    double minDenseAtmosphereTemperatureK = 0.0;
    RotationMode rotationMode = RotationMode::Normal;
    double declinationDegrees = 0.0;
    double orbitalPhaseRadians = 0.0;
    int spinOrbitP = 1;
    int spinOrbitQ = 1;
    bool isRetrograde = false;
    double verticalWindMixingCoefficientKz = 1.0;
    double segmentSolarConstant = 0.0;
    bool hasSolarConstant = false;
    int currentHourIndex = 0;
    int currentAbsoluteHourIndex = 0;
    int logPointIndex = 0;
    bool skipSpinUp = false;
    double minBottomLayerThicknessMeters = 0.0;
    double minTopPressureAtm = 0.0;
};

struct SurfaceGridComputationResult {
    bool cancelled = false;
    PlanetSurfaceGrid grid;
    SubsurfaceModelSettings resolvedSubsurfaceSettings;
    double dayLengthDays = 0.0;
    double yearLengthDays = 0.0;
    double elapsedDays = 0.0;
    bool hasTemperatureRange = false;
    double minTemperatureK = 0.0;
    double maxTemperatureK = 0.0;
    bool hasAirTemperatureRange = false;
    double minAirTemperatureK = 0.0;
    double maxAirTemperatureK = 0.0;
    bool hasPressureRange = false;
    double minPressureAtm = 0.0;
    double maxPressureAtm = 0.0;
    bool hasWindRange = false;
    double minWindSpeedMps = 0.0;
    double maxWindSpeedMps = 0.0;
    double declinationDegrees = 0.0;
    double orbitalPhaseRadians = 0.0;
    int stepsPerDay = 1;
    RotationMode rotationMode = RotationMode::Normal;
    int spinOrbitP = 1;
    int spinOrbitQ = 1;
    bool isRetrograde = false;
};

class SolarCalculatorWidget : public QWidget {
public:
    explicit SolarCalculatorWidget(int precision, bool enableRadiationLog, QWidget *parent = nullptr)
        : QWidget(parent), precision_(precision) {
        auto *validator = new QDoubleValidator(0.0, std::numeric_limits<double>::max(), 10, this);
        validator->setNotation(QDoubleValidator::StandardNotation);

        radiusInput_ = new QLineEdit(this);
        radiusInput_->setPlaceholderText(QStringLiteral("Например, 1.0"));
        radiusInput_->setValidator(validator);

        temperatureInput_ = new QLineEdit(this);
        temperatureInput_->setPlaceholderText(QStringLiteral("Например, 5772"));
        temperatureInput_->setValidator(validator);

        auto *presetsLayout = new QHBoxLayout();

        const auto applyPrimary = [this](const StellarParameters &parameters) {
            setInputValue(radiusInput_, parameters.radiusInSolarRadii);
            setInputValue(temperatureInput_, parameters.temperatureKelvin);
        };

        const auto applySecondary = [this](const std::optional<StellarParameters> &parameters) {
            if (!parameters) {
                secondStarCheckBox_->setChecked(false);
                secondaryRadiusInput_->clear();
                secondaryTemperatureInput_->clear();
                return;
            }

            secondStarCheckBox_->setChecked(true);
            setInputValue(secondaryRadiusInput_, parameters->radiusInSolarRadii);
            setInputValue(secondaryTemperatureInput_, parameters->temperatureKelvin);
        };

        const auto addPresetButton = [this, presetsLayout](const QString &text,
                                                           const std::function<void()> &handler) {
            auto *button = new QPushButton(text, this);
            connect(button, &QPushButton::clicked, this, handler);
            presetsLayout->addWidget(button);
        };

        addPresetButton(QStringLiteral("Солнце"), [this]() {
            resetSolarConstant();
            const auto presets = solarSystemPresets();
            QString selectedPlanet = QStringLiteral("Земля");
            const bool hasEarth =
                std::any_of(presets.begin(), presets.end(), [](const PlanetPreset &preset) {
                    return preset.name == QStringLiteral("Земля");
                });
            if (!hasEarth && !presets.isEmpty()) {
                selectedPlanet = presets.first().name;
            }
            setPlanetPresets(presets, selectedPlanet);
            autoCalculateEnabled_ = true;
            onCalculateRequested();
        });

        addPresetButton(QStringLiteral("Сладкое Небо"), [this]() {
            resetSolarConstant();
            const auto presets = sweetSkyPresets();
            const QString selectedPlanet =
                presets.isEmpty() ? QString() : presets.first().name;
            setPlanetPresets(presets, selectedPlanet);
            autoCalculateEnabled_ = true;
            onCalculateRequested();
        });

        addPresetButton(QStringLiteral("Пусто"), [this, applyPrimary, applySecondary]() {
            radiusInput_->clear();
            temperatureInput_->clear();
            resetSolarConstant();
            clearPlanetPresets();
            applySecondary(std::nullopt);
        });

        auto *primaryFormLayout = new QFormLayout();
        primaryFormLayout->addRow(QStringLiteral("Радиус звезды (в R☉):"), radiusInput_);
        primaryFormLayout->addRow(QStringLiteral("Температура поверхности (K):"), temperatureInput_);
        primaryGroupBox_ = new QGroupBox(QStringLiteral("Звезда 1"), this);
        primaryGroupBox_->setLayout(primaryFormLayout);

        secondStarCheckBox_ = new QCheckBox(QStringLiteral("Добавить вторую звезду"), this);

        secondaryRadiusInput_ = new QLineEdit(this);
        secondaryRadiusInput_->setPlaceholderText(QStringLiteral("Например, 0.9"));
        secondaryRadiusInput_->setValidator(validator);

        secondaryTemperatureInput_ = new QLineEdit(this);
        secondaryTemperatureInput_->setPlaceholderText(QStringLiteral("Например, 5200"));
        secondaryTemperatureInput_->setValidator(validator);

        auto *secondaryFormLayout = new QFormLayout();
        secondaryFormLayout->addRow(QStringLiteral("Радиус второй звезды (в R☉):"), secondaryRadiusInput_);
        secondaryFormLayout->addRow(QStringLiteral("Температура второй звезды (K):"), secondaryTemperatureInput_);

        secondaryGroupBox_ = new QGroupBox(QStringLiteral("Звезда 2"), this);
        secondaryGroupBox_->setLayout(secondaryFormLayout);
        secondaryGroupBox_->setEnabled(false);
        secondaryGroupBox_->setVisible(false);

        connect(secondStarCheckBox_, &QCheckBox::toggled, this, [this](bool checked) {
            secondaryGroupBox_->setEnabled(checked);
            secondaryGroupBox_->setVisible(checked);
        });

        auto *calculateButton = new QPushButton(QStringLiteral("Рассчитать"), this);

        resultLabel_ = new QLabel(
            QStringLiteral("Введите параметры и нажмите \"Рассчитать\"."), this);
        resultLabel_->setWordWrap(true);

        planetComboBox_ = new QComboBox(this);
        planetSemiMajorAxisLabel_ = new QLabel(QStringLiteral("—"), this);
        planetDayLengthLabel_ = new QLabel(QStringLiteral("—"), this);
        planetYearLengthLabel_ = new QLabel(QStringLiteral("—"), this);
        planetMassLabel_ = new QLabel(QStringLiteral("—"), this);
        planetRadiusLabel_ = new QLabel(QStringLiteral("—"), this);
        planetSurfaceGravityLabel_ = new QLabel(QStringLiteral("—"), this);
        planetSurfaceAreaLabel_ = new QLabel(QStringLiteral("—"), this);
        planetEccentricityLabel_ = new QLabel(QStringLiteral("—"), this);
        planetObliquityLabel_ = new QLabel(QStringLiteral("—"), this);
        planetPerihelionArgumentLabel_ = new QLabel(QStringLiteral("—"), this);
        materialComboBox_ = new QComboBox(this);
        populateMaterials();
        rotationModeComboBox_ = new QComboBox(this);
        rotationModeComboBox_->addItem(QStringLiteral("Обычное вращение (широта)"),
                                       static_cast<int>(RotationMode::Normal));
        rotationModeComboBox_->addItem(
            QStringLiteral("Приливная синхронизация (угол от подсолнечной точки)"),
            static_cast<int>(RotationMode::TidalLocked));
        spinOrbitPSpinBox_ = new QSpinBox(this);
        spinOrbitPSpinBox_->setRange(1, 20);
        spinOrbitPSpinBox_->setValue(1);
        spinOrbitQSpinBox_ = new QSpinBox(this);
        spinOrbitQSpinBox_->setRange(1, 20);
        spinOrbitQSpinBox_->setValue(1);
        const QString spinOrbitTooltip =
            QStringLiteral("Резонанс p:q — отношение спиновой и орбитальной частот.\n"
                           "Безразмерные коэффициенты, 1:1 соответствует синхронному вращению.");
        spinOrbitPSpinBox_->setToolTip(spinOrbitTooltip);
        spinOrbitQSpinBox_->setToolTip(spinOrbitTooltip);
        heightSeedSpinBox_ = new QSpinBox(this);
        heightSeedSpinBox_->setRange(0, std::numeric_limits<int>::max());
        heightmapImportButton_ = new QPushButton(QStringLiteral("Импорт heightmap..."), this);
        heightmapScaleSpinBox_ = new QDoubleSpinBox(this);
        heightmapScaleSpinBox_->setRange(0.0, 100.0);
        heightmapScaleSpinBox_->setDecimals(3);
        heightmapScaleSpinBox_->setSingleStep(0.1);
        heightmapScaleSpinBox_->setSuffix(QStringLiteral(" км"));
        heightmapScaleSpinBox_->setToolTip(
            QStringLiteral("Масштаб рельефа: 1.0 = диапазон (-1..+1) км."));
        flatHeightButton_ = new QPushButton(QStringLiteral("Обнулить высоту"), this);
        flatHeightButton_->setCheckable(true);
        flatHeightButton_->setEnabled(false);
        updateFlatHeightButtonText(false);
        disableAtmosphereButton_ = new QPushButton(QStringLiteral("Отключить атмосферу"), this);
        disableAtmosphereButton_->setCheckable(true);
        disableAtmosphereButton_->setEnabled(false);
        updateAtmosphereButtonText(false);
        cloudAlbedoSpinBox_ = new QDoubleSpinBox(this);
        cloudAlbedoSpinBox_->setRange(0.0, 1.0);
        cloudAlbedoSpinBox_->setDecimals(2);
        cloudAlbedoSpinBox_->setSingleStep(0.05);
        diurnalCoolingBiasSpinBox_ = new QDoubleSpinBox(this);
        diurnalCoolingBiasSpinBox_->setRange(0.0, 200.0);
        diurnalCoolingBiasSpinBox_->setDecimals(1);
        diurnalCoolingBiasSpinBox_->setSingleStep(1.0);
        diurnalCoolingBiasSpinBox_->setSuffix(QStringLiteral(" K"));
        diurnalCoolingBiasSpinBox_->setToolTip(
            QStringLiteral("Грубая поправка на суточное охлаждение поверхности.\n"
                           "Сдвигает стартовую температуру вниз при инициализации."));
        manualGreenhouseOnTopCheckBox_ = new QCheckBox(
            QStringLiteral("Добавлять ручную парниковую непрозрачность поверх атмосферы"), this);
        manualGreenhouseOnTopCheckBox_->setToolTip(
            QStringLiteral("Если включено, значение парниковой непрозрачности добавляется\n"
                           "к атмосферной модели и используется как дополнительный фактор."));
        advancedRadiationCheckBox_ = new QCheckBox(
            QStringLiteral("Точная радиационная модель (многослойная)"), this);
        advancedRadiationCheckBox_->setToolTip(
            QStringLiteral("Уточняет передачу коротковолнового и длинноволнового излучения\n"
                           "через атмосферу в многослойной аппроксимации.\n"
                           "По умолчанию точная модель включена."));
        debugLogCheckBox_ = new QCheckBox(
            QStringLiteral("Показывать отладочные логи в консоли"), this);
        debugLogCheckBox_->setChecked(enableRadiationLog);
        debugLogCheckBox_->setToolTip(
            QStringLiteral("Выводит подробные расчётные параметры в консоль/дебаг.\n"
                           "По умолчанию вывод отключён."));
        debugLogEnabled_ = debugLogCheckBox_->isChecked();
        setSolarDebugLoggingEnabled(debugLogEnabled_);
        modeIllustrationWidget_ = new ModeIllustrationWidget(this);
        modeIllustrationWidget_->setRotationMode(
            static_cast<RotationMode>(rotationModeComboBox_->currentData().toInt()));
        addPlanetButton_ = new QPushButton(QStringLiteral("Добавить"), this);
        deletePlanetButton_ = new QPushButton(this);
        deletePlanetButton_->setIcon(style()->standardIcon(QStyle::SP_TrashIcon));
        deletePlanetButton_->setToolTip(QStringLiteral("Удалить планету"));
        deletePlanetButton_->setVisible(false);

        auto *planetSelectorLayout = new QHBoxLayout();
        planetSelectorLayout->addWidget(planetComboBox_);
        planetSelectorLayout->addWidget(deletePlanetButton_);

        auto *planetHeaderLayout = new QFormLayout();
        planetHeaderLayout->addRow(QStringLiteral("Планета:"), planetSelectorLayout);

        auto *planetLeftFormLayout = new QFormLayout();
        planetLeftFormLayout->addRow(QStringLiteral("Большая полуось (а.е.):"), planetSemiMajorAxisLabel_);
        planetLeftFormLayout->addRow(QStringLiteral("Длина суток (земн. дни):"), planetDayLengthLabel_);
        planetLeftFormLayout->addRow(QStringLiteral("Длина года (земн. дни):"), planetYearLengthLabel_);
        planetLeftFormLayout->addRow(QStringLiteral("Масса (M⊕):"), planetMassLabel_);
        planetLeftFormLayout->addRow(QStringLiteral("Радиус (км):"), planetRadiusLabel_);
        planetLeftFormLayout->addRow(QStringLiteral("G на поверхности (g⊕):"),
                                     planetSurfaceGravityLabel_);
        planetLeftFormLayout->addRow(QStringLiteral("Площадь поверхности (км²):"),
                                     planetSurfaceAreaLabel_);

        auto *planetRightFormLayout = new QFormLayout();
        planetRightFormLayout->addRow(QStringLiteral("Эксцентриситет орбиты:"), planetEccentricityLabel_);
        planetRightFormLayout->addRow(QStringLiteral("Наклон оси (°):"), planetObliquityLabel_);
        planetRightFormLayout->addRow(QStringLiteral("Аргумент перицентра (°):"),
                                      planetPerihelionArgumentLabel_);
        auto *rotationModeLayout = new QHBoxLayout();
        rotationModeLayout->addWidget(rotationModeComboBox_);
        auto *rotationModeWidget = new QWidget(this);
        rotationModeWidget->setLayout(rotationModeLayout);
        auto *spinOrbitWidget = new QWidget(this);
        auto *spinOrbitLayout = new QHBoxLayout(spinOrbitWidget);
        spinOrbitLayout->addWidget(spinOrbitPSpinBox_);
        spinOrbitLayout->addWidget(new QLabel(QStringLiteral(":"), spinOrbitWidget));
        spinOrbitLayout->addWidget(spinOrbitQSpinBox_);
        spinOrbitLayout->addStretch();
        spinOrbitLayout->setContentsMargins(0, 0, 0, 0);
        auto *planetColumnsLayout = new QHBoxLayout();
        planetColumnsLayout->addLayout(planetLeftFormLayout);
        planetColumnsLayout->addLayout(planetRightFormLayout);

        auto *planetControlsLayout = new QFormLayout();
        planetControlsLayout->addRow(QStringLiteral("Материал поверхности:"), materialComboBox_);
        planetControlsLayout->addRow(QStringLiteral("Режим вращения:"), rotationModeWidget);
        planetControlsLayout->addRow(QStringLiteral("Резонанс вращения p:q:"), spinOrbitWidget);
        planetControlsLayout->addRow(QStringLiteral("Семя рельефа:"), heightSeedSpinBox_);
        auto *heightmapControls = new QWidget(this);
        auto *heightmapControlsLayout = new QHBoxLayout(heightmapControls);
        heightmapControlsLayout->setContentsMargins(0, 0, 0, 0);
        heightmapControlsLayout->addWidget(heightmapImportButton_);
        heightmapControlsLayout->addWidget(heightmapScaleSpinBox_);
        heightmapControlsLayout->addStretch();
        planetControlsLayout->addRow(QStringLiteral("Heightmap:"), heightmapControls);
        auto *heightAtmosphereButtons = new QWidget(this);
        auto *heightAtmosphereLayout = new QHBoxLayout(heightAtmosphereButtons);
        heightAtmosphereLayout->setContentsMargins(0, 0, 0, 0);
        heightAtmosphereLayout->addWidget(flatHeightButton_);
        heightAtmosphereLayout->addWidget(disableAtmosphereButton_);
        heightAtmosphereLayout->addStretch();
        planetControlsLayout->addRow(QStringLiteral("Высота поверхности:"), heightAtmosphereButtons);
        planetControlsLayout->addRow(QStringLiteral("Альбедо облаков (0..1):"), cloudAlbedoSpinBox_);
        planetControlsLayout->addRow(QStringLiteral("Поправка суточного охлаждения (К):"),
                                     diurnalCoolingBiasSpinBox_);
        planetControlsLayout->addRow(QString(), manualGreenhouseOnTopCheckBox_);
        planetControlsLayout->addRow(QString(), advancedRadiationCheckBox_);
        planetControlsLayout->addRow(QString(), debugLogCheckBox_);
        planetControlsLayout->addRow(QStringLiteral("Солнечная постоянная (Вт/м²):"), resultLabel_);

        auto *planetActionsLayout = new QHBoxLayout();
        planetActionsLayout->addStretch();
        planetActionsLayout->addWidget(addPlanetButton_);

        auto *planetFormLayout = new QVBoxLayout();
        planetFormLayout->addLayout(planetHeaderLayout);
        planetFormLayout->addLayout(planetColumnsLayout);
        planetFormLayout->addLayout(planetControlsLayout);
        planetFormLayout->addLayout(planetActionsLayout);
        auto *planetGroupBox = new QGroupBox(QStringLiteral("Планеты"), this);
        planetGroupBox->setLayout(planetFormLayout);

        atmosphereWidget_ = new AtmosphereWidget(this, false);
        atmosphereWidget_->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        connect(atmosphereWidget_, &AtmosphereWidget::compositionChanged, this,
                [this](const AtmosphereComposition &composition) {
                    if (!planetComboBox_) {
                        return;
                    }
                    const int currentIndex = planetComboBox_->currentIndex();
                    if (currentIndex < 0) {
                        return;
                    }
                    planetComboBox_->setItemData(
                        currentIndex, QVariant::fromValue(composition), kRoleAtmosphere);
                    updateSurfaceGridTemperatures();
                });
        connect(atmosphereWidget_, &AtmosphereWidget::minBottomLayerThicknessChanged, this,
                [this](double meters) {
                    if (!planetComboBox_) {
                        return;
                    }
                    const int currentIndex = planetComboBox_->currentIndex();
                    if (currentIndex < 0) {
                        return;
                    }
                    planetComboBox_->setItemData(currentIndex,
                                                 meters,
                                                 kRoleAtmosphereBottomLayerThickness);
                    updateSurfaceGridTemperatures();
                });
        connect(atmosphereWidget_, &AtmosphereWidget::minTopPressureChanged, this,
                [this](double minTopPressureAtm) {
                    if (!planetComboBox_) {
                        return;
                    }
                    const int currentIndex = planetComboBox_->currentIndex();
                    if (currentIndex < 0) {
                        return;
                    }
                    planetComboBox_->setItemData(currentIndex,
                                                 minTopPressureAtm,
                                                 kRoleMinTopPressureAtm);
                    updateSurfaceGridTemperatures();
                });
        connect(atmosphereWidget_, &AtmosphereWidget::verticalWindMixingCoefficientChanged, this,
                [this](double coefficient) {
                    if (!planetComboBox_) {
                        return;
                    }
                    const int currentIndex = planetComboBox_->currentIndex();
                    if (currentIndex < 0) {
                        return;
                    }
                    planetComboBox_->setItemData(currentIndex,
                                                 coefficient,
                                                 kRoleVerticalWindMixingCoefficient);
                    updateSurfaceGridTemperatures();
                });

        auto *starsPanelLayout = new QVBoxLayout();
        starsPanelLayout->addWidget(primaryGroupBox_);
        starsPanelLayout->addWidget(secondStarCheckBox_);
        starsPanelLayout->addWidget(secondaryGroupBox_);
        auto *starsPanel = new QGroupBox(QStringLiteral("Панель звезд"), this);
        starsPanel->setLayout(starsPanelLayout);

        temperaturePlot_ = new SurfaceTemperaturePlot(this);
        starSystemTopView_ = new StarSystemTopView(this);
        starSystemTopView_->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        connect(starSystemTopView_, &StarSystemTopView::planetSelected, this, [this](int index) {
            if (!planetComboBox_ || index < 0 || index >= planetComboBox_->count()) {
                return;
            }
            planetComboBox_->setCurrentIndex(index);
        });
        surfaceMapWidget_ = new SurfaceMapWidget(this);
        surfaceMapWidget_->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        surfaceMapWidget_->setGrid(&surfaceGrid_);
        surfaceGlobeWidget_ = new SurfaceGlobeWidget(this);
        surfaceGlobeWidget_->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        surfaceGlobeWidget_->setGrid(&surfaceGrid_);
        connect(surfaceGlobeWidget_, &SurfaceGlobeWidget::pointClicked, this,
                [this](int pointIndex) {
                    const SurfacePoint *point = surfaceGrid_.pointAt(pointIndex);
                    if (!point) {
                        return;
                    }
                    if (!surfacePointStatusDialog_) {
                        surfacePointStatusDialog_ = new SurfacePointStatusDialog(this);
                    }
                    selectedSurfacePointIndex_ = pointIndex;
                    const double tileAreaKm2 = surfaceGrid_.tileAreaKm2(pointIndex);
                    const double tileEdgeLengthKm = surfaceGrid_.tileEdgeLengthKm(pointIndex);
                    surfacePointStatusDialog_->setPoint(*point, tileAreaKm2, tileEdgeLengthKm);
                    surfacePointStatusDialog_->show();
                    surfacePointStatusDialog_->raise();
                    surfacePointStatusDialog_->activateWindow();
                });
        surfaceViewStack_ = new QStackedWidget(this);
        surfaceViewStack_->addWidget(surfaceMapWidget_);
        surfaceViewStack_->addWidget(surfaceGlobeWidget_);
        surfaceViewStack_->setCurrentWidget(surfaceGlobeWidget_);
        surfaceMinTemperatureLabel_ = new QLabel(QStringLiteral("Мин: —"), this);
        surfaceMaxTemperatureLabel_ = new QLabel(QStringLiteral("Макс: —"), this);
        surfaceSolarConstantLabel_ = new QLabel(QStringLiteral("S0: —"), this);
        surfaceSolarConstantLabel_->setToolTip(
            QStringLiteral("Солнечная постоянная, пересчитанная на текущее расстояние."));
        temperaturePauseButton_ = new QPushButton(QStringLiteral("Пауза"), this);
        surfaceSimToggleButton_ = new QPushButton(QStringLiteral("Течение времени"), this);
        surfaceSimToggleButton_->setToolTip(
            QStringLiteral("Запустить или остановить течение времени"));
        surfaceSimSpeedComboBox_ = new QComboBox(this);
        surfaceSimSpeedComboBox_->addItem(QStringLiteral("1x"), 1.0);
        surfaceSimSpeedComboBox_->addItem(QStringLiteral("10x"), 10.0);
        surfaceSimSpeedComboBox_->addItem(QStringLiteral("Без ограничений"), 0.0);
        surfaceSimTimeLabel_ = new QLabel(QStringLiteral("t = —"), this);
        surfaceMapModeComboBox_ = new QComboBox(this);
        surfaceMapModeComboBox_->addItem(QStringLiteral("Температура поверхности"),
                                         static_cast<int>(SurfaceMapMode::Temperature));
        surfaceMapModeComboBox_->addItem(QStringLiteral("Температура нижнего слоя атмосферы"),
                                         static_cast<int>(SurfaceMapMode::AirTemperature));
        surfaceMapModeComboBox_->addItem(QStringLiteral("Высота"),
                                         static_cast<int>(SurfaceMapMode::Height));
        surfaceMapModeComboBox_->addItem(QStringLiteral("Ветер"),
                                         static_cast<int>(SurfaceMapMode::Wind));
        surfaceMapModeComboBox_->addItem(QStringLiteral("Давление"),
                                         static_cast<int>(SurfaceMapMode::Pressure));
        surfaceMapModeComboBox_->addItem(QStringLiteral("Реалистичный"),
                                         static_cast<int>(SurfaceMapMode::Realistic));
        subsurfaceLayersSpinBox_ = new QSpinBox(this);
        subsurfaceLayersSpinBox_->setRange(1, 200);
        subsurfaceLayersSpinBox_->setValue(24);
        surfaceViewToggleButton_ = new QPushButton(QStringLiteral("3D вид"), this);
        surfaceViewToggleButton_->setCheckable(true);
        surfaceMarkupCheckBox_ = new QCheckBox(QStringLiteral("Разметка"), this);
        surfaceMarkupCheckBox_->setChecked(true);
        surfaceGlobeWidget_->setMarkupVisible(surfaceMarkupCheckBox_->isChecked());
        connect(surfaceMarkupCheckBox_, &QCheckBox::toggled, this, [this](bool checked) {
            if (surfaceGlobeWidget_) {
                surfaceGlobeWidget_->setMarkupVisible(checked);
            }
        });
        subsurfaceTopThicknessSpinBox_ = new QDoubleSpinBox(this);
        subsurfaceTopThicknessSpinBox_->setRange(0.001, 10.0);
        subsurfaceTopThicknessSpinBox_->setDecimals(3);
        subsurfaceTopThicknessSpinBox_->setSingleStep(0.01);
        subsurfaceTopThicknessSpinBox_->setValue(0.05);
        subsurfaceTopThicknessSpinBox_->setSuffix(QStringLiteral(" м"));
        subsurfaceDepthSpinBox_ = new QDoubleSpinBox(this);
        subsurfaceDepthSpinBox_->setRange(0.1, 200.0);
        subsurfaceDepthSpinBox_->setDecimals(2);
        subsurfaceDepthSpinBox_->setSingleStep(0.1);
        subsurfaceDepthSpinBox_->setValue(2.0);
        subsurfaceDepthSpinBox_->setSuffix(QStringLiteral(" м"));
        subsurfaceGeothermalFluxSpinBox_ = new QDoubleSpinBox(this);
        subsurfaceGeothermalFluxSpinBox_->setRange(0.0, 1.0);
        subsurfaceGeothermalFluxSpinBox_->setDecimals(4);
        subsurfaceGeothermalFluxSpinBox_->setSingleStep(0.005);
        // Типичный диапазон геотермального потока для безатмосферных тел: 0.01–0.02 Вт/м².
        subsurfaceGeothermalFluxSpinBox_->setValue(0.015);
        subsurfaceGeothermalFluxSpinBox_->setSuffix(QStringLiteral(" Вт/м²"));
        subsurfaceGeothermalFluxSpinBox_->setToolTip(
            QStringLiteral("Типичный диапазон для безатмосферных тел: 0.01–0.02 Вт/м²."));
        subsurfaceBoundaryComboBox_ = new QComboBox(this);
        subsurfaceBoundaryComboBox_->addItem(QStringLiteral("Поток = 0"),
                                             static_cast<int>(SubsurfaceBottomBoundaryCondition::Insulating));
        subsurfaceBoundaryComboBox_->addItem(QStringLiteral("Фиксированная температура"),
                                             static_cast<int>(SubsurfaceBottomBoundaryCondition::FixedTemperature));
        temperatureScaleWidget_ = new SurfaceTemperatureScaleWidget(this);
        temperatureScaleWidget_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        temperatureScaleWidget_->setMinimumHeight(18);
        heightScaleWidget_ = new SurfaceHeightScaleWidget(this);
        heightScaleWidget_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        heightScaleWidget_->setMinimumHeight(18);
        windScaleWidget_ = new SurfaceWindScaleWidget(this);
        windScaleWidget_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        windScaleWidget_->setMinimumHeight(18);
        pressureScaleWidget_ = new SurfacePressureScaleWidget(this);
        pressureScaleWidget_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        pressureScaleWidget_->setMinimumHeight(18);
        realisticScaleWidget_ = new SurfaceRealisticScaleWidget(this);
        realisticScaleWidget_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        realisticScaleWidget_->setMinimumHeight(18);
        surfaceLegendScaleStack_ = new QStackedWidget(this);
        surfaceLegendScaleStack_->addWidget(temperatureScaleWidget_);
        surfaceLegendScaleStack_->addWidget(heightScaleWidget_);
        surfaceLegendScaleStack_->addWidget(windScaleWidget_);
        surfaceLegendScaleStack_->addWidget(pressureScaleWidget_);
        surfaceLegendScaleStack_->addWidget(realisticScaleWidget_);
        surfaceLegendScaleStack_->setCurrentWidget(temperatureScaleWidget_);
        auto *surfaceLegendTopLayout = new QHBoxLayout();
        surfaceLegendTopLayout->addWidget(surfaceMinTemperatureLabel_);
        surfaceLegendTopLayout->addStretch();
        surfaceLegendTopLayout->addWidget(surfaceSolarConstantLabel_);
        surfaceLegendTopLayout->addStretch();
        surfaceLegendTopLayout->addWidget(surfaceMaxTemperatureLabel_);
        auto *surfaceControlLayout = new QHBoxLayout();
        // Управление дублирует диалог, чтобы расчётом можно было управлять без модального окна.
        surfaceControlLayout->addWidget(temperaturePauseButton_);
        surfaceControlLayout->addWidget(new QLabel(QStringLiteral("Карта:"), this));
        surfaceControlLayout->addWidget(surfaceMapModeComboBox_);
        surfaceControlLayout->addWidget(surfaceViewToggleButton_);
        surfaceControlLayout->addWidget(surfaceMarkupCheckBox_);
        surfaceControlLayout->addStretch();
        surfaceCalculationLabel_ = new QLabel(QStringLiteral("Пересчёт поверхности..."), this);
        surfaceCalculationProgress_ = new QProgressBar(this);
        surfaceCalculationProgress_->setRange(0, 0);
        surfaceCalculationProgress_->setTextVisible(false);
        surfaceCalculationLabel_->setVisible(false);
        surfaceCalculationProgress_->setVisible(false);
        auto *surfaceCalculationLayout = new QHBoxLayout();
        surfaceCalculationLayout->addWidget(surfaceCalculationLabel_);
        surfaceCalculationLayout->addWidget(surfaceCalculationProgress_, 1);
        surfaceCalculationLayout->addStretch();
        auto *surfaceLegendBottomLayout = new QHBoxLayout();
        surfaceLegendBottomLayout->addStretch();
        surfaceLegendBottomLayout->addWidget(surfaceLegendScaleStack_, 1);
        surfaceLegendBottomLayout->addStretch();
        auto *subsurfaceFormLayout = new QFormLayout();
        subsurfaceFormLayout->addRow(QStringLiteral("Слои подповерхности:"), subsurfaceLayersSpinBox_);
        subsurfaceFormLayout->addRow(QStringLiteral("Верхняя толщина:"), subsurfaceTopThicknessSpinBox_);
        subsurfaceFormLayout->addRow(QStringLiteral("Глубина модели:"), subsurfaceDepthSpinBox_);
        subsurfaceFormLayout->addRow(QStringLiteral("Граница снизу:"), subsurfaceBoundaryComboBox_);
        subsurfaceFormLayout->addRow(QStringLiteral("Геотермальный поток:"),
                                     subsurfaceGeothermalFluxSpinBox_);
        auto *subsurfaceGroupBox = new QGroupBox(QStringLiteral("Подповерхностная модель"), this);
        subsurfaceGroupBox->setLayout(subsurfaceFormLayout);
        subsurfaceGroupBox->setVisible(false);
        auto *subsurfaceToggleButton = new QToolButton(this);
        subsurfaceToggleButton->setText(QStringLiteral("Подповерхностная модель"));
        subsurfaceToggleButton->setCheckable(true);
        subsurfaceToggleButton->setChecked(false);
        subsurfaceToggleButton->setArrowType(Qt::RightArrow);
        subsurfaceToggleButton->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
        connect(subsurfaceToggleButton, &QToolButton::toggled, this,
                [subsurfaceGroupBox, subsurfaceToggleButton](bool checked) {
                    subsurfaceGroupBox->setVisible(checked);
                    subsurfaceToggleButton->setArrowType(checked ? Qt::DownArrow
                                                                : Qt::RightArrow);
                });
        auto *surfaceMapLayout = new QVBoxLayout();
        surfaceMapLayout->addLayout(surfaceLegendTopLayout);
        surfaceMapLayout->addLayout(surfaceControlLayout);
        surfaceMapLayout->addLayout(surfaceCalculationLayout);
        surfaceMapLayout->addWidget(subsurfaceToggleButton);
        surfaceMapLayout->addWidget(subsurfaceGroupBox);
        surfaceMapLayout->addWidget(surfaceViewStack_, 1);
        surfaceMapLayout->addLayout(surfaceLegendBottomLayout);
        surfaceMapContainer_ = new QWidget(this);
        surfaceMapContainer_->setLayout(surfaceMapLayout);
        auto *plotGroupBox = new QGroupBox(QStringLiteral("Температурный профиль"), this);
        auto *plotLayout = new QVBoxLayout(plotGroupBox);
        auto *segmentLayout = new QHBoxLayout();
        segmentSelectorWidget_ = new SegmentSelectorWidget(plotGroupBox);
        segmentSelectorWidget_->setEnabled(false);
        segmentLayout->addWidget(new QLabel(QStringLiteral("Сегмент орбиты:"), plotGroupBox));
        segmentLayout->addWidget(segmentSelectorWidget_, 1);
        plotLayout->addLayout(segmentLayout);
        temperaturePlotSourceComboBox_ = new QComboBox(plotGroupBox);
        temperaturePlotSourceComboBox_->addItem(QStringLiteral("Поверхность (сетка)"),
                                                static_cast<int>(TemperaturePlotSource::SurfaceGrid));
        // temperaturePlotSourceComboBox_->addItem(QStringLiteral("1D модель по широте"),
        //                                         static_cast<int>(TemperaturePlotSource::LatitudeModel1D));
        // Широтная модель устарела и не отражает ночную сторону при приливном захвате.
        auto *plotSourceLayout = new QHBoxLayout();
        plotSourceLayout->addWidget(new QLabel(QStringLiteral("Источник:"), plotGroupBox));
        plotSourceLayout->addWidget(temperaturePlotSourceComboBox_);
        plotLayout->addLayout(plotSourceLayout);
        // auto *smoothingCheckBox = new QCheckBox(QStringLiteral("Сглаживать график"), plotGroupBox);
        // plotLayout->addWidget(smoothingCheckBox);
        plotLayout->addWidget(temperaturePlot_);
        plotLayout->addWidget(modeIllustrationWidget_, 0, Qt::AlignHCenter);
        plotGroupBox->setLayout(plotLayout);
        plotGroupBox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        connect(temperaturePlotSourceComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int) { ensureTemperaturePlotSourceSurfaceGrid(); });

        auto *leftLayout = new QVBoxLayout();
        leftLayout->addLayout(presetsLayout);
        leftLayout->addWidget(starsPanel);
        leftLayout->addWidget(planetGroupBox);
        leftLayout->addWidget(calculateButton);
        leftLayout->addStretch();

        rightTabs_ = new QTabWidget(this);
        rightTabs_->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        rightTabs_->addTab(starSystemTopView_, tr("Система"));
        rightTabs_->addTab(plotGroupBox, tr("Температура"));
        rightTabs_->addTab(atmosphereWidget_, tr("Атмосфера"));
        rightTabs_->addTab(surfaceMapContainer_, tr("Поверхность"));
        connect(rightTabs_, &QTabWidget::currentChanged, this,
                [this](int) {
                    updateStarSystemTimerState();
                    if (!isSurfaceTabActive()) {
                        updateSurfaceCalculationIndicator();
                        return;
                    }
                    // Пересчитываем поверхность только при первом открытии вкладки.
                    if (surfaceGrid_.points().isEmpty() || surfaceGridDirty_) {
                        updateSurfaceGridTemperatures();
                    }
                    updateSurfaceCalculationIndicator();
                });

        auto *rightLayout = new QVBoxLayout();
        auto *timeFlowLayout = new QHBoxLayout();
        timeFlowLayout->addWidget(surfaceSimToggleButton_);
        timeFlowLayout->addWidget(surfaceSimSpeedComboBox_);
        timeFlowLayout->addWidget(surfaceSimTimeLabel_);
        timeFlowLayout->addStretch();
        rightLayout->addLayout(timeFlowLayout);
        rightLayout->addWidget(rightTabs_, 1);

        auto *layout = new QHBoxLayout(this);
        layout->addLayout(leftLayout, 0);
        layout->addLayout(rightLayout, 1);

        setLayout(layout);
        resize(480, 360);
        // Оставляем один поток под UI, чтобы параллельные вычисления не блокировали интерфейс.
        QThreadPool::globalInstance()->setMaxThreadCount(
            qMax(1, QThread::idealThreadCount() - 1));

        connect(planetComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
                [this](int) {
            cancelTemperatureCalculation();
            updatePlanetSemiMajorAxisLabel();
            updatePlanetDayLengthLabel();
            updatePlanetYearLengthLabel();
            updatePlanetMassLabel();
            updatePlanetRadiusLabel();
            updatePlanetDerivedLabels();
            updateAtmospherePlanetParameters();
            updateAtmosphereComposition();
            updatePlanetOrbitLabels();
            updateSurfaceTemperatureAggregationHistoryDays();
            syncMaterialWithPlanet();
            syncGeothermalFluxWithPlanet();
            // При смене планеты пропускаем спин-ап, чтобы быстро показать оценку температуры.
            surfaceSkipSpinUp_ = true;
            if (isSurfaceTabActive()) {
                scheduleSurfaceGridTemperatureUpdate(true);
            } else {
                // Помечаем сетку, чтобы при открытии вкладки пересчитать
                // температурный профиль с учётом новой подповерхностной модели.
                surfaceGridDirty_ = true;
            }
            syncRotationModeWithPlanet();
            syncSpinOrbitWithPlanet();
            syncHeightSeedWithPlanet();
            syncHeightmapScaleWithPlanet();
            syncFlatHeightWithPlanet();
            syncAtmosphereDisableWithPlanet();
            syncCloudAlbedoWithPlanet();
            syncDiurnalCoolingBiasWithPlanet();
            syncManualGreenhouseOnTopWithPlanet();
            syncAdvancedRadiationWithPlanet();
            applyStellarPresetForCurrentPlanet();
            refreshSubsurfaceSettingsUi();
            updatePlanetActions();
            if (autoCalculateEnabled_ && hasPrimaryInputs() &&
                (!secondStarCheckBox_->isChecked() || hasSecondaryInputs())) {
                onCalculateRequested();
            } else {
                updateTemperaturePlot();
            }
            updateStarSystemOrbits();
        });

        connect(addPlanetButton_, &QPushButton::clicked, this, [this]() { onAddPlanetRequested(); });
        connect(deletePlanetButton_, &QPushButton::clicked, this, [this]() {
            const int index = planetComboBox_->currentIndex();
            if (index < 0 || !isCustomPlanetIndex(index)) {
                return;
            }

            const QString planetName = planetComboBox_->itemData(index, kRolePlanetName).toString();
            const auto result = QMessageBox::question(
                this,
                QStringLiteral("Удаление планеты"),
                QStringLiteral("Удалить планету \"%1\"?").arg(planetName));
            if (result != QMessageBox::Yes) {
                return;
            }
            planetComboBox_->removeItem(index);
            updatePlanetSemiMajorAxisLabel();
            updatePlanetDayLengthLabel();
            updatePlanetYearLengthLabel();
            updatePlanetMassLabel();
            updatePlanetRadiusLabel();
            updatePlanetDerivedLabels();
            updateAtmospherePlanetParameters();
            updateAtmosphereComposition();
            updateSurfaceTemperatureAggregationHistoryDays();
            updatePlanetActions();
            updateTemperaturePlot();
            refreshStarSystemView();
        });

        connect(materialComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int) {
            syncPlanetMaterialWithSelection();
            updateTemperaturePlot();
        });

        connect(rotationModeComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
                [this](int) {
            syncPlanetRotationModeWithSelection();
            updateRotationModeIllustration();
            updateTemperaturePlot();
        });

        connect(spinOrbitPSpinBox_, QOverload<int>::of(&QSpinBox::valueChanged), this,
                [this](int) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetSpinOrbitWithSelection();
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        });

        connect(spinOrbitQSpinBox_, QOverload<int>::of(&QSpinBox::valueChanged), this,
                [this](int) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetSpinOrbitWithSelection();
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        });

        connect(heightSeedSpinBox_, QOverload<int>::of(&QSpinBox::valueChanged), this, [this](int) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetHeightSeedWithSelection();
            updateSurfaceGridTemperatures();
        });

        connect(heightmapScaleSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double) {
                    if (planetComboBox_->currentIndex() < 0) {
                        return;
                    }
                    syncPlanetHeightmapScaleWithSelection();
                    updateSurfaceGridTemperatures();
                });

        connect(heightmapImportButton_, &QPushButton::clicked, this,
                [this]() { onImportHeightmapRequested(); });

        connect(flatHeightButton_, &QPushButton::toggled, this, [this](bool checked) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetFlatHeightWithSelection(checked);
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        });

        connect(disableAtmosphereButton_, &QPushButton::toggled, this, [this](bool checked) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetAtmosphereDisabledWithSelection(checked);
            updateAtmosphereComposition();
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        });

        connect(cloudAlbedoSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
                [this](double) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetCloudAlbedoWithSelection();
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        });

        connect(diurnalCoolingBiasSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetDiurnalCoolingBiasWithSelection();
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        });

        connect(manualGreenhouseOnTopCheckBox_, &QCheckBox::toggled, this, [this](bool) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetManualGreenhouseOnTopWithSelection();
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        });

        connect(advancedRadiationCheckBox_, &QCheckBox::toggled, this, [this](bool) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetAdvancedRadiationWithSelection();
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        });

        connect(debugLogCheckBox_, &QCheckBox::toggled, this, [this](bool checked) {
            debugLogEnabled_ = checked;
            setSolarDebugLoggingEnabled(checked);
        });

        connect(temperaturePauseButton_, &QPushButton::clicked, this, [this]() {
            const bool shouldPause = !temperaturePauseFlag_.load();
            temperaturePauseFlag_.store(shouldPause);
            updateTemperaturePauseUi(shouldPause);
            if (surfaceSimTimer_) {
                if (shouldPause) {
                    surfaceSimTimer_->stop();
                } else if (surfaceSimRunning_) {
                    surfaceSimTimer_->start();
                }
            }
        });

        connect(surfaceSimToggleButton_, &QPushButton::clicked, this,
                [this]() { toggleSurfaceSimulation(); });

        connect(surfaceSimSpeedComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
                [this](int) {
                    const double speedValue =
                        surfaceSimSpeedComboBox_->currentData().toDouble();
                    if (speedValue > 0.0) {
                        surfaceSimSpeedMultiplier_ = speedValue;
                        surfaceSimUnlimited_ = false;
                    } else {
                        surfaceSimUnlimited_ = true;
                    }
                    updateSurfaceSimulationTimerInterval();
                    updateSurfaceSimulationUi();
                });

        connect(surfaceMapModeComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
                [this](int) {
                    const SurfaceMapMode mode =
                        static_cast<SurfaceMapMode>(surfaceMapModeComboBox_->currentData().toInt());
                    applySurfaceMapMode(mode);
                });

        const auto onSubsurfaceChanged = [this]() {
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        };
        connect(subsurfaceLayersSpinBox_, QOverload<int>::of(&QSpinBox::valueChanged), this,
                [this, onSubsurfaceChanged](int) {
                    subsurfaceLayerCountAuto_ = false;
                    onSubsurfaceChanged();
                });
        connect(subsurfaceTopThicknessSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this, onSubsurfaceChanged](double) {
                    subsurfaceTopThicknessAuto_ = false;
                    onSubsurfaceChanged();
                });
        connect(subsurfaceDepthSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
                [this, onSubsurfaceChanged](double) {
                    subsurfaceBottomDepthAuto_ = false;
                    onSubsurfaceChanged();
                });
        connect(subsurfaceGeothermalFluxSpinBox_,
                QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this,
                [onSubsurfaceChanged](double) { onSubsurfaceChanged(); });
        connect(subsurfaceBoundaryComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
                [onSubsurfaceChanged](int) { onSubsurfaceChanged(); });

        connect(surfaceViewToggleButton_, &QPushButton::toggled, this, [this](bool checked) {
            if (!surfaceViewStack_) {
                return;
            }
            surfaceViewStack_->setCurrentWidget(checked
                                                    ? static_cast<QWidget *>(surfaceGlobeWidget_)
                                                    : static_cast<QWidget *>(surfaceMapWidget_));
            surfaceViewToggleButton_->setText(checked ? QStringLiteral("2D вид")
                                                      : QStringLiteral("3D вид"));
        });
        surfaceViewToggleButton_->setChecked(true);

        // Сглаживание влияет только на отображение кривых, а не на физический расчет.
        // connect(smoothingCheckBox, &QCheckBox::toggled, temperaturePlot_,
        //         &SurfaceTemperaturePlot::setSmoothingEnabled);

        connect(calculateButton, &QPushButton::clicked, this, [this]() {
            autoCalculateEnabled_ = true;
            onCalculateRequested();
        });

        applyPrimary(StellarParameters{1.0, 5772.0, 1.0});
        applySecondary(std::nullopt);
        setPlanetPresets(solarSystemPresets(), QStringLiteral("Венера"));

        starSystemTimer_ = new QTimer(this);
        starSystemTimer_->setInterval(50);
        connect(starSystemTimer_, &QTimer::timeout, this, [this]() { advanceStarSystemTime(); });
        updateStarSystemTimerState();
    }

    ~SolarCalculatorWidget() override {
        cancelSurfaceGridTemperatureUpdate();
        if (surfaceGridWatcher_) {
            surfaceGridWatcher_->future().waitForFinished();
        }
    }

private:
    struct SurfaceFluxBreakdown {
        double totalFluxWPerM2 = 0.0;
        double primaryFluxWPerM2 = 0.0;
        double secondaryFluxWPerM2 = 0.0;
        double primaryDistanceAu = 0.0;
        double secondaryDistanceAu = 0.0;
        bool hasPrimary = false;
        bool hasSecondary = false;
    };

    void onCalculateRequested() {
        BinarySystemParameters parameters{};

        if (!readStellarParameters(radiusInput_, temperatureInput_,
                                   QStringLiteral("первой звезды"), parameters.primary)) {
            return;
        }

        if (secondStarCheckBox_->isChecked()) {
            StellarParameters secondary{};
            if (!readStellarParameters(secondaryRadiusInput_, secondaryTemperatureInput_,
                                       QStringLiteral("второй звезды"), secondary)) {
                return;
            }
            parameters.secondary = secondary;
        }

        const StellarCacheKey stellarKey{
            parameters.primary.radiusInSolarRadii,
            parameters.primary.temperatureKelvin,
            static_cast<bool>(parameters.secondary),
            parameters.secondary ? parameters.secondary->radiusInSolarRadii : 0.0,
            parameters.secondary ? parameters.secondary->temperatureKelvin : 0.0};
        if (!lastStellarKey_ || *lastStellarKey_ != stellarKey) {
            lastStellarKey_ = stellarKey;
        }

        double semiMajorAxis = 0.0;
        if (!readSemiMajorAxis(semiMajorAxis)) {
            return;
        }

        parameters.primary.distanceInAU = semiMajorAxis;
        if (parameters.secondary) {
            parameters.secondary->distanceInAU = semiMajorAxis;
        }

        const double primaryFlux = SolarCalculator::solarConstant(parameters.primary);
        double totalFlux = primaryFlux;
        QString details;

        if (parameters.secondary) {
            const double secondaryFlux = SolarCalculator::solarConstant(*parameters.secondary);
            totalFlux += secondaryFlux;
            details = QStringLiteral(" (первая: %1 Вт/м², вторая: %2 Вт/м²)")
                          .arg(primaryFlux, 0, 'g', precision_)
                          .arg(secondaryFlux, 0, 'g', precision_);
        }

        resultLabel_->setText(
            QStringLiteral("Солнечная постоянная у планеты: %1 Вт/м²%2")
                .arg(totalFlux, 0, 'g', precision_)
                .arg(details));

        lastSolarConstant_ = totalFlux;
        lastSolarConstantDistanceAU_ = semiMajorAxis;
        hasSolarConstant_ = true;
        refreshSubsurfaceSettingsUi();
        updateTemperaturePlot();
        updateSurfaceOrbitLighting();
    }

    bool readStellarParameters(QLineEdit *radiusInput, QLineEdit *temperatureInput,
                               const QString &label, StellarParameters &parameters) {
        bool ok = false;
        const double radius = radiusInput->text().toDouble(&ok);
        if (!ok || radius <= 0.0) {
            showInputError(QStringLiteral("Укажите положительный радиус %1.").arg(label));
            return false;
        }

        const double temperature = temperatureInput->text().toDouble(&ok);
        if (!ok || temperature <= 0.0) {
            showInputError(QStringLiteral("Укажите положительную температуру %1.").arg(label));
            return false;
        }

        parameters.radiusInSolarRadii = radius;
        parameters.temperatureKelvin = temperature;
        return true;
    }

    bool readSemiMajorAxis(double &semiMajorAxis) {
        const QVariant value = planetComboBox_->currentData(kRoleSemiMajorAxis);
        if (!value.isValid()) {
            showInputError(QStringLiteral("Выберите планету из списка."));
            return false;
        }
        semiMajorAxis = value.toDouble();
        if (semiMajorAxis <= 0.0) {
            showInputError(QStringLiteral("Укажите положительную большую полуось орбиты планеты."));
            return false;
        }
        return true;
    }

    void showInputError(const QString &message) {
        QMessageBox::warning(this, QStringLiteral("Некорректный ввод"), message);
    }

    bool shouldLogRadiationForPoint(int pointIndex) const {
        if (!debugLogEnabled_) {
            return false;
        }
        if (surfaceGrid_.points().isEmpty()) {
            return false;
        }
        int targetIndex = selectedSurfacePointIndex_;
        if (targetIndex < 0 || targetIndex >= surfaceGrid_.points().size()) {
            // Если ячейка не выбрана, логируем нулевую — так лог остаётся однострочным и предсказуемым.
            targetIndex = 0;
        }
        return pointIndex == targetIndex;
    }

    void ensureTemperaturePlotSourceSurfaceGrid() {
        if (!temperaturePlotSourceComboBox_) {
            return;
        }
        const int surfaceIndex = temperaturePlotSourceComboBox_->findData(
            static_cast<int>(TemperaturePlotSource::SurfaceGrid));
        if (surfaceIndex < 0 || temperaturePlotSourceComboBox_->currentIndex() == surfaceIndex) {
            return;
        }
        const QSignalBlocker blocker(temperaturePlotSourceComboBox_);
        temperaturePlotSourceComboBox_->setCurrentIndex(surfaceIndex);
    }

    QLineEdit *radiusInput_ = nullptr;
    QLineEdit *temperatureInput_ = nullptr;

    QCheckBox *secondStarCheckBox_ = nullptr;
    QLineEdit *secondaryRadiusInput_ = nullptr;
    QLineEdit *secondaryTemperatureInput_ = nullptr;

    QGroupBox *primaryGroupBox_ = nullptr;
    QGroupBox *secondaryGroupBox_ = nullptr;
    QComboBox *planetComboBox_ = nullptr;
    QLabel *planetSemiMajorAxisLabel_ = nullptr;
    QLabel *planetDayLengthLabel_ = nullptr;
    QLabel *planetYearLengthLabel_ = nullptr;
    QLabel *planetMassLabel_ = nullptr;
    QLabel *planetRadiusLabel_ = nullptr;
    QLabel *planetSurfaceGravityLabel_ = nullptr;
    QLabel *planetSurfaceAreaLabel_ = nullptr;
    QLabel *planetEccentricityLabel_ = nullptr;
    QLabel *planetObliquityLabel_ = nullptr;
    QLabel *planetPerihelionArgumentLabel_ = nullptr;
    QComboBox *materialComboBox_ = nullptr;
    QComboBox *rotationModeComboBox_ = nullptr;
    QSpinBox *spinOrbitPSpinBox_ = nullptr;
    QSpinBox *spinOrbitQSpinBox_ = nullptr;
    QSpinBox *heightSeedSpinBox_ = nullptr;
    QPushButton *heightmapImportButton_ = nullptr;
    QDoubleSpinBox *heightmapScaleSpinBox_ = nullptr;
    QPushButton *flatHeightButton_ = nullptr;
    QPushButton *disableAtmosphereButton_ = nullptr;
    QDoubleSpinBox *cloudAlbedoSpinBox_ = nullptr;
    QDoubleSpinBox *diurnalCoolingBiasSpinBox_ = nullptr;
    QCheckBox *manualGreenhouseOnTopCheckBox_ = nullptr;
    QCheckBox *advancedRadiationCheckBox_ = nullptr;
    QCheckBox *debugLogCheckBox_ = nullptr;
    ModeIllustrationWidget *modeIllustrationWidget_ = nullptr;
    QPushButton *addPlanetButton_ = nullptr;
    QPushButton *deletePlanetButton_ = nullptr;
    AtmosphereWidget *atmosphereWidget_ = nullptr;
    QTabWidget *rightTabs_ = nullptr;
    QWidget *surfaceMapContainer_ = nullptr;

    QLabel *resultLabel_ = nullptr;
    SurfaceTemperaturePlot *temperaturePlot_ = nullptr;
    StarSystemTopView *starSystemTopView_ = nullptr;
    SurfaceMapWidget *surfaceMapWidget_ = nullptr;
    SurfaceGlobeWidget *surfaceGlobeWidget_ = nullptr;
    QPointer<SurfacePointStatusDialog> surfacePointStatusDialog_;
    int selectedSurfacePointIndex_ = -1;
    QStackedWidget *surfaceViewStack_ = nullptr;
    QLabel *surfaceMinTemperatureLabel_ = nullptr;
    QLabel *surfaceMaxTemperatureLabel_ = nullptr;
    QLabel *surfaceSolarConstantLabel_ = nullptr;
    QPushButton *temperaturePauseButton_ = nullptr;
    QPushButton *surfaceSimToggleButton_ = nullptr;
    QComboBox *surfaceSimSpeedComboBox_ = nullptr;
    QLabel *surfaceSimTimeLabel_ = nullptr;
    QComboBox *surfaceMapModeComboBox_ = nullptr;
    QPushButton *surfaceViewToggleButton_ = nullptr;
    QCheckBox *surfaceMarkupCheckBox_ = nullptr;
    bool debugLogEnabled_ = false;
    QLabel *surfaceCalculationLabel_ = nullptr;
    QProgressBar *surfaceCalculationProgress_ = nullptr;
    QSpinBox *subsurfaceLayersSpinBox_ = nullptr;
    QDoubleSpinBox *subsurfaceTopThicknessSpinBox_ = nullptr;
    QDoubleSpinBox *subsurfaceDepthSpinBox_ = nullptr;
    QComboBox *subsurfaceBoundaryComboBox_ = nullptr;
    QDoubleSpinBox *subsurfaceGeothermalFluxSpinBox_ = nullptr;
    bool subsurfaceLayerCountAuto_ = true;
    bool subsurfaceTopThicknessAuto_ = true;
    bool subsurfaceBottomDepthAuto_ = true;
    SurfaceTemperatureScaleWidget *temperatureScaleWidget_ = nullptr;
    SurfaceHeightScaleWidget *heightScaleWidget_ = nullptr;
    SurfaceWindScaleWidget *windScaleWidget_ = nullptr;
    SurfacePressureScaleWidget *pressureScaleWidget_ = nullptr;
    SurfaceRealisticScaleWidget *realisticScaleWidget_ = nullptr;
    QStackedWidget *surfaceLegendScaleStack_ = nullptr;
    SegmentSelectorWidget *segmentSelectorWidget_ = nullptr;
    QComboBox *temperaturePlotSourceComboBox_ = nullptr;
    QProgressDialog *temperatureProgressDialog_ = nullptr;
    QFutureWatcher<SurfaceGridComputationResult> *surfaceGridWatcher_ = nullptr;
    QElapsedTimer temperatureElapsed_;
    QTimer *temperatureUiTimer_ = nullptr;
    QTimer *surfaceSimTimer_ = nullptr;
    QTimer *surfaceGridDebounceTimer_ = nullptr;
    QElapsedTimer starSystemElapsedTimer_;
    QTimer *starSystemTimer_ = nullptr;
    PlanetSurfaceGrid surfaceGrid_;
    std::shared_ptr<std::atomic_bool> temperatureCancelFlag_;
    std::shared_ptr<std::atomic_bool> surfaceGridCancelFlag_;
    std::atomic_bool temperaturePauseFlag_{false};
    int temperatureRequestId_ = 0;
    int precision_ = kDefaultPrecision;
    std::shared_ptr<std::atomic<int>> surfaceGridRequestId_ =
        std::make_shared<std::atomic<int>>(0);
    bool surfaceGridCalculationRunning_ = false;
    QSet<QString> presetPlanetNames_;
    double lastSolarConstant_ = 0.0;
    double lastSolarConstantDistanceAU_ = 0.0;
    bool hasSolarConstant_ = false;
    double surfaceMinTemperatureK_ = 0.0;
    double surfaceMaxTemperatureK_ = 0.0;
    bool hasSurfaceTemperatureRange_ = false;
    double surfaceMinAirTemperatureK_ = 0.0;
    double surfaceMaxAirTemperatureK_ = 0.0;
    bool hasSurfaceAirTemperatureRange_ = false;
    double surfaceMinHeightKm_ = 0.0;
    double surfaceMaxHeightKm_ = 0.0;
    bool hasSurfaceHeightRange_ = false;
    double surfaceMinWindSpeedMps_ = 0.0;
    double surfaceMaxWindSpeedMps_ = 0.0;
    bool hasSurfaceWindRange_ = false;
    double surfaceMinPressureAtm_ = 0.0;
    double surfaceMaxPressureAtm_ = 0.0;
    bool hasSurfacePressureRange_ = false;
    SurfaceMapMode surfaceMapMode_ = SurfaceMapMode::Temperature;
    bool autoCalculateEnabled_ = false;
    double surfaceSimSpeedMultiplier_ = 1.0;
    bool surfaceSimUnlimited_ = false;
    std::optional<StellarCacheKey> lastStellarKey_;
    bool surfaceSimRunning_ = false;
    bool surfaceSimResumeAfterGridUpdate_ = false;
    bool surfaceGridDirty_ = false;
    bool surfaceSkipSpinUp_ = false;
    struct SurfaceTimeFlow {
        int dayIndex = 0;
        int hourIndex = 0;
        // Орбитальное время ведём в земных часах (24 ч/сутки), отдельно от планетных суток.
        qint64 orbitHourIndex = 0;
        double orbitSegmentAccumulator = 0.0;
        int stepsPerDay = 24;
        double orbitSegmentStep = 0.0;

        void reset() {
            dayIndex = 0;
            hourIndex = 0;
            orbitHourIndex = 0;
            orbitSegmentAccumulator = 0.0;
            stepsPerDay = 24;
            orbitSegmentStep = 0.0;
        }

        int surfaceTotalHourIndex() const {
            return (stepsPerDay > 0) ? dayIndex * stepsPerDay + hourIndex : 0;
        }

        qint64 orbitTotalHourIndex() const {
            return orbitHourIndex;
        }

        double surfaceElapsedDays() const {
            if (stepsPerDay <= 0) {
                return 0.0;
            }
            return static_cast<double>(dayIndex) +
                static_cast<double>(hourIndex) / static_cast<double>(stepsPerDay);
        }

        double orbitElapsedDays() const {
            return static_cast<double>(orbitHourIndex) / 24.0;
        }

        void configure(double dayLengthDays, double yearLengthDays) {
            const int previousSteps = stepsPerDay;
            const int totalHours =
                (previousSteps > 0) ? dayIndex * previousSteps + hourIndex : 0;
            stepsPerDay = qMax(1, qRound(dayLengthDays * 24.0));
            if (stepsPerDay > 0) {
                dayIndex = totalHours / stepsPerDay;
                hourIndex = totalHours % stepsPerDay;
            } else {
                dayIndex = 0;
                hourIndex = 0;
            }
            const double resolvedYearLengthDays =
                (yearLengthDays > 0.0)
                    ? yearLengthDays
                    : static_cast<double>(kSurfaceOrbitSegmentsPerYear);
            const double yearLengthHours = qMax(1.0, resolvedYearLengthDays * 24.0);
            // Один час даёт долю орбиты: (segments / orbit) / (hours / orbit).
            orbitSegmentStep =
                static_cast<double>(kSurfaceOrbitSegmentsPerYear) / yearLengthHours;
            const double exactSegments =
                static_cast<double>(orbitHourIndex) * orbitSegmentStep;
            orbitSegmentAccumulator = exactSegments - std::floor(exactSegments);
        }

        int resolvedOrbitSegmentIndex() const {
            if (orbitSegmentStep <= 0.0) {
                return 0;
            }
            const double exactSegments =
                static_cast<double>(orbitHourIndex) * orbitSegmentStep;
            return static_cast<int>(std::floor(exactSegments));
        }

        void advanceHour(OrbitAnimationModel *orbitAnimation) {
            ++orbitHourIndex;
            orbitSegmentAccumulator += orbitSegmentStep;
            const int segmentsToAdvance =
                static_cast<int>(std::floor(orbitSegmentAccumulator));
            if (orbitAnimation) {
                for (int i = 0; i < segmentsToAdvance; ++i) {
                    orbitAnimation->advanceSegment();
                }
            }
            orbitSegmentAccumulator -= segmentsToAdvance;
            ++hourIndex;
            if (hourIndex >= stepsPerDay) {
                hourIndex = 0;
                ++dayIndex;
            }
        }
    } surfaceTime_;
    struct SurfaceTemperatureAggregationState {
        double targetLongitudeRadians = 0.0;
        bool hasTargetLongitude = false;
        QVector<int> pointIndicesByLatitude;
        QVector<double> dailyMinimums;
        QVector<double> dailyMaximums;
        QVector<QVector<double>> minimumHistory;
        QVector<QVector<double>> maximumHistory;
        int historyDays = kSurfaceTemperatureHistoryDays;
    } surfaceTemperatureAggregation_;
    OrbitAnimationModel surfaceOrbitAnimation_;
    bool surfaceOrbitAnimationInitialized_ = false;
    double starSystemElapsedDays_ = 0.0;
    struct SurfaceOrbitConfig {
        double semiMajorAxis = 1.0;
        double eccentricity = 0.0;
        double obliquity = 0.0;
        double perihelionArgument = 0.0;
        bool isValid = false;
    } surfaceOrbitConfig_;

    bool isStarSystemTabActive() const {
        return rightTabs_ && starSystemTopView_ &&
            rightTabs_->currentWidget() == starSystemTopView_;
    }

    void updateStarSystemTimerState() {
        if (isStarSystemTabActive()) {
            startStarSystemTimer();
            updateStarSystemOrbits();
        } else {
            stopStarSystemTimer();
        }
    }

    void startStarSystemTimer() {
        if (!starSystemTimer_ || starSystemTimer_->isActive()) {
            return;
        }
        starSystemElapsedTimer_.start();
        starSystemTimer_->start();
    }

    void stopStarSystemTimer() {
        if (!starSystemTimer_ || !starSystemTimer_->isActive()) {
            return;
        }
        starSystemTimer_->stop();
    }

    void advanceStarSystemTime() {
        if (!starSystemTimer_) {
            return;
        }
        if (!starSystemElapsedTimer_.isValid()) {
            starSystemElapsedTimer_.start();
            return;
        }
        const qint64 elapsedMs = starSystemElapsedTimer_.restart();
        const double elapsedSeconds = static_cast<double>(elapsedMs) / 1000.0;
        starSystemElapsedDays_ += elapsedSeconds * kStarSystemDaysPerSecond;
        updateStarSystemOrbits();
    }

    void updateFlatHeightButtonText(bool useFlatHeight) {
        if (!flatHeightButton_) {
            return;
        }
        flatHeightButton_->setText(useFlatHeight ? QStringLiteral("Вернуть высоту")
                                                 : QStringLiteral("Обнулить высоту"));
    }

    void updateAtmosphereButtonText(bool atmosphereDisabled) {
        if (!disableAtmosphereButton_) {
            return;
        }
        disableAtmosphereButton_->setText(atmosphereDisabled
                                              ? QStringLiteral("Вернуть атмосферу")
                                              : QStringLiteral("Отключить атмосферу"));
    }

    void updateTemperaturePauseUi(bool paused) {
        if (temperaturePauseButton_) {
            temperaturePauseButton_->setText(paused ? QStringLiteral("Продолжить")
                                                    : QStringLiteral("Пауза"));
        }
        if (temperatureProgressDialog_) {
            auto *pauseButton = temperatureProgressDialog_->findChild<QPushButton *>(
                QStringLiteral("temperaturePauseButton"));
            if (pauseButton) {
                pauseButton->setText(paused ? QStringLiteral("Продолжить")
                                            : QStringLiteral("Пауза"));
            }
        }
    }

    void updateSurfaceSimulationUi() {
        if (surfaceSimToggleButton_) {
            surfaceSimToggleButton_->setText(
                surfaceSimRunning_
                    ? QStringLiteral("Остановить течение")
                    : QStringLiteral("Течение времени"));
        }
        if (surfaceSimTimeLabel_) {
            const QString speedLabel = surfaceSimUnlimited_
                                           ? QStringLiteral("∞")
                                           : QStringLiteral("%1x")
                                                 .arg(surfaceSimSpeedMultiplier_, 0, 'g', 3);
            if (!surfaceSimRunning_ && surfaceTime_.dayIndex == 0 &&
                surfaceTime_.hourIndex == 0) {
                surfaceSimTimeLabel_->setText(
                    QStringLiteral("t = — (%1)").arg(speedLabel));
            } else {
                surfaceSimTimeLabel_->setText(
                    QStringLiteral("t = День %1, Час %2 (%3)")
                        .arg(surfaceTime_.dayIndex + 1)
                        .arg(surfaceTime_.hourIndex + 1)
                        .arg(speedLabel));
            }
        }
    }

    void updateSurfaceSimulationTimerInterval() {
        if (!surfaceSimTimer_) {
            return;
        }
        // Интервал 0 мс в QTimer означает максимально возможную частоту событий.
        const int intervalMs = surfaceSimUnlimited_
                                   ? 0
                                   : qMax(1, qRound(1000.0 / surfaceSimSpeedMultiplier_));
        surfaceSimTimer_->setInterval(intervalMs);
    }

    void syncSurfaceTimeSettingsFromPlanet() {
        if (!planetComboBox_) {
            surfaceTime_.configure(1.0, static_cast<double>(kSurfaceOrbitSegmentsPerYear));
            return;
        }
        const double dayLengthDays = planetComboBox_->currentData(kRoleDayLength).toDouble();
        const double yearLengthDays =
            planetComboBox_->currentData(kRoleYearLengthDays).toDouble();
        surfaceTime_.configure(dayLengthDays, yearLengthDays);
    }

    int resolveSurfaceStepsPerDay() {
        syncSurfaceTimeSettingsFromPlanet();
        return surfaceTime_.stepsPerDay;
    }

    double resolveSurfaceElapsedDays() {
        syncSurfaceTimeSettingsFromPlanet();
        return surfaceTime_.surfaceElapsedDays();
    }

    OrbitSegment resolvePlanetOrbitSegment(double elapsedDays,
                                           double yearLengthDays,
                                           double semiMajorAxisAu,
                                           double eccentricity) const {
        OrbitSegment segment;
        segment.distanceAU = semiMajorAxisAu;
        if (yearLengthDays <= 0.0 || semiMajorAxisAu <= 0.0) {
            return segment;
        }
        double fraction = std::fmod(elapsedDays / yearLengthDays, 1.0);
        if (fraction < 0.0) {
            fraction += 1.0;
        }
        const double meanAnomaly = 2.0 * M_PI * fraction;
        OrbitSegmentCalculator calculator(semiMajorAxisAu, eccentricity);
        return calculator.orbitAtMeanAnomaly(meanAnomaly);
    }

    void resetSurfaceSimulation(bool resetTime) {
        surfaceSimRunning_ = false;
        if (resetTime) {
            surfaceTime_.reset();
            starSystemElapsedDays_ = 0.0;
        }
        resetSurfaceTemperatureAggregation();
        resetSurfaceOrbitAnimation();
        if (surfaceSimTimer_) {
            surfaceSimTimer_->stop();
        }
        updateSurfaceSimulationUi();
        updateStarSystemOrbits();
    }

    void resumeSurfaceSimulationAfterGridUpdate() {
        if (!surfaceSimResumeAfterGridUpdate_) {
            return;
        }
        surfaceSimResumeAfterGridUpdate_ = false;
        if (!surfaceSimTimer_) {
            surfaceSimTimer_ = new QTimer(this);
            connect(surfaceSimTimer_, &QTimer::timeout, this,
                    [this]() { advanceSurfaceSimulation(); });
        }
        updateSurfaceSimulationTimerInterval();
        surfaceSimRunning_ = true;
        updateSurfaceSimulationUi();
        if (!temperaturePauseFlag_.load()) {
            surfaceSimTimer_->start();
        }
        updateStarSystemOrbits();
    }

    void advanceSurfaceOrbitAnimationToCurrentTime() {
        if (!surfaceOrbitAnimationInitialized_) {
            return;
        }
        const int segmentsPerYear = qMax(1, kSurfaceOrbitSegmentsPerYear);
        const int segmentsToAdvance =
            surfaceTime_.resolvedOrbitSegmentIndex() % segmentsPerYear;
        for (int i = 0; i < segmentsToAdvance; ++i) {
            surfaceOrbitAnimation_.advanceSegment();
        }
    }

    void resetSurfaceOrbitAnimation() {
        if (!planetComboBox_) {
            surfaceOrbitAnimation_.reset(1.0, 0.0, 0.0, 0.0, kSurfaceOrbitSegmentsPerYear);
            surfaceOrbitConfig_ = {};
            surfaceOrbitConfig_.isValid = true;
            surfaceOrbitAnimationInitialized_ = true;
            return;
        }
        const double semiMajorAxis = planetComboBox_->currentData(kRoleSemiMajorAxis).toDouble();
        const double eccentricity = planetComboBox_->currentData(kRoleEccentricity).toDouble();
        const double obliquity = planetComboBox_->currentData(kRoleObliquity).toDouble();
        const double perihelionArgument =
            planetComboBox_->currentData(kRolePerihelionArgument).toDouble();
        surfaceOrbitAnimation_.reset(semiMajorAxis,
                                     eccentricity,
                                     obliquity,
                                     perihelionArgument,
                                     kSurfaceOrbitSegmentsPerYear);
        surfaceOrbitConfig_.semiMajorAxis = semiMajorAxis;
        surfaceOrbitConfig_.eccentricity = eccentricity;
        surfaceOrbitConfig_.obliquity = obliquity;
        surfaceOrbitConfig_.perihelionArgument = perihelionArgument;
        surfaceOrbitConfig_.isValid = true;
        surfaceOrbitAnimationInitialized_ = true;
    }

    void syncSurfaceOrbitAnimationToTime() {
        syncSurfaceTimeSettingsFromPlanet();
        if (!planetComboBox_) {
            if (!surfaceOrbitAnimationInitialized_) {
                resetSurfaceOrbitAnimation();
            }
            return;
        }
        const double semiMajorAxis = planetComboBox_->currentData(kRoleSemiMajorAxis).toDouble();
        const double eccentricity = planetComboBox_->currentData(kRoleEccentricity).toDouble();
        const double obliquity = planetComboBox_->currentData(kRoleObliquity).toDouble();
        const double perihelionArgument =
            planetComboBox_->currentData(kRolePerihelionArgument).toDouble();
        const bool configMatches =
            surfaceOrbitConfig_.isValid &&
            qFuzzyCompare(surfaceOrbitConfig_.semiMajorAxis + 1.0, semiMajorAxis + 1.0) &&
            qFuzzyCompare(surfaceOrbitConfig_.eccentricity + 1.0, eccentricity + 1.0) &&
            qFuzzyCompare(surfaceOrbitConfig_.obliquity + 1.0, obliquity + 1.0) &&
            qFuzzyCompare(surfaceOrbitConfig_.perihelionArgument + 1.0,
                          perihelionArgument + 1.0);
        if (!configMatches) {
            resetSurfaceOrbitAnimation();
            advanceSurfaceOrbitAnimationToCurrentTime();
            return;
        }
        if (!surfaceOrbitAnimationInitialized_) {
            resetSurfaceOrbitAnimation();
            advanceSurfaceOrbitAnimationToCurrentTime();
        }
    }

    void ensureSurfaceOrbitAnimationReady() {
        syncSurfaceOrbitAnimationToTime();
        if (!surfaceOrbitAnimationInitialized_) {
            resetSurfaceOrbitAnimation();
        }
    }

    void advanceSurfaceTimeFlow() {
        syncSurfaceOrbitAnimationToTime();
        surfaceTime_.advanceHour(&surfaceOrbitAnimation_);
        updateSurfaceSimulationUi();
        updateStarSystemOrbits();
    }

    void resetSurfaceTemperatureAggregation() {
        surfaceTemperatureAggregation_ = {};
        updateSurfaceTemperatureAggregationHistoryDays();
        if (temperaturePlot_) {
            temperaturePlot_->clearSeries();
        }
    }

    int resolveSurfaceTemperatureHistoryDays() const {
        if (!planetComboBox_) {
            return kSurfaceTemperatureHistoryDays;
        }
        const QVariant value = planetComboBox_->currentData(kRoleYearLengthDays);
        if (!value.isValid()) {
            return kSurfaceTemperatureHistoryDays;
        }
        const double yearLengthDays = value.toDouble();
        if (yearLengthDays <= 0.0) {
            return kSurfaceTemperatureHistoryDays;
        }
        return qMax(1, static_cast<int>(qCeil(yearLengthDays)));
    }

    void updateSurfaceTemperatureAggregationHistoryDays() {
        const int historyDays = resolveSurfaceTemperatureHistoryDays();
        surfaceTemperatureAggregation_.historyDays = historyDays;
        const int historySize = surfaceTemperatureAggregation_.minimumHistory.size();
        for (int i = 0; i < historySize; ++i) {
            auto &minHistory = surfaceTemperatureAggregation_.minimumHistory[i];
            auto &maxHistory = surfaceTemperatureAggregation_.maximumHistory[i];
            if (minHistory.size() > historyDays) {
                minHistory.remove(0, minHistory.size() - historyDays);
            }
            if (maxHistory.size() > historyDays) {
                maxHistory.remove(0, maxHistory.size() - historyDays);
            }
        }
    }

    double selectSurfaceAggregationLongitude(double substellarLongitudeRadians) const {
        if (selectedSurfacePointIndex_ >= 0 &&
            selectedSurfacePointIndex_ < surfaceGrid_.points().size()) {
            // Если пользователь выбрал конкретную точку, фиксируем меридиан по ней,
            // чтобы график отражал знакомое положение на карте.
            return surfaceGrid_.points().at(selectedSurfacePointIndex_).longitudeRadians;
        }
        // Иначе берём подсолнечный меридиан: он физически осмысленный и не зависит от шумов сетки.
        return substellarLongitudeRadians;
    }

    void ensureSurfaceTemperatureAggregationTargets(double targetLongitudeRadians) {
        const int binCount = latitudePoints();
        if (!surfaceTemperatureAggregation_.hasTargetLongitude ||
            surfaceTemperatureAggregation_.pointIndicesByLatitude.size() != binCount ||
            !qFuzzyCompare(surfaceTemperatureAggregation_.targetLongitudeRadians + 1.0,
                           targetLongitudeRadians + 1.0)) {
            surfaceTemperatureAggregation_.targetLongitudeRadians = targetLongitudeRadians;
            surfaceTemperatureAggregation_.hasTargetLongitude = true;
            surfaceTemperatureAggregation_.pointIndicesByLatitude.resize(binCount);
            surfaceTemperatureAggregation_.dailyMinimums.resize(binCount);
            surfaceTemperatureAggregation_.dailyMaximums.resize(binCount);
            surfaceTemperatureAggregation_.minimumHistory.resize(binCount);
            surfaceTemperatureAggregation_.maximumHistory.resize(binCount);
            surfaceTemperatureAggregation_.dailyMinimums.fill(0.0);
            surfaceTemperatureAggregation_.dailyMaximums.fill(0.0);
            for (int binIndex = 0; binIndex < binCount; ++binIndex) {
                surfaceTemperatureAggregation_.minimumHistory[binIndex].clear();
                surfaceTemperatureAggregation_.maximumHistory[binIndex].clear();
            }

            const double step = latitudeStepDegrees();
            for (int binIndex = 0; binIndex < binCount; ++binIndex) {
                const double binLatitude = -90.0 + step * static_cast<double>(binIndex);
                double bestLatDiff = std::numeric_limits<double>::max();
                double bestLonDiff = std::numeric_limits<double>::max();
                int bestIndex = -1;
                for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
                    const auto &point = surfaceGrid_.points().at(i);
                    const double latDiff = std::abs(point.latitudeDeg - binLatitude);
                    const double lonDiff =
                        longitudeDistanceRadians(point.longitudeRadians, targetLongitudeRadians);
                    // Фиксируем меридиан: сначала выбираем ближайшую широту,
                    // а при равенстве берём точку, лежащую ближе к нужной долготе.
                    if (latDiff < bestLatDiff - 1e-6 ||
                        (std::abs(latDiff - bestLatDiff) <= 1e-6 && lonDiff < bestLonDiff)) {
                        bestLatDiff = latDiff;
                        bestLonDiff = lonDiff;
                        bestIndex = i;
                    }
                }
                surfaceTemperatureAggregation_.pointIndicesByLatitude[binIndex] = bestIndex;
            }
        }
    }

    void updateSurfaceTemperatureAggregation(bool publishDaily,
                                             RotationMode rotationMode,
                                             bool hasAtmosphere) {
        if (!temperaturePlot_ || surfaceGrid_.points().isEmpty() ||
            surfaceTemperatureAggregation_.pointIndicesByLatitude.isEmpty()) {
            return;
        }

        const int binCount = surfaceTemperatureAggregation_.pointIndicesByLatitude.size();
        for (int binIndex = 0; binIndex < binCount; ++binIndex) {
            const int pointIndex =
                surfaceTemperatureAggregation_.pointIndicesByLatitude.at(binIndex);
            if (pointIndex < 0 || pointIndex >= surfaceGrid_.points().size()) {
                continue;
            }
            const double temperatureK = surfaceGrid_.points().at(pointIndex).temperatureK;
            double &minValue = surfaceTemperatureAggregation_.dailyMinimums[binIndex];
            double &maxValue = surfaceTemperatureAggregation_.dailyMaximums[binIndex];
            if (minValue == 0.0 && maxValue == 0.0) {
                minValue = temperatureK;
                maxValue = temperatureK;
            } else {
                minValue = qMin(minValue, temperatureK);
                maxValue = qMax(maxValue, temperatureK);
            }
        }

        if (!publishDaily) {
            return;
        }

        QVector<TemperatureRangePoint> dailyPoints;
        QVector<TemperatureSummaryPoint> summaryPoints;
        dailyPoints.reserve(binCount);
        summaryPoints.reserve(binCount);

        const double step = latitudeStepDegrees();
        const int historyLimit = qMax(1, surfaceTemperatureAggregation_.historyDays);
        for (int binIndex = 0; binIndex < binCount; ++binIndex) {
            const double latitude = -90.0 + step * static_cast<double>(binIndex);
            double dailyMin = surfaceTemperatureAggregation_.dailyMinimums.at(binIndex);
            double dailyMax = surfaceTemperatureAggregation_.dailyMaximums.at(binIndex);
            if (dailyMin == 0.0 && dailyMax == 0.0) {
                dailyMin = 0.0;
                dailyMax = 0.0;
            }

            auto &minHistory = surfaceTemperatureAggregation_.minimumHistory[binIndex];
            auto &maxHistory = surfaceTemperatureAggregation_.maximumHistory[binIndex];
            minHistory.push_back(dailyMin);
            maxHistory.push_back(dailyMax);
            if (minHistory.size() > historyLimit) {
                minHistory.remove(0, minHistory.size() - historyLimit);
            }
            if (maxHistory.size() > historyLimit) {
                maxHistory.remove(0, maxHistory.size() - historyLimit);
            }

            const double dailyMean = 0.5 * (dailyMin + dailyMax);
            TemperatureRangePoint dailyPoint;
            dailyPoint.latitudeDegrees = latitude;
            dailyPoint.hasInsolation = true;
            dailyPoint.minimumKelvin = dailyMin;
            dailyPoint.maximumKelvin = dailyMax;
            dailyPoint.meanDailyKelvin = dailyMean;
            dailyPoint.meanDayKelvin = dailyMax;
            dailyPoint.meanNightKelvin = dailyMin;
            dailyPoint.minimumCelsius = dailyMin - kKelvinOffset;
            dailyPoint.maximumCelsius = dailyMax - kKelvinOffset;
            dailyPoint.meanDailyCelsius = dailyMean - kKelvinOffset;
            dailyPoint.meanDayCelsius = dailyMax - kKelvinOffset;
            dailyPoint.meanNightCelsius = dailyMin - kKelvinOffset;
            dailyPoints.push_back(dailyPoint);

            double minAnnual = std::numeric_limits<double>::max();
            double maxAnnual = std::numeric_limits<double>::lowest();
            double meanAnnualSum = 0.0;
            double meanAnnualDaySum = 0.0;
            double meanAnnualNightSum = 0.0;
            const int historyCount = qMin(minHistory.size(), maxHistory.size());
            for (int i = 0; i < historyCount; ++i) {
                const double historyMin = minHistory.at(i);
                const double historyMax = maxHistory.at(i);
                minAnnual = qMin(minAnnual, historyMin);
                maxAnnual = qMax(maxAnnual, historyMax);
                meanAnnualSum += 0.5 * (historyMin + historyMax);
                meanAnnualDaySum += historyMax;
                meanAnnualNightSum += historyMin;
            }
            if (historyCount == 0) {
                minAnnual = dailyMin;
                maxAnnual = dailyMax;
            }
            const double meanAnnual = (historyCount > 0) ? meanAnnualSum / historyCount : dailyMean;
            const double meanAnnualDay =
                (historyCount > 0) ? meanAnnualDaySum / historyCount : dailyMax;
            const double meanAnnualNight =
                (historyCount > 0) ? meanAnnualNightSum / historyCount : dailyMin;

            TemperatureSummaryPoint summaryPoint;
            summaryPoint.latitudeDegrees = latitude;
            summaryPoint.minimumKelvin = minAnnual;
            summaryPoint.maximumKelvin = maxAnnual;
            summaryPoint.meanAnnualKelvin = meanAnnual;
            summaryPoint.meanAnnualDayKelvin = meanAnnualDay;
            summaryPoint.meanAnnualNightKelvin = meanAnnualNight;
            summaryPoint.minimumCelsius = minAnnual - kKelvinOffset;
            summaryPoint.maximumCelsius = maxAnnual - kKelvinOffset;
            summaryPoint.meanAnnualCelsius = meanAnnual - kKelvinOffset;
            summaryPoint.meanAnnualDayCelsius = meanAnnualDay - kKelvinOffset;
            summaryPoint.meanAnnualNightCelsius = meanAnnualNight - kKelvinOffset;
            summaryPoints.push_back(summaryPoint);
        }

        temperaturePlot_->setTemperatureSeries(dailyPoints,
                                               summaryPoints,
                                               QString(),
                                               rotationMode,
                                               hasAtmosphere);

        surfaceTemperatureAggregation_.dailyMinimums.fill(0.0);
        surfaceTemperatureAggregation_.dailyMaximums.fill(0.0);
    }

    void updateSurfacePointStatusDialog() {
        if (!surfacePointStatusDialog_) {
            return;
        }
        const SurfacePoint *point = surfaceGrid_.pointAt(selectedSurfacePointIndex_);
        if (!point) {
            surfacePointStatusDialog_->clearPoint();
            return;
        }
        const double tileAreaKm2 = surfaceGrid_.tileAreaKm2(selectedSurfacePointIndex_);
        const double tileEdgeLengthKm = surfaceGrid_.tileEdgeLengthKm(selectedSurfacePointIndex_);
        surfacePointStatusDialog_->setPoint(*point, tileAreaKm2, tileEdgeLengthKm);
        const auto &atmosphereGrid = surfaceGrid_.atmosphericGrid();
        const AtmosphericColumn *column = nullptr;
        if (selectedSurfacePointIndex_ >= 0 &&
            selectedSurfacePointIndex_ < atmosphereGrid.columnCount()) {
            column = &atmosphereGrid.columns().at(selectedSurfacePointIndex_);
        }
        surfacePointStatusDialog_->setAtmosphereProfile(column);
    }

    void setInputValue(QLineEdit *input, double value) {
        input->setText(QString::number(value));
    }

    void setPlanetPresets(const QVector<PlanetPreset> &planets,
                          const QString &selectedPlanetName = QString()) {
        const QSignalBlocker blocker(planetComboBox_);
        planetComboBox_->clear();
        presetPlanetNames_.clear();
        for (const auto &planet : planets) {
            presetPlanetNames_.insert(planet.name);
            addPlanetItem(planet, false);
        }
        if (selectedPlanetName.isEmpty()) {
            planetComboBox_->setCurrentIndex(-1);
        } else {
            const int selectedIndex = findPlanetIndexByName(selectedPlanetName);
            planetComboBox_->setCurrentIndex(selectedIndex);
        }
        updatePlanetSemiMajorAxisLabel();
        updatePlanetDayLengthLabel();
        updatePlanetYearLengthLabel();
        updatePlanetMassLabel();
        updatePlanetRadiusLabel();
        updatePlanetDerivedLabels();
        updateAtmospherePlanetParameters();
        updateAtmosphereComposition();
        updatePlanetOrbitLabels();
        syncMaterialWithPlanet();
        syncRotationModeWithPlanet();
        syncHeightSeedWithPlanet();
        syncHeightmapScaleWithPlanet();
        syncFlatHeightWithPlanet();
        syncAtmosphereDisableWithPlanet();
        syncCloudAlbedoWithPlanet();
        syncDiurnalCoolingBiasWithPlanet();
        syncManualGreenhouseOnTopWithPlanet();
        syncAdvancedRadiationWithPlanet();
        applyStellarPresetForCurrentPlanet();
        updatePlanetActions();
        refreshStarSystemView();
    }

    void clearPlanetPresets() {
        const QSignalBlocker blocker(planetComboBox_);
        planetComboBox_->clear();
        presetPlanetNames_.clear();
        planetSemiMajorAxisLabel_->setText(QStringLiteral("—"));
        planetDayLengthLabel_->setText(QStringLiteral("—"));
        planetYearLengthLabel_->setText(QStringLiteral("—"));
        planetMassLabel_->setText(QStringLiteral("—"));
        planetRadiusLabel_->setText(QStringLiteral("—"));
        planetSurfaceGravityLabel_->setText(QStringLiteral("—"));
        planetSurfaceAreaLabel_->setText(QStringLiteral("—"));
        planetEccentricityLabel_->setText(QStringLiteral("—"));
        planetObliquityLabel_->setText(QStringLiteral("—"));
        planetPerihelionArgumentLabel_->setText(QStringLiteral("—"));
        updateAtmospherePlanetParameters();
        updateAtmosphereComposition();
        {
            const QSignalBlocker rotationBlocker(rotationModeComboBox_);
            rotationModeComboBox_->setCurrentIndex(-1);
        }
        updateRotationModeIllustration();
        if (heightSeedSpinBox_) {
            const QSignalBlocker seedBlocker(heightSeedSpinBox_);
            heightSeedSpinBox_->setValue(0);
        }
        if (heightmapScaleSpinBox_) {
            const QSignalBlocker scaleBlocker(heightmapScaleSpinBox_);
            heightmapScaleSpinBox_->setValue(0.0);
            heightmapScaleSpinBox_->setEnabled(false);
        }
        if (heightmapImportButton_) {
            heightmapImportButton_->setEnabled(false);
        }
        if (flatHeightButton_) {
            const QSignalBlocker flatBlocker(flatHeightButton_);
            flatHeightButton_->setChecked(false);
            flatHeightButton_->setEnabled(false);
            updateFlatHeightButtonText(false);
        }
        if (disableAtmosphereButton_) {
            const QSignalBlocker atmosphereBlocker(disableAtmosphereButton_);
            disableAtmosphereButton_->setChecked(false);
            disableAtmosphereButton_->setEnabled(false);
            updateAtmosphereButtonText(false);
        }
        if (cloudAlbedoSpinBox_) {
            const QSignalBlocker cloudBlocker(cloudAlbedoSpinBox_);
            cloudAlbedoSpinBox_->setValue(0.0);
        }
        if (diurnalCoolingBiasSpinBox_) {
            const QSignalBlocker diurnalBlocker(diurnalCoolingBiasSpinBox_);
            diurnalCoolingBiasSpinBox_->setValue(0.0);
        }
        if (manualGreenhouseOnTopCheckBox_) {
            const QSignalBlocker greenhouseBlocker(manualGreenhouseOnTopCheckBox_);
            manualGreenhouseOnTopCheckBox_->setChecked(false);
        }
        if (advancedRadiationCheckBox_) {
            const QSignalBlocker advancedBlocker(advancedRadiationCheckBox_);
            advancedRadiationCheckBox_->setChecked(false);
        }
        updatePlanetActions();
        updateTemperaturePlot();
        refreshStarSystemView();
    }

    void refreshStarSystemView() {
        updateStarSystemOrbits();
    }

    void updateStarSystemSelection() {
        if (!starSystemTopView_) {
            return;
        }
        starSystemTopView_->setSelectedIndex(planetComboBox_ ? planetComboBox_->currentIndex() : -1);
    }

    double resolveStarSystemElapsedDays() {
        // Если есть прогресс поверхностной симуляции, синхронизируем орбитальный виджет
        // с этой шкалой времени, чтобы положение планет отражало текущую «симуляционную» дату.
        if (surfaceSimRunning_ || surfaceTime_.orbitTotalHourIndex() > 0) {
            return surfaceTime_.orbitElapsedDays();
        }
        return starSystemElapsedDays_;
    }

    void updateStarSystemOrbits() {
        if (!starSystemTopView_ || !planetComboBox_) {
            return;
        }

        StellarParameters primary{};
        bool hasPrimary = readStellarParametersForRender(radiusInput_, temperatureInput_, primary);
        if (!hasPrimary && planetComboBox_->currentIndex() >= 0) {
            // Используем параметры пресета, если поля ввода звезды пустые или некорректны.
            bool radiusOk = false;
            bool temperatureOk = false;
            const double radius =
                planetComboBox_->currentData(kRolePrimaryStarRadius).toDouble(&radiusOk);
            const double temperature =
                planetComboBox_->currentData(kRolePrimaryStarTemperature).toDouble(&temperatureOk);
            if (radiusOk && temperatureOk && radius > 0.0 && temperature > 0.0) {
                primary.radiusInSolarRadii = radius;
                primary.temperatureKelvin = temperature;
                hasPrimary = true;
            }
        }
        if (hasPrimary) {
            starSystemTopView_->setStarParameters(primary.temperatureKelvin,
                                                  primary.radiusInSolarRadii);
        }
        StellarParameters secondary{};
        bool hasSecondary = false;
        if (secondStarCheckBox_ && secondStarCheckBox_->isChecked()) {
            hasSecondary = readStellarParametersForRender(secondaryRadiusInput_,
                                                          secondaryTemperatureInput_,
                                                          secondary);
        }
        if (!hasSecondary && planetComboBox_->currentIndex() >= 0) {
            bool radiusOk = false;
            bool temperatureOk = false;
            const double radius =
                planetComboBox_->currentData(kRoleSecondaryStarRadius).toDouble(&radiusOk);
            const double temperature =
                planetComboBox_->currentData(kRoleSecondaryStarTemperature).toDouble(&temperatureOk);
            if (radiusOk && temperatureOk && radius > 0.0 && temperature > 0.0) {
                secondary.radiusInSolarRadii = radius;
                secondary.temperatureKelvin = temperature;
                hasSecondary = true;
            }
        }
        if (hasSecondary) {
            starSystemTopView_->setSecondaryStarParameters(secondary.temperatureKelvin,
                                                           secondary.radiusInSolarRadii);
        } else {
            starSystemTopView_->setSecondaryStarParameters(0.0, 0.0);
        }

        const double elapsedDays = resolveStarSystemElapsedDays();
        if (surfaceSimRunning_ || surfaceTime_.orbitTotalHourIndex() > 0) {
            // Держим локальный таймер синхронизированным с симуляцией,
            // чтобы после остановки времени орбиты не перескакивали назад.
            starSystemElapsedDays_ = elapsedDays;
        }

        QVector<StarSystemTopView::PlanetOrbit> planets;
        planets.reserve(planetComboBox_->count());
        for (int i = 0; i < planetComboBox_->count(); ++i) {
            StarSystemTopView::PlanetOrbit orbit;
            orbit.name = planetComboBox_->itemData(i, kRolePlanetName).toString();
            orbit.semiMajorAxisAu = planetComboBox_->itemData(i, kRoleSemiMajorAxis).toDouble();
            orbit.eccentricity = planetComboBox_->itemData(i, kRoleEccentricity).toDouble();
            orbit.perihelionArgumentDegrees =
                planetComboBox_->itemData(i, kRolePerihelionArgument).toDouble();
            const double yearLengthDays =
                planetComboBox_->itemData(i, kRoleYearLengthDays).toDouble();
            const OrbitSegment segment = resolvePlanetOrbitSegment(elapsedDays,
                                                                  yearLengthDays,
                                                                  orbit.semiMajorAxisAu,
                                                                  orbit.eccentricity);
            orbit.trueAnomalyRadians = segment.trueAnomalyRadians;
            planets.push_back(orbit);
        }

        starSystemTopView_->setPlanets(planets);
        StarSystemTopView::BinaryOrbitParameters binaryParameters;
        const bool hasBinaryOrbit =
            planetComboBox_->currentData(kRoleHasBinaryOrbit).toBool();
        if (hasBinaryOrbit) {
            binaryParameters.semiMajorAxisAu =
                planetComboBox_->currentData(kRoleBinarySemiMajorAxis).toDouble();
            binaryParameters.periodDays =
                planetComboBox_->currentData(kRoleBinaryPeriodDays).toDouble();
            binaryParameters.eccentricity =
                planetComboBox_->currentData(kRoleBinaryEccentricity).toDouble();
            binaryParameters.inclinationDegrees =
                planetComboBox_->currentData(kRoleBinaryInclinationDegrees).toDouble();
            binaryParameters.argumentPericenterDegreesA =
                planetComboBox_->currentData(kRoleBinaryArgumentPericenterA).toDouble();
            binaryParameters.argumentPericenterDegreesB =
                planetComboBox_->currentData(kRoleBinaryArgumentPericenterB).toDouble();
        }
        starSystemTopView_->setBinarySystemParameters(binaryParameters);
        starSystemTopView_->setBinarySystemElapsedDays(elapsedDays);
        updateStarSystemSelection();
        updateSurfaceOrbitLighting();
    }

    QString formatPlanetName(const PlanetPreset &planet) const {
        return QStringLiteral("%1 (%2 а.е.)")
            .arg(planet.name, formatSemiMajorAxis(planet.semiMajorAxis));
    }

    QString formatSemiMajorAxis(double value) const {
        return QLocale().toString(value, 'f', 2);
    }

    QString formatDistance(double value) const {
        return QLocale().toString(value, 'f', 3);
    }

    QString formatDayLength(double value) const {
        return QLocale().toString(value, 'f', 2);
    }

    QString formatMass(double value) const {
        return QLocale().toString(value, 'f', 3);
    }

    QString formatRadius(double value) const {
        return QLocale().toString(value, 'f', 1);
    }

    QString formatSurfaceGravity(double value) const {
        return QLocale().toString(value, 'f', 2);
    }

    QString formatSurfaceArea(double value) const {
        return QLocale().toString(value, 'f', 0);
    }

    QString formatEccentricity(double value) const {
        return QLocale().toString(value, 'f', 3);
    }

    QString formatAngle(double value) const {
        return QLocale().toString(value, 'f', 1);
    }

    void updatePlanetSemiMajorAxisLabel() {
        const QVariant value = planetComboBox_->currentData(kRoleSemiMajorAxis);
        if (!value.isValid()) {
            planetSemiMajorAxisLabel_->setText(QStringLiteral("—"));
            return;
        }
        planetSemiMajorAxisLabel_->setText(formatSemiMajorAxis(value.toDouble()));
    }

    void updatePlanetDayLengthLabel() {
        const QVariant value = planetComboBox_->currentData(kRoleDayLength);
        if (!value.isValid()) {
            planetDayLengthLabel_->setText(QStringLiteral("—"));
            return;
        }
        planetDayLengthLabel_->setText(formatDayLength(value.toDouble()));
    }

    void updatePlanetYearLengthLabel() {
        const QVariant value = planetComboBox_->currentData(kRoleYearLengthDays);
        if (!value.isValid()) {
            planetYearLengthLabel_->setText(QStringLiteral("—"));
            return;
        }
        planetYearLengthLabel_->setText(formatDayLength(value.toDouble()));
    }

    void updatePlanetMassLabel() {
        const QVariant value = planetComboBox_->currentData(kRoleMassEarths);
        if (!value.isValid()) {
            planetMassLabel_->setText(QStringLiteral("—"));
            return;
        }
        planetMassLabel_->setText(formatMass(value.toDouble()));
    }

    void updatePlanetRadiusLabel() {
        const QVariant value = planetComboBox_->currentData(kRoleRadiusKm);
        if (!value.isValid()) {
            planetRadiusLabel_->setText(QStringLiteral("—"));
            return;
        }
        planetRadiusLabel_->setText(formatRadius(value.toDouble()));
    }

    void updatePlanetDerivedLabels() {
        updatePlanetSurfaceGravityLabel();
        updatePlanetSurfaceAreaLabel();
    }

    void updatePlanetSurfaceGravityLabel() {
        const QVariant massValue = planetComboBox_->currentData(kRoleMassEarths);
        const QVariant radiusValue = planetComboBox_->currentData(kRoleRadiusKm);
        if (!massValue.isValid() || !radiusValue.isValid()) {
            planetSurfaceGravityLabel_->setText(QStringLiteral("—"));
            return;
        }
        const double massEarths = massValue.toDouble();
        const double radiusKm = radiusValue.toDouble();
        if (massEarths <= 0.0 || radiusKm <= 0.0) {
            planetSurfaceGravityLabel_->setText(QStringLiteral("—"));
            return;
        }
        // Относительное g в земных единицах: g/g⊕ = (M/M⊕) / (R/R⊕)².
        const double radiusEarths = radiusKm / kEarthRadiusKm;
        const double surfaceGravityEarths = massEarths / (radiusEarths * radiusEarths);
        planetSurfaceGravityLabel_->setText(formatSurfaceGravity(surfaceGravityEarths));
    }

    void updatePlanetSurfaceAreaLabel() {
        const QVariant radiusValue = planetComboBox_->currentData(kRoleRadiusKm);
        if (!radiusValue.isValid()) {
            planetSurfaceAreaLabel_->setText(QStringLiteral("—"));
            return;
        }
        const double radiusKm = radiusValue.toDouble();
        if (radiusKm <= 0.0) {
            planetSurfaceAreaLabel_->setText(QStringLiteral("—"));
            return;
        }
        // Площадь поверхности сферы: S = 4πR².
        const double surfaceArea = 4.0 * M_PI * radiusKm * radiusKm;
        planetSurfaceAreaLabel_->setText(formatSurfaceArea(surfaceArea));
    }

    void updateAtmospherePlanetParameters() {
        if (!atmosphereWidget_) {
            return;
        }
        const QVariant massValue = planetComboBox_->currentData(kRoleMassEarths);
        const QVariant radiusValue = planetComboBox_->currentData(kRoleRadiusKm);
        if (!massValue.isValid() || !radiusValue.isValid()) {
            atmosphereWidget_->clearPlanetParameters();
            return;
        }
        atmosphereWidget_->setPlanetParameters(massValue.toDouble(), radiusValue.toDouble());
    }

    void updateAtmosphereComposition() {
        if (!atmosphereWidget_) {
            return;
        }
        const int index = planetComboBox_ ? planetComboBox_->currentIndex() : -1;
        const bool atmosphereDisabled = isAtmosphereDisabledForIndex(index);
        const QVariant compositionValue = planetComboBox_->currentData(kRoleAtmosphere);
        const QVariant bottomLayerValue =
            planetComboBox_->currentData(kRoleAtmosphereBottomLayerThickness);
        const QVariant minTopPressureValue =
            planetComboBox_->currentData(kRoleMinTopPressureAtm);
        const QVariant mixingValue =
            planetComboBox_->currentData(kRoleVerticalWindMixingCoefficient);
        const double bottomLayerThickness = bottomLayerValue.isValid()
            ? bottomLayerValue.toDouble()
            : kDefaultAtmosphereBottomLayerThicknessMeters;
        const double minTopPressureAtm = minTopPressureValue.isValid()
            ? minTopPressureValue.toDouble()
            : kDefaultMinTopPressureAtm;
        const double verticalWindMixingCoefficient = mixingValue.isValid()
            ? mixingValue.toDouble()
            : kDefaultVerticalWindMixingCoefficientKz;
        if (planetComboBox_ && index >= 0 && !bottomLayerValue.isValid()) {
            planetComboBox_->setItemData(index,
                                         bottomLayerThickness,
                                         kRoleAtmosphereBottomLayerThickness);
        }
        if (planetComboBox_ && index >= 0 && !minTopPressureValue.isValid()) {
            planetComboBox_->setItemData(index,
                                         minTopPressureAtm,
                                         kRoleMinTopPressureAtm);
        }
        if (planetComboBox_ && index >= 0 && !mixingValue.isValid()) {
            planetComboBox_->setItemData(index,
                                         verticalWindMixingCoefficient,
                                         kRoleVerticalWindMixingCoefficient);
        }
        if (!compositionValue.isValid()) {
            atmosphereWidget_->setComposition(AtmosphereComposition{});
            atmosphereWidget_->setMinBottomLayerThicknessMeters(bottomLayerThickness);
            atmosphereWidget_->setMinTopPressureAtm(minTopPressureAtm);
            atmosphereWidget_->setVerticalWindMixingCoefficient(verticalWindMixingCoefficient);
            atmosphereWidget_->setEnabled(!atmosphereDisabled);
            return;
        }
        atmosphereWidget_->setComposition(compositionValue.value<AtmosphereComposition>());
        atmosphereWidget_->setMinBottomLayerThicknessMeters(bottomLayerThickness);
        atmosphereWidget_->setMinTopPressureAtm(minTopPressureAtm);
        atmosphereWidget_->setVerticalWindMixingCoefficient(verticalWindMixingCoefficient);
        atmosphereWidget_->setEnabled(!atmosphereDisabled);
    }

    bool isAtmosphereDisabledForIndex(int index) const {
        if (!planetComboBox_ || index < 0) {
            return false;
        }
        return planetComboBox_->itemData(index, kRoleAtmosphereDisabled).toBool();
    }

    AtmosphereComposition currentAtmosphereForCalculations() const {
        AtmosphereComposition atmosphere;
        if (!planetComboBox_) {
            return atmosphere;
        }
        const int index = planetComboBox_->currentIndex();
        if (isAtmosphereDisabledForIndex(index)) {
            // Атмосфера может быть временно отключена кнопкой, не теряя заданный состав.
            return atmosphere;
        }
        const QVariant atmosphereValue = planetComboBox_->currentData(kRoleAtmosphere);
        if (atmosphereValue.isValid()) {
            atmosphere = atmosphereValue.value<AtmosphereComposition>();
        }
        return atmosphere;
    }

    double currentAtmosphereBottomLayerThicknessMeters() const {
        if (!planetComboBox_) {
            return kDefaultAtmosphereBottomLayerThicknessMeters;
        }
        const QVariant thicknessValue =
            planetComboBox_->currentData(kRoleAtmosphereBottomLayerThickness);
        if (!thicknessValue.isValid()) {
            return kDefaultAtmosphereBottomLayerThicknessMeters;
        }
        return qMax(0.0, thicknessValue.toDouble());
    }

    double currentVerticalWindMixingCoefficientKz() const {
        if (!planetComboBox_) {
            return kDefaultVerticalWindMixingCoefficientKz;
        }
        const QVariant mixingValue =
            planetComboBox_->currentData(kRoleVerticalWindMixingCoefficient);
        if (!mixingValue.isValid()) {
            return kDefaultVerticalWindMixingCoefficientKz;
        }
        return qMax(0.0, mixingValue.toDouble());
    }

    double currentMinDenseAtmosphereTemperatureKelvin() const {
        if (!planetComboBox_) {
            return 0.0;
        }
        const QVariant temperatureValue =
            planetComboBox_->currentData(kRoleMinDenseAtmosphereTemperature);
        if (!temperatureValue.isValid()) {
            return 0.0;
        }
        return qMax(0.0, temperatureValue.toDouble());
    }

    double currentMinTopPressureAtm() const {
        if (!planetComboBox_) {
            return kDefaultMinTopPressureAtm;
        }
        const QVariant pressureValue =
            planetComboBox_->currentData(kRoleMinTopPressureAtm);
        if (!pressureValue.isValid()) {
            return kDefaultMinTopPressureAtm;
        }
        return qMax(0.0, pressureValue.toDouble());
    }

    int latitudeStepDegrees() const {
        return 1;
    }

    int latitudePoints() const {
        return 180 / latitudeStepDegrees() + 1;
    }

    void updatePlanetOrbitLabels() {
        const QVariant eccentricity = planetComboBox_->currentData(kRoleEccentricity);
        const QVariant obliquity = planetComboBox_->currentData(kRoleObliquity);
        const QVariant perihelionArgument = planetComboBox_->currentData(kRolePerihelionArgument);
        if (!eccentricity.isValid() || !obliquity.isValid() || !perihelionArgument.isValid()) {
            planetEccentricityLabel_->setText(QStringLiteral("—"));
            planetObliquityLabel_->setText(QStringLiteral("—"));
            planetPerihelionArgumentLabel_->setText(QStringLiteral("—"));
            if (surfaceGlobeWidget_) {
                surfaceGlobeWidget_->setAxisTiltDegrees(0.0);
            }
            return;
        }
        planetEccentricityLabel_->setText(formatEccentricity(eccentricity.toDouble()));
        planetObliquityLabel_->setText(formatAngle(obliquity.toDouble()));
        planetPerihelionArgumentLabel_->setText(formatAngle(perihelionArgument.toDouble()));
        if (surfaceGlobeWidget_) {
            surfaceGlobeWidget_->setAxisTiltDegrees(obliquity.toDouble());
        }
    }

    bool hasPrimaryInputs() const {
        return !radiusInput_->text().trimmed().isEmpty() &&
               !temperatureInput_->text().trimmed().isEmpty();
    }

    bool hasSecondaryInputs() const {
        return !secondaryRadiusInput_->text().trimmed().isEmpty() &&
               !secondaryTemperatureInput_->text().trimmed().isEmpty();
    }

    void applyStellarPresetForCurrentPlanet() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || isCustomPlanetIndex(index)) {
            return;
        }

        const QVariant primaryRadiusValue =
            planetComboBox_->itemData(index, kRolePrimaryStarRadius);
        const QVariant primaryTemperatureValue =
            planetComboBox_->itemData(index, kRolePrimaryStarTemperature);
        if (!primaryRadiusValue.isValid() || !primaryTemperatureValue.isValid()) {
            return;
        }

        setInputValue(radiusInput_, primaryRadiusValue.toDouble());
        setInputValue(temperatureInput_, primaryTemperatureValue.toDouble());

        const bool hasSecondary =
            planetComboBox_->itemData(index, kRoleHasSecondaryStar).toBool();
        if (!hasSecondary) {
            secondStarCheckBox_->setChecked(false);
            secondaryRadiusInput_->clear();
            secondaryTemperatureInput_->clear();
            return;
        }

        const QVariant secondaryRadiusValue =
            planetComboBox_->itemData(index, kRoleSecondaryStarRadius);
        const QVariant secondaryTemperatureValue =
            planetComboBox_->itemData(index, kRoleSecondaryStarTemperature);
        secondStarCheckBox_->setChecked(true);
        if (secondaryRadiusValue.isValid()) {
            setInputValue(secondaryRadiusInput_, secondaryRadiusValue.toDouble());
        } else {
            secondaryRadiusInput_->clear();
        }
        if (secondaryTemperatureValue.isValid()) {
            setInputValue(secondaryTemperatureInput_, secondaryTemperatureValue.toDouble());
        } else {
            secondaryTemperatureInput_->clear();
        }
    }

    void addPlanetItem(const PlanetPreset &planet, bool isCustom) {
        planetComboBox_->addItem(formatPlanetName(planet), planet.semiMajorAxis);
        const int index = planetComboBox_->count() - 1;
        planetComboBox_->setItemData(index, planet.semiMajorAxis, kRoleSemiMajorAxis);
        planetComboBox_->setItemData(index, planet.dayLengthDays, kRoleDayLength);
        planetComboBox_->setItemData(index, planet.yearLengthDays, kRoleYearLengthDays);
        planetComboBox_->setItemData(index, planet.eccentricity, kRoleEccentricity);
        planetComboBox_->setItemData(index, planet.obliquityDegrees, kRoleObliquity);
        planetComboBox_->setItemData(index, planet.perihelionArgumentDegrees, kRolePerihelionArgument);
        planetComboBox_->setItemData(index, planet.massEarths, kRoleMassEarths);
        planetComboBox_->setItemData(index, planet.radiusKm, kRoleRadiusKm);
        const RotationMode rotationMode =
            planet.tidallyLocked ? RotationMode::TidalLocked : RotationMode::Normal;
        planetComboBox_->setItemData(index, static_cast<int>(rotationMode), kRoleRotationMode);
        planetComboBox_->setItemData(index, planet.rotationModeOverride, kRoleRotationModeOverride);
        planetComboBox_->setItemData(index, isCustom, kRoleIsCustom);
        planetComboBox_->setItemData(index, planet.name, kRolePlanetName);
        planetComboBox_->setItemData(index, planet.spinOrbitP, kRoleSpinOrbitP);
        planetComboBox_->setItemData(index, planet.spinOrbitQ, kRoleSpinOrbitQ);
        planetComboBox_->setItemData(index, planet.surfaceMaterialId, kRoleMaterialId);
        planetComboBox_->setItemData(index, QVariant::fromValue(planet.atmosphere), kRoleAtmosphere);
        planetComboBox_->setItemData(index, planet.greenhouseOpacity, kRoleGreenhouseOpacity);
        planetComboBox_->setItemData(index,
                                     planet.manualGreenhouseOnTopOfAtmosphere,
                                     kRoleManualGreenhouseOnTopOfAtmosphere);
        planetComboBox_->setItemData(index,
                                     static_cast<int>(RadiationModelType::Layered),
                                     kRoleAdvancedRadiationModel);
        planetComboBox_->setItemData(index, planet.cloudAlbedo, kRoleCloudAlbedo);
        planetComboBox_->setItemData(index, planet.disableWaterAndClouds, kRoleDisableWaterAndClouds);
        planetComboBox_->setItemData(index, 0.0, kRoleDiurnalCoolingBiasK);
        planetComboBox_->setItemData(index, planet.geothermalFluxWPerM2, kRoleGeothermalFlux);
        planetComboBox_->setItemData(index, static_cast<int>(planet.heightSourceType),
                                     kRoleHeightSourceType);
        planetComboBox_->setItemData(index, planet.primaryStar.radiusInSolarRadii,
                                     kRolePrimaryStarRadius);
        planetComboBox_->setItemData(index, planet.primaryStar.temperatureKelvin,
                                     kRolePrimaryStarTemperature);
        planetComboBox_->setItemData(index, static_cast<bool>(planet.secondaryStar),
                                     kRoleHasSecondaryStar);
        if (planet.secondaryStar) {
            planetComboBox_->setItemData(index,
                                         planet.secondaryStar->radiusInSolarRadii,
                                         kRoleSecondaryStarRadius);
            planetComboBox_->setItemData(index,
                                         planet.secondaryStar->temperatureKelvin,
                                         kRoleSecondaryStarTemperature);
        } else {
            planetComboBox_->setItemData(index, QVariant(), kRoleSecondaryStarRadius);
            planetComboBox_->setItemData(index, QVariant(), kRoleSecondaryStarTemperature);
        }
        planetComboBox_->setItemData(index, planet.heightmapPath, kRoleHeightmapPath);
        planetComboBox_->setItemData(index, planet.heightmapScaleKm, kRoleHeightmapScaleKm);
        planetComboBox_->setItemData(index, planet.heightSeed, kRoleHeightSeed);
        planetComboBox_->setItemData(index, planet.useContinentsHeight, kRoleUseContinentsHeight);
        planetComboBox_->setItemData(index, planet.hasSeaLevel, kRoleHasSeaLevel);
        planetComboBox_->setItemData(index, false, kRoleFlatHeight);
        planetComboBox_->setItemData(index, false, kRoleAtmosphereDisabled);
        planetComboBox_->setItemData(index, planet.starPresetId, kRoleStarPresetId);
        planetComboBox_->setItemData(index, planet.hasBinaryOrbit, kRoleHasBinaryOrbit);
        planetComboBox_->setItemData(index, planet.binarySemiMajorAxisAu, kRoleBinarySemiMajorAxis);
        planetComboBox_->setItemData(index, planet.binaryPeriodDays, kRoleBinaryPeriodDays);
        planetComboBox_->setItemData(index, planet.binaryEccentricity, kRoleBinaryEccentricity);
        planetComboBox_->setItemData(index,
                                     planet.binaryInclinationDegrees,
                                     kRoleBinaryInclinationDegrees);
        planetComboBox_->setItemData(index,
                                     planet.binaryArgumentPericenterDegreesA,
                                     kRoleBinaryArgumentPericenterA);
        planetComboBox_->setItemData(index,
                                     planet.binaryArgumentPericenterDegreesB,
                                     kRoleBinaryArgumentPericenterB);
        planetComboBox_->setItemData(index,
                                     planet.minBottomLayerThicknessMeters,
                                     kRoleAtmosphereBottomLayerThickness);
        const QVariant minTopPressureValue =
            (planet.minTopPressureAtm > 0.0)
                ? QVariant(planet.minTopPressureAtm)
                : QVariant();
        planetComboBox_->setItemData(index,
                                     minTopPressureValue,
                                     kRoleMinTopPressureAtm);
        planetComboBox_->setItemData(index,
                                     planet.minDenseAtmosphereTemperatureK,
                                     kRoleMinDenseAtmosphereTemperature);
        planetComboBox_->setItemData(index,
                                     planet.verticalWindMixingCoefficient,
                                     kRoleVerticalWindMixingCoefficient);
    }

    bool isCustomPlanetIndex(int index) const {
        return planetComboBox_->itemData(index, kRoleIsCustom).toBool();
    }

    int findPlanetIndexByName(const QString &name) const {
        for (int i = 0; i < planetComboBox_->count(); ++i) {
            if (planetComboBox_->itemData(i, kRolePlanetName).toString() == name) {
                return i;
            }
        }
        return -1;
    }

    void updatePlanetActions() {
        const int index = planetComboBox_->currentIndex();
        deletePlanetButton_->setVisible(index >= 0 && isCustomPlanetIndex(index));
    }

    void onAddPlanetRequested() {
        QDialog dialog(this);
        dialog.setWindowTitle(QStringLiteral("Добавить планету"));

        auto *nameInput = new QLineEdit(&dialog);
        nameInput->setPlaceholderText(QStringLiteral("Название"));

        auto *axisInput = new QLineEdit(&dialog);
        axisInput->setPlaceholderText(QStringLiteral("Например, 1.0"));
        auto *validator = new QDoubleValidator(0.0, std::numeric_limits<double>::max(), 10, &dialog);
        validator->setNotation(QDoubleValidator::StandardNotation);
        validator->setLocale(QLocale::C);
        axisInput->setValidator(validator);

        auto *dayLengthInput = new QLineEdit(&dialog);
        dayLengthInput->setPlaceholderText(QStringLiteral("Например, 1.0"));
        dayLengthInput->setValidator(validator);

        auto *yearLengthInput = new QLineEdit(&dialog);
        yearLengthInput->setPlaceholderText(QStringLiteral("Например, 365.25"));
        yearLengthInput->setValidator(validator);

        auto *massInput = new QLineEdit(&dialog);
        massInput->setPlaceholderText(QStringLiteral("Например, 1.0"));
        massInput->setValidator(validator);

        auto *radiusInput = new QLineEdit(&dialog);
        radiusInput->setPlaceholderText(QStringLiteral("Например, 6371"));
        radiusInput->setValidator(validator);

        auto *eccentricityInput = new QLineEdit(&dialog);
        eccentricityInput->setPlaceholderText(QStringLiteral("Например, 0.0167"));
        auto *eccentricityValidator = new QDoubleValidator(0.0, 0.999, 6, &dialog);
        eccentricityValidator->setNotation(QDoubleValidator::StandardNotation);
        eccentricityValidator->setLocale(QLocale::C);
        eccentricityInput->setValidator(eccentricityValidator);

        auto *obliquityInput = new QLineEdit(&dialog);
        obliquityInput->setPlaceholderText(QStringLiteral("Например, 23.44"));
        auto *obliquityValidator = new QDoubleValidator(0.0, 180.0, 4, &dialog);
        obliquityValidator->setNotation(QDoubleValidator::StandardNotation);
        obliquityValidator->setLocale(QLocale::C);
        obliquityInput->setValidator(obliquityValidator);

        auto *perihelionArgumentInput = new QLineEdit(&dialog);
        perihelionArgumentInput->setPlaceholderText(QStringLiteral("Например, 102.94"));
        auto *perihelionValidator = new QDoubleValidator(0.0, 360.0, 4, &dialog);
        perihelionValidator->setNotation(QDoubleValidator::StandardNotation);
        perihelionValidator->setLocale(QLocale::C);
        perihelionArgumentInput->setValidator(perihelionValidator);

        auto *greenhouseOpacityInput = new QLineEdit(&dialog);
        greenhouseOpacityInput->setPlaceholderText(QStringLiteral("Например, 0.0"));
        auto *greenhouseValidator = new QDoubleValidator(0.0, 0.999, 4, &dialog);
        greenhouseValidator->setNotation(QDoubleValidator::StandardNotation);
        greenhouseValidator->setLocale(QLocale::C);
        greenhouseOpacityInput->setValidator(greenhouseValidator);

        auto *cloudAlbedoInput = new QDoubleSpinBox(&dialog);
        cloudAlbedoInput->setRange(0.0, 1.0);
        cloudAlbedoInput->setDecimals(2);
        cloudAlbedoInput->setSingleStep(0.05);

        auto *diurnalCoolingBiasInput = new QDoubleSpinBox(&dialog);
        diurnalCoolingBiasInput->setRange(0.0, 200.0);
        diurnalCoolingBiasInput->setDecimals(1);
        diurnalCoolingBiasInput->setSingleStep(1.0);
        diurnalCoolingBiasInput->setSuffix(QStringLiteral(" K"));
        diurnalCoolingBiasInput->setToolTip(
            QStringLiteral("Грубая поправка на суточное охлаждение поверхности."));

        auto *manualGreenhouseOnTopInput = new QCheckBox(
            QStringLiteral("Добавлять поверх атмосферной модели"), &dialog);
        manualGreenhouseOnTopInput->setToolTip(
            QStringLiteral("Использовать значение парниковой непрозрачности как дополнительный\n"
                           "фактор поверх атмосферной модели (например, для проверки гипотез)."));

        auto *heightSeedInput = new QSpinBox(&dialog);
        heightSeedInput->setRange(0, std::numeric_limits<int>::max());
        heightSeedInput->setValue(0);

        auto *formLayout = new QFormLayout();
        formLayout->addRow(QStringLiteral("Имя:"), nameInput);
        formLayout->addRow(QStringLiteral("Большая полуось (а.е.):"), axisInput);
        formLayout->addRow(QStringLiteral("Длина суток (земн. дни):"), dayLengthInput);
        formLayout->addRow(QStringLiteral("Длина года (земн. дни):"), yearLengthInput);
        formLayout->addRow(QStringLiteral("Масса (в массах Земли):"), massInput);
        formLayout->addRow(QStringLiteral("Радиус (км):"), radiusInput);
        formLayout->addRow(QStringLiteral("Эксцентриситет:"), eccentricityInput);
        formLayout->addRow(QStringLiteral("Наклон оси (°):"), obliquityInput);
        formLayout->addRow(QStringLiteral("Аргумент перицентра (°):"), perihelionArgumentInput);
        formLayout->addRow(QStringLiteral("Парниковая непрозрачность (0..1):"),
                           greenhouseOpacityInput);
        formLayout->addRow(QStringLiteral("Парниковая непрозрачность поверх атмосферы:"),
                           manualGreenhouseOnTopInput);
        formLayout->addRow(QStringLiteral("Альбедо облаков (0..1):"), cloudAlbedoInput);
        formLayout->addRow(QStringLiteral("Поправка суточного охлаждения (К):"),
                           diurnalCoolingBiasInput);
        formLayout->addRow(QStringLiteral("Семя рельефа:"), heightSeedInput);

        auto *materialInput = new QComboBox(&dialog);
        for (const auto &material : surfaceMaterials()) {
            materialInput->addItem(material.name, material.id);
        }
        formLayout->addRow(QStringLiteral("Материал поверхности:"), materialInput);

        auto *rotationModeInput = new QComboBox(&dialog);
        rotationModeInput->addItem(QStringLiteral("Обычное вращение (широта)"),
                                   static_cast<int>(RotationMode::Normal));
        rotationModeInput->addItem(QStringLiteral("Приливная синхронизация (угол от подсолнечной точки)"),
                                   static_cast<int>(RotationMode::TidalLocked));
        auto *spinOrbitPInput = new QSpinBox(&dialog);
        spinOrbitPInput->setRange(1, 20);
        spinOrbitPInput->setValue(1);
        auto *spinOrbitQInput = new QSpinBox(&dialog);
        spinOrbitQInput->setRange(1, 20);
        spinOrbitQInput->setValue(1);
        const QString spinOrbitTooltip =
            QStringLiteral("Резонанс p:q — отношение спиновой и орбитальной частот.\n"
                           "Безразмерные коэффициенты, 1:1 соответствует синхронному вращению.");
        spinOrbitPInput->setToolTip(spinOrbitTooltip);
        spinOrbitQInput->setToolTip(spinOrbitTooltip);
        auto *spinOrbitWidget = new QWidget(&dialog);
        auto *spinOrbitLayout = new QHBoxLayout(spinOrbitWidget);
        spinOrbitLayout->addWidget(spinOrbitPInput);
        spinOrbitLayout->addWidget(new QLabel(QStringLiteral(":"), spinOrbitWidget));
        spinOrbitLayout->addWidget(spinOrbitQInput);
        spinOrbitLayout->addStretch();
        spinOrbitLayout->setContentsMargins(0, 0, 0, 0);
        auto *atmosphereInput = new AtmosphereWidget(&dialog, true);
        formLayout->addRow(QStringLiteral("Режим вращения:"), rotationModeInput);
        formLayout->addRow(QStringLiteral("Резонанс вращения p:q:"), spinOrbitWidget);

        auto *formWidget = new QWidget(&dialog);
        formWidget->setLayout(formLayout);

        auto *contentLayout = new QHBoxLayout();
        contentLayout->addWidget(atmosphereInput, 1);
        contentLayout->addWidget(formWidget, 0);

        auto *dialogLayout = new QVBoxLayout(&dialog);
        dialogLayout->addLayout(contentLayout);

        const auto updateAtmosphereParameters = [massInput, radiusInput, atmosphereInput]() {
            bool massOk = false;
            bool radiusOk = false;
            const double massEarths = massInput->text().toDouble(&massOk);
            const double radiusKm = radiusInput->text().toDouble(&radiusOk);
            if (massOk && radiusOk && massEarths > 0.0 && radiusKm > 0.0) {
                atmosphereInput->setPlanetParameters(massEarths, radiusKm);
            } else {
                atmosphereInput->clearPlanetParameters();
            }
        };

        connect(massInput, &QLineEdit::textChanged, &dialog, updateAtmosphereParameters);
        connect(radiusInput, &QLineEdit::textChanged, &dialog, updateAtmosphereParameters);
        updateAtmosphereParameters();

        auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
        dialogLayout->addWidget(buttons);

        connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
        connect(buttons, &QDialogButtonBox::accepted, &dialog,
                [&dialog, nameInput, axisInput, dayLengthInput, yearLengthInput, massInput, radiusInput,
                 eccentricityInput, obliquityInput, perihelionArgumentInput,
                 greenhouseOpacityInput, manualGreenhouseOnTopInput, cloudAlbedoInput,
                 diurnalCoolingBiasInput,
                 heightSeedInput, materialInput,
                 rotationModeInput, spinOrbitPInput, spinOrbitQInput, atmosphereInput, this]() {
            const QString name = nameInput->text().trimmed();
            if (name.isEmpty()) {
                showInputError(QStringLiteral("Введите имя планеты."));
                return;
            }

            if (presetPlanetNames_.contains(name)) {
                showInputError(QStringLiteral("Нельзя добавлять планеты с именем из пресета."));
                return;
            }

            bool ok = false;
            const double axis = axisInput->text().toDouble(&ok);
            if (!ok || axis <= 0.0) {
                showInputError(QStringLiteral("Укажите положительную большую полуось орбиты планеты."));
                return;
            }

            const double dayLength = dayLengthInput->text().toDouble(&ok);
            if (!ok || dayLength <= 0.0) {
                showInputError(QStringLiteral("Укажите положительную длину суток планеты."));
                return;
            }

            const double yearLength = yearLengthInput->text().toDouble(&ok);
            if (!ok || yearLength <= 0.0) {
                showInputError(QStringLiteral("Укажите положительную длину года планеты."));
                return;
            }

            const double massEarths = massInput->text().toDouble(&ok);
            if (!ok || massEarths <= 0.0) {
                showInputError(QStringLiteral("Укажите положительную массу планеты в массах Земли."));
                return;
            }

            const double radiusKm = radiusInput->text().toDouble(&ok);
            if (!ok || radiusKm <= 0.0) {
                showInputError(QStringLiteral("Укажите положительный радиус планеты в километрах."));
                return;
            }

            const double eccentricity = eccentricityInput->text().toDouble(&ok);
            if (!ok || eccentricity < 0.0 || eccentricity >= 1.0) {
                showInputError(QStringLiteral("Укажите эксцентриситет от 0 до 1 (не включая 1)."));
                return;
            }

            const double obliquity = obliquityInput->text().toDouble(&ok);
            if (!ok || obliquity < 0.0 || obliquity > 180.0) {
                showInputError(QStringLiteral("Укажите наклон оси от 0 до 180 градусов."));
                return;
            }

            const double perihelionArgument = perihelionArgumentInput->text().toDouble(&ok);
            if (!ok || perihelionArgument < 0.0 || perihelionArgument >= 360.0) {
                showInputError(QStringLiteral("Укажите аргумент перицентра от 0 до 360 градусов."));
                return;
            }

            double greenhouseOpacity = 0.0;
            const QString greenhouseText = greenhouseOpacityInput->text().trimmed();
            if (!greenhouseText.isEmpty()) {
                greenhouseOpacity = greenhouseText.toDouble(&ok);
                if (!ok || greenhouseOpacity < 0.0 || greenhouseOpacity >= 1.0) {
                    showInputError(QStringLiteral("Укажите непрозрачность парникового слоя от 0 до 1."));
                    return;
                }
            }
            const bool manualGreenhouseOnTop = manualGreenhouseOnTopInput->isChecked();
            const double cloudAlbedo = cloudAlbedoInput->value();
            const double diurnalCoolingBiasK = diurnalCoolingBiasInput->value();
            const double geothermalFlux =
                subsurfaceGeothermalFluxSpinBox_ ? subsurfaceGeothermalFluxSpinBox_->value() : 0.0;

            const int existingIndex = findPlanetIndexByName(name);
            const QString materialId = materialInput->currentData().toString();
            const RotationMode rotationMode =
                static_cast<RotationMode>(rotationModeInput->currentData().toInt());
            const bool tidallyLocked = (rotationMode == RotationMode::TidalLocked);
            const int spinOrbitP = spinOrbitPInput->value();
            const int spinOrbitQ = spinOrbitQInput->value();
            const quint32 heightSeed = static_cast<quint32>(heightSeedInput->value());
            const AtmosphereComposition composition = atmosphereInput->composition(false);
            const double minBottomLayerThickness =
                atmosphereInput->minBottomLayerThicknessMeters();
            const double minTopPressureAtm =
                atmosphereInput->minTopPressureAtm();
            const double verticalWindMixingCoefficient =
                atmosphereInput->verticalWindMixingCoefficient();
            StellarParameters primaryStar{1.0, 5772.0, 1.0};
            readStellarParametersForRender(radiusInput_, temperatureInput_, primaryStar);
            std::optional<StellarParameters> secondaryStar = std::nullopt;
            if (secondStarCheckBox_->isChecked()) {
                StellarParameters secondaryCandidate{1.0, 5772.0, 1.0};
                if (readStellarParametersForRender(secondaryRadiusInput_,
                                                   secondaryTemperatureInput_,
                                                   secondaryCandidate)) {
                    secondaryStar = secondaryCandidate;
                }
            }
            const bool existingHasSeaLevel =
                existingIndex >= 0
                    ? planetComboBox_->itemData(existingIndex, kRoleHasSeaLevel).toBool()
                    : false;
            const double minDenseAtmosphereTemperatureK =
                existingIndex >= 0
                    ? planetComboBox_->itemData(existingIndex,
                                                kRoleMinDenseAtmosphereTemperature).toDouble()
                    : 0.0;
            const bool disableWaterAndClouds =
                existingIndex >= 0
                    ? planetComboBox_->itemData(existingIndex, kRoleDisableWaterAndClouds).toBool()
                    : false;
            PlanetPreset preset{name, axis, dayLength, yearLength, eccentricity, obliquity,
                                perihelionArgument, massEarths, radiusKm, materialId,
                                composition, QStringLiteral("custom"), primaryStar, secondaryStar,
                                greenhouseOpacity, manualGreenhouseOnTop, cloudAlbedo,
                                disableWaterAndClouds,
                                geothermalFlux, spinOrbitP, spinOrbitQ, tidallyLocked};
            preset.rotationModeOverride = true;
            preset.heightSeed = heightSeed;
            preset.hasSeaLevel = existingHasSeaLevel;
            preset.minBottomLayerThicknessMeters = minBottomLayerThickness;
            preset.minTopPressureAtm = minTopPressureAtm;
            preset.minDenseAtmosphereTemperatureK = minDenseAtmosphereTemperatureK;
            preset.verticalWindMixingCoefficient = verticalWindMixingCoefficient;
            if (existingIndex >= 0) {
                if (!isCustomPlanetIndex(existingIndex)) {
                    showInputError(QStringLiteral("Нельзя заменить планету из пресета."));
                    return;
                }
                const QVariant atmosphereValue = QVariant::fromValue(composition);
                planetComboBox_->setItemText(existingIndex, formatPlanetName(preset));
                planetComboBox_->setItemData(existingIndex, axis, kRoleSemiMajorAxis);
                planetComboBox_->setItemData(existingIndex, dayLength, kRoleDayLength);
                planetComboBox_->setItemData(existingIndex, yearLength, kRoleYearLengthDays);
                planetComboBox_->setItemData(existingIndex, eccentricity, kRoleEccentricity);
                planetComboBox_->setItemData(existingIndex, obliquity, kRoleObliquity);
                planetComboBox_->setItemData(existingIndex, perihelionArgument,
                                             kRolePerihelionArgument);
                planetComboBox_->setItemData(existingIndex, massEarths, kRoleMassEarths);
                planetComboBox_->setItemData(existingIndex, radiusKm, kRoleRadiusKm);
                planetComboBox_->setItemData(existingIndex, static_cast<int>(rotationMode),
                                             kRoleRotationMode);
                planetComboBox_->setItemData(existingIndex, true, kRoleRotationModeOverride);
                planetComboBox_->setItemData(existingIndex, true, kRoleIsCustom);
                planetComboBox_->setItemData(existingIndex, name, kRolePlanetName);
                planetComboBox_->setItemData(existingIndex, spinOrbitP, kRoleSpinOrbitP);
                planetComboBox_->setItemData(existingIndex, spinOrbitQ, kRoleSpinOrbitQ);
                planetComboBox_->setItemData(existingIndex, materialId, kRoleMaterialId);
                planetComboBox_->setItemData(existingIndex, atmosphereValue, kRoleAtmosphere);
                planetComboBox_->setItemData(existingIndex, preset.greenhouseOpacity,
                                             kRoleGreenhouseOpacity);
                planetComboBox_->setItemData(existingIndex,
                                             preset.manualGreenhouseOnTopOfAtmosphere,
                                             kRoleManualGreenhouseOnTopOfAtmosphere);
                const QVariant advancedRadiationValue =
                    planetComboBox_->itemData(existingIndex, kRoleAdvancedRadiationModel);
                const QVariant fallbackRadiation =
                    static_cast<int>(RadiationModelType::Layered);
                planetComboBox_->setItemData(
                    existingIndex,
                    advancedRadiationValue.isValid() ? advancedRadiationValue : fallbackRadiation,
                    kRoleAdvancedRadiationModel);
                planetComboBox_->setItemData(existingIndex, preset.cloudAlbedo, kRoleCloudAlbedo);
                planetComboBox_->setItemData(existingIndex,
                                             disableWaterAndClouds,
                                             kRoleDisableWaterAndClouds);
                planetComboBox_->setItemData(existingIndex,
                                             diurnalCoolingBiasK,
                                             kRoleDiurnalCoolingBiasK);
                planetComboBox_->setItemData(existingIndex,
                                             preset.geothermalFluxWPerM2,
                                             kRoleGeothermalFlux);
                planetComboBox_->setItemData(existingIndex,
                                             static_cast<int>(preset.heightSourceType),
                                             kRoleHeightSourceType);
                planetComboBox_->setItemData(existingIndex,
                                             preset.primaryStar.radiusInSolarRadii,
                                             kRolePrimaryStarRadius);
                planetComboBox_->setItemData(existingIndex,
                                             preset.primaryStar.temperatureKelvin,
                                             kRolePrimaryStarTemperature);
                planetComboBox_->setItemData(existingIndex,
                                             static_cast<bool>(preset.secondaryStar),
                                             kRoleHasSecondaryStar);
                if (preset.secondaryStar) {
                    planetComboBox_->setItemData(existingIndex,
                                                 preset.secondaryStar->radiusInSolarRadii,
                                                 kRoleSecondaryStarRadius);
                    planetComboBox_->setItemData(existingIndex,
                                                 preset.secondaryStar->temperatureKelvin,
                                                 kRoleSecondaryStarTemperature);
                } else {
                    planetComboBox_->setItemData(existingIndex,
                                                 QVariant(),
                                                 kRoleSecondaryStarRadius);
                    planetComboBox_->setItemData(existingIndex,
                                                 QVariant(),
                                                 kRoleSecondaryStarTemperature);
                }
                planetComboBox_->setItemData(existingIndex, preset.heightmapPath, kRoleHeightmapPath);
                planetComboBox_->setItemData(existingIndex, preset.heightmapScaleKm,
                                             kRoleHeightmapScaleKm);
                planetComboBox_->setItemData(existingIndex, preset.heightSeed, kRoleHeightSeed);
                planetComboBox_->setItemData(existingIndex, preset.useContinentsHeight,
                                             kRoleUseContinentsHeight);
                planetComboBox_->setItemData(existingIndex, preset.hasSeaLevel, kRoleHasSeaLevel);
                planetComboBox_->setItemData(existingIndex,
                                             planetComboBox_->itemData(existingIndex, kRoleFlatHeight),
                                             kRoleFlatHeight);
                planetComboBox_->setItemData(
                    existingIndex,
                    planetComboBox_->itemData(existingIndex, kRoleAtmosphereDisabled),
                    kRoleAtmosphereDisabled);
                planetComboBox_->setItemData(existingIndex, preset.starPresetId, kRoleStarPresetId);
                planetComboBox_->setItemData(existingIndex, false, kRoleHasBinaryOrbit);
                planetComboBox_->setItemData(existingIndex, 0.0, kRoleBinarySemiMajorAxis);
                planetComboBox_->setItemData(existingIndex, 0.0, kRoleBinaryPeriodDays);
                planetComboBox_->setItemData(existingIndex, 0.0, kRoleBinaryEccentricity);
                planetComboBox_->setItemData(existingIndex, 0.0, kRoleBinaryInclinationDegrees);
                planetComboBox_->setItemData(existingIndex, 0.0, kRoleBinaryArgumentPericenterA);
                planetComboBox_->setItemData(existingIndex, 0.0, kRoleBinaryArgumentPericenterB);
                planetComboBox_->setItemData(existingIndex,
                                             preset.minBottomLayerThicknessMeters,
                                             kRoleAtmosphereBottomLayerThickness);
                planetComboBox_->setItemData(existingIndex,
                                             preset.minTopPressureAtm,
                                             kRoleMinTopPressureAtm);
                planetComboBox_->setItemData(existingIndex,
                                             preset.minDenseAtmosphereTemperatureK,
                                             kRoleMinDenseAtmosphereTemperature);
                planetComboBox_->setItemData(existingIndex,
                                             preset.verticalWindMixingCoefficient,
                                             kRoleVerticalWindMixingCoefficient);
                planetComboBox_->setCurrentIndex(existingIndex);
            } else {
                addPlanetItem(preset, true);
                const int newIndex = planetComboBox_->count() - 1;
                planetComboBox_->setItemData(newIndex,
                                             diurnalCoolingBiasK,
                                             kRoleDiurnalCoolingBiasK);
                planetComboBox_->setCurrentIndex(newIndex);
            }

            updatePlanetSemiMajorAxisLabel();
            updatePlanetDayLengthLabel();
            updatePlanetYearLengthLabel();
            updatePlanetMassLabel();
            updatePlanetRadiusLabel();
            updatePlanetDerivedLabels();
            updateAtmosphereComposition();
            updatePlanetOrbitLabels();
            syncMaterialWithPlanet();
            syncRotationModeWithPlanet();
            updatePlanetActions();
            refreshStarSystemView();
            dialog.accept();
        });

        dialog.exec();
    }

    void populateMaterials() {
        for (const auto &material : surfaceMaterials()) {
            materialComboBox_->addItem(material.name, material.id);
        }
    }

    std::optional<SurfaceMaterial> currentMaterial() const {
        const QString id = materialComboBox_->currentData().toString();
        for (const auto &material : surfaceMaterials()) {
            if (material.id == id) {
                return material;
            }
        }
        return std::nullopt;
    }

    void syncMaterialWithPlanet() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }

        const QString materialId = planetComboBox_->itemData(index, kRoleMaterialId).toString();
        const int materialIndex = materialComboBox_->findData(materialId);
        if (materialIndex >= 0) {
            materialComboBox_->setCurrentIndex(materialIndex);
        }
    }

    void syncRotationModeWithPlanet() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }

        const RotationMode rotationMode = resolveRotationModeForPlanetIndex(index);
        const int modeIndex = rotationModeComboBox_->findData(static_cast<int>(rotationMode));
        if (modeIndex >= 0) {
            const QSignalBlocker blocker(rotationModeComboBox_);
            rotationModeComboBox_->setCurrentIndex(modeIndex);
        }
        updateRotationModeIllustration();
    }

    void syncSpinOrbitWithPlanet() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !spinOrbitPSpinBox_ || !spinOrbitQSpinBox_) {
            return;
        }

        const int spinOrbitP = planetComboBox_->itemData(index, kRoleSpinOrbitP).toInt();
        const int spinOrbitQ = planetComboBox_->itemData(index, kRoleSpinOrbitQ).toInt();
        const QSignalBlocker blockerP(spinOrbitPSpinBox_);
        const QSignalBlocker blockerQ(spinOrbitQSpinBox_);
        spinOrbitPSpinBox_->setValue(qMax(1, spinOrbitP));
        spinOrbitQSpinBox_->setValue(qMax(1, spinOrbitQ));
    }

    void syncHeightSeedWithPlanet() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !heightSeedSpinBox_) {
            return;
        }

        const int heightSeed = planetComboBox_->itemData(index, kRoleHeightSeed).toInt();
        const QSignalBlocker blocker(heightSeedSpinBox_);
        heightSeedSpinBox_->setValue(heightSeed);
    }

    void syncHeightmapScaleWithPlanet() {
        if (!heightmapScaleSpinBox_ || !heightmapImportButton_) {
            return;
        }

        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            const QSignalBlocker blocker(heightmapScaleSpinBox_);
            heightmapScaleSpinBox_->setValue(0.0);
            heightmapScaleSpinBox_->setEnabled(false);
            heightmapImportButton_->setEnabled(false);
            return;
        }

        const QVariant scaleValue = planetComboBox_->itemData(index, kRoleHeightmapScaleKm);
        const double scaleKm = scaleValue.isValid() ? qMax(0.0, scaleValue.toDouble()) : 0.0;
        const QSignalBlocker blocker(heightmapScaleSpinBox_);
        heightmapScaleSpinBox_->setValue(scaleKm);
        heightmapScaleSpinBox_->setEnabled(true);
        heightmapImportButton_->setEnabled(true);
    }

    void syncFlatHeightWithPlanet() {
        if (!flatHeightButton_) {
            return;
        }

        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            const QSignalBlocker blocker(flatHeightButton_);
            flatHeightButton_->setChecked(false);
            flatHeightButton_->setEnabled(false);
            updateFlatHeightButtonText(false);
            return;
        }

        const bool useFlatHeight = planetComboBox_->itemData(index, kRoleFlatHeight).toBool();
        const QSignalBlocker blocker(flatHeightButton_);
        flatHeightButton_->setChecked(useFlatHeight);
        flatHeightButton_->setEnabled(true);
        updateFlatHeightButtonText(useFlatHeight);
    }

    void syncAtmosphereDisableWithPlanet() {
        if (!disableAtmosphereButton_) {
            return;
        }

        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            const QSignalBlocker blocker(disableAtmosphereButton_);
            disableAtmosphereButton_->setChecked(false);
            disableAtmosphereButton_->setEnabled(false);
            updateAtmosphereButtonText(false);
            return;
        }

        const bool atmosphereDisabled =
            planetComboBox_->itemData(index, kRoleAtmosphereDisabled).toBool();
        const QSignalBlocker blocker(disableAtmosphereButton_);
        disableAtmosphereButton_->setChecked(atmosphereDisabled);
        disableAtmosphereButton_->setEnabled(true);
        updateAtmosphereButtonText(atmosphereDisabled);
        if (atmosphereWidget_) {
            atmosphereWidget_->setEnabled(!atmosphereDisabled);
        }
    }

    void syncCloudAlbedoWithPlanet() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !cloudAlbedoSpinBox_) {
            return;
        }

        const double cloudAlbedo = planetComboBox_->itemData(index, kRoleCloudAlbedo).toDouble();
        const QSignalBlocker blocker(cloudAlbedoSpinBox_);
        cloudAlbedoSpinBox_->setValue(cloudAlbedo);
    }

    void syncDiurnalCoolingBiasWithPlanet() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !diurnalCoolingBiasSpinBox_) {
            return;
        }

        const double diurnalCoolingBiasK =
            planetComboBox_->itemData(index, kRoleDiurnalCoolingBiasK).toDouble();
        const QSignalBlocker blocker(diurnalCoolingBiasSpinBox_);
        diurnalCoolingBiasSpinBox_->setValue(qMax(0.0, diurnalCoolingBiasK));
    }

    void syncGeothermalFluxWithPlanet() {
        if (!subsurfaceGeothermalFluxSpinBox_) {
            return;
        }
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }

        const double geothermalFlux =
            planetComboBox_->itemData(index, kRoleGeothermalFlux).toDouble();
        const QSignalBlocker blocker(subsurfaceGeothermalFluxSpinBox_);
        subsurfaceGeothermalFluxSpinBox_->setValue(geothermalFlux);
    }

    void syncManualGreenhouseOnTopWithPlanet() {
        if (!manualGreenhouseOnTopCheckBox_) {
            return;
        }
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            const QSignalBlocker blocker(manualGreenhouseOnTopCheckBox_);
            manualGreenhouseOnTopCheckBox_->setChecked(false);
            return;
        }

        const bool manualOnTop =
            planetComboBox_->itemData(index, kRoleManualGreenhouseOnTopOfAtmosphere).toBool();
        const QSignalBlocker blocker(manualGreenhouseOnTopCheckBox_);
        manualGreenhouseOnTopCheckBox_->setChecked(manualOnTop);
    }

    void syncAdvancedRadiationWithPlanet() {
        if (!advancedRadiationCheckBox_) {
            return;
        }
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            const QSignalBlocker blocker(advancedRadiationCheckBox_);
            advancedRadiationCheckBox_->setChecked(false);
            return;
        }

        const int modelTypeValue =
            planetComboBox_->itemData(index, kRoleAdvancedRadiationModel).toInt();
        const auto modelType = static_cast<RadiationModelType>(modelTypeValue);
        const bool useAdvanced = modelType == RadiationModelType::Layered;
        const QSignalBlocker blocker(advancedRadiationCheckBox_);
        advancedRadiationCheckBox_->setChecked(useAdvanced);
    }

    void syncPlanetFlatHeightWithSelection(bool useFlatHeight) {
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }
        planetComboBox_->setItemData(index, useFlatHeight, kRoleFlatHeight);
        updateFlatHeightButtonText(useFlatHeight);
    }

    void syncPlanetAtmosphereDisabledWithSelection(bool atmosphereDisabled) {
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }
        planetComboBox_->setItemData(index, atmosphereDisabled, kRoleAtmosphereDisabled);
        updateAtmosphereButtonText(atmosphereDisabled);
        if (atmosphereWidget_) {
            atmosphereWidget_->setEnabled(!atmosphereDisabled);
        }
    }

    void syncPlanetMaterialWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }
        planetComboBox_->setItemData(index, materialComboBox_->currentData(), kRoleMaterialId);
    }

    void syncPlanetRotationModeWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }
        planetComboBox_->setItemData(index, rotationModeComboBox_->currentData(), kRoleRotationMode);
        planetComboBox_->setItemData(index, true, kRoleRotationModeOverride);
    }

    void syncPlanetSpinOrbitWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !spinOrbitPSpinBox_ || !spinOrbitQSpinBox_) {
            return;
        }
        planetComboBox_->setItemData(index, spinOrbitPSpinBox_->value(), kRoleSpinOrbitP);
        planetComboBox_->setItemData(index, spinOrbitQSpinBox_->value(), kRoleSpinOrbitQ);
    }

    void syncPlanetHeightSeedWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !heightSeedSpinBox_) {
            return;
        }
        planetComboBox_->setItemData(index, heightSeedSpinBox_->value(), kRoleHeightSeed);
    }

    void syncPlanetHeightmapScaleWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !heightmapScaleSpinBox_) {
            return;
        }
        planetComboBox_->setItemData(index,
                                     heightmapScaleSpinBox_->value(),
                                     kRoleHeightmapScaleKm);
    }

    void onImportHeightmapRequested() {
        if (!planetComboBox_) {
            return;
        }
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            showInputError(QStringLiteral("Сначала выберите планету."));
            return;
        }

        const QString filePath = QFileDialog::getOpenFileName(
            this,
            QStringLiteral("Импорт heightmap"),
            QString(),
            QStringLiteral("Heightmap (*.tif *.tiff *.png *.jpg *.jpeg *.bmp);;Все файлы (*.*)"));
        if (filePath.isEmpty()) {
            return;
        }

        const double scaleKm = heightmapScaleSpinBox_ ? heightmapScaleSpinBox_->value() : 0.0;
        SurfaceHeightmap validator;
        if (!validator.loadFromFile(filePath, scaleKm)) {
            showInputError(
                QStringLiteral("Heightmap должен быть equirectangular (2:1) и 16-битным."));
            return;
        }

        QImage image(filePath);
        if (image.isNull()) {
            showInputError(QStringLiteral("Не удалось загрузить файл heightmap."));
            return;
        }

        const QString planetName =
            planetComboBox_->itemData(index, kRolePlanetName).toString();
        HeightmapStorage storage;
        constexpr int kHeightmapSubdivisionLevel = 0;
        // Сохраняем нормализованный 16-битный height-файл в папке heightmaps:
        // исходный .tif нужен только для импорта и не обязателен при старте.
        if (!storage.saveHeightmap(image, scaleKm, planetName, kHeightmapSubdivisionLevel)) {
            showInputError(QStringLiteral("Не удалось сохранить heightmap в локальном хранилище."));
            return;
        }

        const QString storedPath = storage.makeFilePath(planetName, kHeightmapSubdivisionLevel);
        planetComboBox_->setItemData(index, storedPath, kRoleHeightmapPath);
        planetComboBox_->setItemData(index, scaleKm, kRoleHeightmapScaleKm);
        planetComboBox_->setItemData(index,
                                     static_cast<int>(HeightSourceType::HeightmapEquirectangular),
                                     kRoleHeightSourceType);
        updateSurfaceGridTemperatures();
        updateTemperaturePlot();
    }

    void syncPlanetCloudAlbedoWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !cloudAlbedoSpinBox_) {
            return;
        }
        planetComboBox_->setItemData(index, cloudAlbedoSpinBox_->value(), kRoleCloudAlbedo);
    }

    void syncPlanetDiurnalCoolingBiasWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !diurnalCoolingBiasSpinBox_) {
            return;
        }
        planetComboBox_->setItemData(index,
                                     diurnalCoolingBiasSpinBox_->value(),
                                     kRoleDiurnalCoolingBiasK);
    }

    void syncPlanetManualGreenhouseOnTopWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !manualGreenhouseOnTopCheckBox_) {
            return;
        }
        planetComboBox_->setItemData(
            index,
            manualGreenhouseOnTopCheckBox_->isChecked(),
            kRoleManualGreenhouseOnTopOfAtmosphere);
    }

    void syncPlanetAdvancedRadiationWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !advancedRadiationCheckBox_) {
            return;
        }
        const auto modelType = advancedRadiationCheckBox_->isChecked()
                                   ? RadiationModelType::Layered
                                   : RadiationModelType::Fast;
        planetComboBox_->setItemData(index,
                                     static_cast<int>(modelType),
                                     kRoleAdvancedRadiationModel);
    }

    void updateRotationModeIllustration() {
        if (!modeIllustrationWidget_) {
            return;
        }
        const QVariant modeData = rotationModeComboBox_->currentData();
        const auto mode = modeData.isValid()
            ? static_cast<RotationMode>(modeData.toInt())
            : RotationMode::Normal;
        modeIllustrationWidget_->setRotationMode(mode);
    }

    RotationMode resolveRotationModeForPlanetIndex(int index) const {
        if (index < 0) {
            return RotationMode::Normal;
        }
        const double dayLengthDays = planetComboBox_->itemData(index, kRoleDayLength).toDouble();
        const double yearLengthDays =
            planetComboBox_->itemData(index, kRoleYearLengthDays).toDouble();
        const bool hasOverride =
            planetComboBox_->itemData(index, kRoleRotationModeOverride).toBool();
        const RotationMode manualMode =
            static_cast<RotationMode>(planetComboBox_->itemData(index, kRoleRotationMode).toInt());
        return resolveRotationModeFromPeriods(dayLengthDays, yearLengthDays, hasOverride, manualMode);
    }

    void clearTemperatureSegments() {
        segmentSelectorWidget_->setSegments({});
        segmentSelectorWidget_->setEnabled(false);
        temperaturePlot_->clearSeries();
    }

    void rebuildSurfaceGrid() {
        const double radiusKm = planetComboBox_->currentData(kRoleRadiusKm).toDouble();
        surfaceGrid_.setRadiusKm(radiusKm);

        const int latitudePointCount = latitudePoints();
        const int pointsPerLatitude = 6;
        const int targetPointCount = qMax(1, latitudePointCount * pointsPerLatitude);
        // Число ячеек в геодезической сетке равно 20 * 4^n, подбираем n под желаемое
        // количество точек, чтобы сохранить приблизительную плотность сетки.
        const double ratio = qMax(1.0, static_cast<double>(targetPointCount) / 20.0);
        const int subdivisionLevel = qMax(0, static_cast<int>(qRound(qLn(ratio) / qLn(4.0))));

        HeightSourceType heightSource = HeightSourceType::Procedural;
        QString heightmapPath;
        double heightmapScaleKm = 0.0;
        const QString planetName =
            planetComboBox_->currentData(kRolePlanetName).toString();
        HeightmapStorage storage;
        if (!planetName.isEmpty() &&
            storage.loadHeightmap(&heightmapPath, &heightmapScaleKm, planetName, subdivisionLevel)) {
            heightSource = HeightSourceType::HeightmapEquirectangular;
        }
        const int index = planetComboBox_->currentIndex();
        if (index >= 0) {
            planetComboBox_->setItemData(index, static_cast<int>(heightSource),
                                         kRoleHeightSourceType);
            planetComboBox_->setItemData(index, heightmapPath, kRoleHeightmapPath);
            planetComboBox_->setItemData(index, heightmapScaleKm, kRoleHeightmapScaleKm);
        }

        // Seed влияет только на HeightSourceType::Procedural.
        const quint32 heightSeed =
            planetComboBox_->currentData(kRoleHeightSeed).toUInt();
        const bool useContinentsHeight =
            planetComboBox_->currentData(kRoleUseContinentsHeight).toBool();
        const bool useFlatHeight =
            planetComboBox_->currentData(kRoleFlatHeight).toBool();
        // Уровень моря задаётся пресетом, чтобы отделение суши/океана не зависело
        // от выбранного материала (например, у Венеры остаются "океаны" даже на песке).
        const bool hasSeaLevel =
            planetComboBox_->currentData(kRoleHasSeaLevel).toBool();
        surfaceGrid_.setHeightSource(heightSource, heightmapPath, heightmapScaleKm,
                                     heightSeed, useContinentsHeight, hasSeaLevel);

        if (radiusKm <= 0.0) {
            surfaceGrid_.generateIcosahedronGrid(0);
            return;
        }
        surfaceGrid_.generateIcosahedronGrid(subdivisionLevel);

        if (useFlatHeight) {
            for (auto &point : surfaceGrid_.points()) {
                point.heightKm = 0.0;
            }
        }
    }

    std::optional<SurfacePointStateDefaults> buildSurfacePointStateDefaults() const {
        const auto material = currentMaterial();
        if (!material) {
            return std::nullopt;
        }

        const double manualGreenhouseOpacity =
            planetComboBox_->currentData(kRoleGreenhouseOpacity).toDouble();
        const double diurnalCoolingBiasK =
            planetComboBox_->currentData(kRoleDiurnalCoolingBiasK).toDouble();
        const double albedo = qBound(0.0, material->albedo, 1.0);
        // Не задаем высокий нижний порог: карта поверхности должна стартовать от физического минимума,
        // чтобы при слабой инсоляции температура могла быть значительно ниже 200 K.
        const double minTemperatureKelvin = 3.0;

        const AtmosphereComposition atmosphere = currentAtmosphereForCalculations();
        const double massEarths = planetComboBox_->currentData(kRoleMassEarths).toDouble();
        const double radiusKm = planetComboBox_->currentData(kRoleRadiusKm).toDouble();
        double atmospherePressureAtm = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            atmospherePressureAtm = atmosphere.totalPressureAtm(massEarths, radiusKm);
        }
        const bool hasAtmosphericModel = atmosphere.totalMassGigatons() > 0.0;
        const bool hasAtmospherePressure = atmospherePressureAtm > 0.0;
        const bool atmosphereEnabled = hasAtmosphericModel || hasAtmospherePressure;

        SurfacePointStateDefaults defaults;
        defaults.albedo = albedo;
        defaults.manualGreenhouseOpacity = atmosphereEnabled ? 0.0 : manualGreenhouseOpacity;
        defaults.greenhouseOpacity = defaults.manualGreenhouseOpacity;
        defaults.minTemperatureKelvin = minTemperatureKelvin;
        defaults.diurnalCoolingBiasK = qMax(0.0, diurnalCoolingBiasK);
        defaults.material = *material;
        defaults.subsurfaceSettings = buildSubsurfaceSettings();
        return defaults;
    }

    QString resolveSurfaceMaterialIdForPoint(const AtmosphereComposition &atmosphere,
                                             double massEarths,
                                             double radiusKm) const {
        const QString planetName = planetComboBox_->currentData(kRolePlanetName).toString();
        if (planetName == QStringLiteral("Церрера")) {
            return QStringLiteral("ice");
        }

        double atmospherePressureAtm = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            atmospherePressureAtm = atmosphere.totalPressureAtm(massEarths, radiusKm);
        }
        if (atmospherePressureAtm > 0.0) {
            return QStringLiteral("desert");
        }

        // Без атмосферы выбираем лунный реголит как типичный пылевой покров для безвоздушных тел.
        return QStringLiteral("regolith_moon");
    }

    SubsurfaceModelSettings buildSubsurfaceSettings() const {
        const SubsurfaceModelSettings defaultSettings;
        SubsurfaceModelSettings settings = defaultSettings;

        if (subsurfaceLayersSpinBox_ && !subsurfaceLayerCountAuto_) {
            settings.layerCount = subsurfaceLayersSpinBox_->value();
        }
        double dzFromDelta = settings.topLayerThicknessMeters;
        if (subsurfaceTopThicknessSpinBox_ && !subsurfaceTopThicknessAuto_) {
            dzFromDelta = subsurfaceTopThicknessSpinBox_->value();
            settings.topLayerThicknessMeters = dzFromDelta;
        }
        if (subsurfaceDepthSpinBox_ && !subsurfaceBottomDepthAuto_) {
            settings.bottomDepthMeters = subsurfaceDepthSpinBox_->value();
        }
        if (subsurfaceBoundaryComboBox_) {
            settings.bottomBoundary = static_cast<SubsurfaceBottomBoundaryCondition>(
                subsurfaceBoundaryComboBox_->currentData().toInt());
        }
        if (subsurfaceGeothermalFluxSpinBox_) {
            // Всегда берем актуальное значение из UI при пересчёте модели.
            settings.geothermalFluxWPerM2 = subsurfaceGeothermalFluxSpinBox_->value();
        }

        const bool useAutoLayerCount = subsurfaceLayerCountAuto_ || !subsurfaceLayersSpinBox_;
        const bool useAutoBottomDepth = subsurfaceBottomDepthAuto_ || !subsurfaceDepthSpinBox_;
        const bool useAutoTopThickness =
            subsurfaceTopThicknessAuto_ || !subsurfaceTopThicknessSpinBox_;

        double maxSolarFluxWPerM2 = 0.0;
        if (hasSolarConstant_) {
            maxSolarFluxWPerM2 = lastSolarConstant_;
            if (lastSolarConstantDistanceAU_ > 0.0) {
                const double distanceAu =
                    surfaceOrbitAnimationInitialized_
                        ? surfaceOrbitAnimation_.distanceAU()
                        : (planetComboBox_
                               ? planetComboBox_->currentData(kRoleSemiMajorAxis).toDouble()
                               : 0.0);
                if (distanceAu > 0.0) {
                    maxSolarFluxWPerM2 *=
                        std::pow(lastSolarConstantDistanceAU_ / distanceAu, 2.0);
                }
            }
        }
        if (maxSolarFluxWPerM2 <= 0.0) {
            const double distanceAu =
                planetComboBox_ ? planetComboBox_->currentData(kRoleSemiMajorAxis).toDouble() : 0.0;
            StellarParameters primary{};
            bool hasPrimary = readStellarParametersForRender(radiusInput_, temperatureInput_,
                                                             primary);
            if (!hasPrimary && planetComboBox_) {
                bool radiusOk = false;
                bool temperatureOk = false;
                const double radius =
                    planetComboBox_->currentData(kRolePrimaryStarRadius).toDouble(&radiusOk);
                const double temperature =
                    planetComboBox_->currentData(kRolePrimaryStarTemperature).toDouble(
                        &temperatureOk);
                if (radiusOk && temperatureOk && radius > 0.0 && temperature > 0.0) {
                    primary.radiusInSolarRadii = radius;
                    primary.temperatureKelvin = temperature;
                    hasPrimary = true;
                }
            }
            if (hasPrimary && distanceAu > 0.0) {
                primary.distanceInAU = distanceAu;
                maxSolarFluxWPerM2 += SolarCalculator::solarConstant(primary);
            }

            StellarParameters secondary{};
            bool hasSecondary = false;
            if (secondStarCheckBox_ && secondStarCheckBox_->isChecked()) {
                hasSecondary = readStellarParametersForRender(secondaryRadiusInput_,
                                                              secondaryTemperatureInput_,
                                                              secondary);
            } else if (planetComboBox_ &&
                       planetComboBox_->currentData(kRoleHasSecondaryStar).toBool()) {
                bool radiusOk = false;
                bool temperatureOk = false;
                const double radius =
                    planetComboBox_->currentData(kRoleSecondaryStarRadius).toDouble(&radiusOk);
                const double temperature =
                    planetComboBox_->currentData(kRoleSecondaryStarTemperature).toDouble(
                        &temperatureOk);
                if (radiusOk && temperatureOk && radius > 0.0 && temperature > 0.0) {
                    secondary.radiusInSolarRadii = radius;
                    secondary.temperatureKelvin = temperature;
                    hasSecondary = true;
                }
            }
            if (hasSecondary && distanceAu > 0.0) {
                secondary.distanceInAU = distanceAu;
                maxSolarFluxWPerM2 += SolarCalculator::solarConstant(secondary);
            }
        }

        const auto material = currentMaterial();
        const double thermalConductivity = material ? material->thermalConductivity : 1.0;
        const double density = material ? material->density : 1.0;
        const double specificHeat = material ? material->specificHeat : 1.0;
        const double safeThermalConductivity = qMax(1e-6, thermalConductivity);
        const double safeDensity = qMax(1e-6, density);
        const double safeSpecificHeat = qMax(1e-6, specificHeat);
        const double dayLengthDays =
            planetComboBox_ ? planetComboBox_->currentData(kRoleDayLength).toDouble() : 0.0;
        const double dayLengthSeconds = qMax(0.01, dayLengthDays) * 86400.0;
        const double thermalDiffusivity =
            safeThermalConductivity / (safeDensity * safeSpecificHeat);
        constexpr double kPi = 3.14159265358979323846;
        // Тепловая глубина (skin depth) для суточного периода:
        // δ = sqrt(α * P / π), где α = λ / (ρ·c), P — длительность суток.
        const double thermalSkinDepthMeters =
            std::sqrt(thermalDiffusivity * dayLengthSeconds / kPi);
        // Берём запас по глубине ~ 6δ, чтобы ночные минимумы не занижались
        // из-за слишком близкой нижней границы.
        const double targetBottomDepthMeters =
            qBound(2.5, 6.0 * thermalSkinDepthMeters, 12.0);
        const int targetLayerCount = qBound(
            24, static_cast<int>(qRound(24.0 + (targetBottomDepthMeters - 2.5) * (16.0 / 9.5))),
            40);
        // Глубина модели ~ несколько δ обеспечивает корректную амплитуду ночных температур.
        if (useAutoBottomDepth) {
            settings.bottomDepthMeters = targetBottomDepthMeters;
        }
        if (useAutoLayerCount) {
            settings.layerCount = targetLayerCount;
        }
        if (useAutoTopThickness) {
            const double uniformThickness =
                settings.bottomDepthMeters / qMax(1, settings.layerCount);
            // При авто-режиме не делаем верхний слой тоньше 0.15 м,
            // чтобы избежать слишком агрессивной дискретизации возле поверхности.
            const double minTopThicknessFromSkin =
                qBound(0.15, 0.15 * thermalSkinDepthMeters, 0.2);
            // Верхний слой должен сохранять минимальную толщину, даже если
            // равномерный шаг сетки получается меньше половины верхнего слоя.
            const double maxTopThickness =
                qMax(uniformThickness * 0.5, minTopThicknessFromSkin);
            const double desiredTopThickness =
                qMax(settings.topLayerThicknessMeters, minTopThicknessFromSkin);
            settings.topLayerThicknessMeters =
                qMin(desiredTopThickness, maxTopThickness);
            dzFromDelta = settings.topLayerThicknessMeters;
        }
        const auto adjustBottomDepthForTopThickness = [&](double requiredBottomDepthMeters) {
            if (requiredBottomDepthMeters <= settings.bottomDepthMeters) {
                return;
            }
            // Нужна подстройка, чтобы SubsurfaceGrid не переключался
            // на равномерную сетку и не снижал верхнюю толщину ниже 0.15 м.
            if (useAutoBottomDepth) {
                settings.bottomDepthMeters = requiredBottomDepthMeters;
            } else if (useAutoLayerCount) {
                settings.layerCount = qMax(
                    1,
                    static_cast<int>(std::floor(settings.bottomDepthMeters /
                                                settings.topLayerThicknessMeters)));
                if (useAutoTopThickness) {
                    const double uniformThickness =
                        settings.bottomDepthMeters / qMax(1, settings.layerCount);
                    const double minTopThicknessFromSkin =
                        qBound(0.15, 0.15 * thermalSkinDepthMeters, 0.2);
                    const double maxTopThickness =
                        qMax(uniformThickness * 0.5, minTopThicknessFromSkin);
                    const double desiredTopThickness =
                        qMax(settings.topLayerThicknessMeters, minTopThicknessFromSkin);
                    settings.topLayerThicknessMeters =
                        qMin(desiredTopThickness, maxTopThickness);
                    dzFromDelta = settings.topLayerThicknessMeters;
                }
            }
        };

        adjustBottomDepthForTopThickness(
            settings.topLayerThicknessMeters * settings.layerCount);

        // Тепловая глубина зависит от длины суток; шаг интегрирования влияет только на dz_min.
        const double timeStepSeconds = 3600.0;
        const double deltaTMaxKelvin = 2.0;
        // Ограничиваем нагрев верхнего слоя за шаг: C = ρ·c·dz, ΔT ≈ F_max·dt / C,
        // dz_min = F_max·dt / (ρ·c·ΔT_max) — физический предел на нагрев за шаг.
        const double dzMinFromHeatCapacity =
            (maxSolarFluxWPerM2 > 0.0)
                ? (maxSolarFluxWPerM2 * timeStepSeconds) /
                      (safeDensity * safeSpecificHeat * deltaTMaxKelvin)
                : 0.0;
        const double dzFromHeatCapacity = qMax(dzFromDelta, dzMinFromHeatCapacity);
        // Ограничиваем рост верхнего слоя в авто-режиме, чтобы dz_min не раздувал
        // толщину выше физически желаемого диапазона 0.15–0.2 м.
        const double maxAutoTopThicknessMeters = 0.2;
        // Если пользователь задал толщину вручную, сохраняем её — иначе верхний слой
        // может «раздуваться» из-за грубого шага по времени (1 час).
        settings.topLayerThicknessMeters = useAutoTopThickness
                                               ? qMin(maxAutoTopThicknessMeters,
                                                      dzFromHeatCapacity)
                                               : dzFromDelta;
        adjustBottomDepthForTopThickness(
            settings.topLayerThicknessMeters * settings.layerCount);
        return settings;
    }

    void updateSubsurfaceSettingsUi(const SubsurfaceModelSettings &settings) {
        // Обновляем UI только фактическими параметрами сетки и блокируем сигналы,
        // чтобы избежать рекурсивного пересчёта.
        if (subsurfaceLayersSpinBox_) {
            const QSignalBlocker blocker(subsurfaceLayersSpinBox_);
            subsurfaceLayersSpinBox_->setValue(settings.layerCount);
        }
        if (subsurfaceTopThicknessSpinBox_) {
            const QSignalBlocker blocker(subsurfaceTopThicknessSpinBox_);
            subsurfaceTopThicknessSpinBox_->setValue(settings.topLayerThicknessMeters);
        }
        if (subsurfaceDepthSpinBox_) {
            const QSignalBlocker blocker(subsurfaceDepthSpinBox_);
            subsurfaceDepthSpinBox_->setValue(settings.bottomDepthMeters);
        }
        if (subsurfaceBoundaryComboBox_) {
            const QSignalBlocker blocker(subsurfaceBoundaryComboBox_);
            const int boundaryIndex =
                subsurfaceBoundaryComboBox_->findData(static_cast<int>(settings.bottomBoundary));
            if (boundaryIndex >= 0) {
                subsurfaceBoundaryComboBox_->setCurrentIndex(boundaryIndex);
            }
        }
        if (subsurfaceGeothermalFluxSpinBox_) {
            const QSignalBlocker blocker(subsurfaceGeothermalFluxSpinBox_);
            subsurfaceGeothermalFluxSpinBox_->setValue(settings.geothermalFluxWPerM2);
        }
    }

    void refreshSubsurfaceSettingsUi() {
        if (!planetComboBox_ || planetComboBox_->currentIndex() < 0) {
            return;
        }
        if (!subsurfaceLayersSpinBox_ && !subsurfaceTopThicknessSpinBox_ &&
            !subsurfaceDepthSpinBox_) {
            return;
        }

        // Для автоматических настроек слоёв пересчитываем сетку без полного расчёта карты,
        // чтобы UI отражал реальные параметры после смены планеты.
        SubsurfaceModelSettings settings = buildSubsurfaceSettings();
        const SubsurfaceGrid resolvedGrid(settings.layerCount,
                                          settings.topLayerThicknessMeters,
                                          settings.bottomDepthMeters);
        settings.layerCount = resolvedGrid.layerCount();
        settings.bottomDepthMeters = resolvedGrid.bottomDepthMeters();
        const auto &resolvedThicknesses = resolvedGrid.layerThicknessesMeters();
        if (!resolvedThicknesses.isEmpty()) {
            settings.topLayerThicknessMeters = resolvedThicknesses.front();
        }
        updateSubsurfaceSettingsUi(settings);
    }

    void applySurfaceGridToViews() {
        if (surfaceMapWidget_) {
            surfaceMapWidget_->setGrid(&surfaceGrid_);
        }
        if (surfaceGlobeWidget_) {
            surfaceGlobeWidget_->setGrid(&surfaceGrid_);
        }
        updateSurfacePointStatusDialog();
    }

    void applySurfaceTemperatureRangeToViews(double minTemperature, double maxTemperature) {
        if (surfaceMapWidget_) {
            surfaceMapWidget_->setTemperatureRange(minTemperature, maxTemperature);
        }
        if (surfaceGlobeWidget_) {
            surfaceGlobeWidget_->setTemperatureRange(minTemperature, maxTemperature);
        }
    }

    void applySurfaceWindRangeToViews(double minWindSpeed, double maxWindSpeed) {
        if (surfaceMapWidget_) {
            surfaceMapWidget_->setWindRange(minWindSpeed, maxWindSpeed);
        }
        if (surfaceGlobeWidget_) {
            surfaceGlobeWidget_->setWindRange(minWindSpeed, maxWindSpeed);
        }
    }

    void applySurfacePressureRangeToViews(double minPressureAtm, double maxPressureAtm) {
        if (surfaceMapWidget_) {
            surfaceMapWidget_->setPressureRange(minPressureAtm, maxPressureAtm);
        }
        if (surfaceGlobeWidget_) {
            surfaceGlobeWidget_->setPressureRange(minPressureAtm, maxPressureAtm);
        }
    }

    void applySurfaceMapMode(SurfaceMapMode mode) {
        surfaceMapMode_ = mode;
        if (surfaceMapWidget_) {
            surfaceMapWidget_->setMapMode(mode);
        }
        if (surfaceGlobeWidget_) {
            surfaceGlobeWidget_->setMapMode(mode);
        }
        if (mode == SurfaceMapMode::Temperature && hasSurfaceTemperatureRange_) {
            applySurfaceTemperatureRangeToViews(surfaceMinTemperatureK_, surfaceMaxTemperatureK_);
        } else if (mode == SurfaceMapMode::AirTemperature && hasSurfaceAirTemperatureRange_) {
            applySurfaceTemperatureRangeToViews(surfaceMinAirTemperatureK_,
                                                surfaceMaxAirTemperatureK_);
        }
        refreshSurfaceLegend();
    }

    bool readStellarParametersForRender(QLineEdit *radiusInput,
                                        QLineEdit *temperatureInput,
                                        StellarParameters &parameters) const {
        if (!radiusInput || !temperatureInput) {
            return false;
        }
        bool ok = false;
        const double radius = radiusInput->text().toDouble(&ok);
        if (!ok || radius <= 0.0) {
            return false;
        }
        const double temperature = temperatureInput->text().toDouble(&ok);
        if (!ok || temperature <= 0.0) {
            return false;
        }
        parameters.radiusInSolarRadii = radius;
        parameters.temperatureKelvin = temperature;
        return true;
    }

    void updateSurfaceStarRendering(double declinationDegrees,
                                    double elapsedTimeDays,
                                    double rotationDayLengthDays,
                                    RotationMode rotationMode,
                                    double orbitalPhaseRadians,
                                    int spinOrbitP,
                                    int spinOrbitQ,
                                    bool isRetrograde,
                                    double distanceAu) {
        if (!surfaceGlobeWidget_) {
            return;
        }

        const QVector3D starDirection =
            resolveStarDirection(declinationDegrees,
                                 elapsedTimeDays,
                                 rotationDayLengthDays,
                                 rotationMode,
                                 orbitalPhaseRadians,
                                 spinOrbitP,
                                 spinOrbitQ,
                                 isRetrograde);
        surfaceGlobeWidget_->setStarDirection(starDirection);

        StellarParameters primary{};
        bool hasPrimary = readStellarParametersForRender(radiusInput_, temperatureInput_, primary);
        if (!hasPrimary && planetComboBox_ && planetComboBox_->currentIndex() >= 0) {
            // Запасной путь: когда поля ввода пусты/некорректны, используем параметры звезды
            // из выбранного пресета планеты.
            bool radiusOk = false;
            bool temperatureOk = false;
            const double radius =
                planetComboBox_->currentData(kRolePrimaryStarRadius).toDouble(&radiusOk);
            const double temperature =
                planetComboBox_->currentData(kRolePrimaryStarTemperature).toDouble(&temperatureOk);
            if (radiusOk && temperatureOk && radius > 0.0 && temperature > 0.0) {
                primary.radiusInSolarRadii = radius;
                primary.temperatureKelvin = temperature;
                hasPrimary = true;
            }
        }
        if (hasPrimary) {
            // Для рендера используем основную звезду: направление для вторичной не задано.
            surfaceGlobeWidget_->setStarTemperature(primary.temperatureKelvin);
            surfaceGlobeWidget_->setStarRadiusSolar(primary.radiusInSolarRadii);
        }

        if (distanceAu > 0.0) {
            surfaceGlobeWidget_->setStarDistanceAu(distanceAu);
        }
    }

    SurfaceFluxBreakdown computeSurfaceFluxBreakdown() {
        SurfaceFluxBreakdown breakdown;
        if (!planetComboBox_ || planetComboBox_->currentIndex() < 0) {
            return breakdown;
        }

        syncSurfaceOrbitAnimationToTime();
        const double orbitalPhaseRadians = surfaceOrbitAnimation_.orbitalPhaseRadians();
        const double planetDistanceAu = surfaceOrbitAnimation_.distanceAU();
        if (planetDistanceAu <= 0.0) {
            return breakdown;
        }

        const QPointF planetPositionAu = orbitPointFromPhase(planetDistanceAu,
                                                             orbitalPhaseRadians);

        StellarParameters primary{};
        bool hasPrimary = readStellarParametersForRender(radiusInput_, temperatureInput_, primary);
        if (!hasPrimary && planetComboBox_->currentIndex() >= 0) {
            bool radiusOk = false;
            bool temperatureOk = false;
            const double radius =
                planetComboBox_->currentData(kRolePrimaryStarRadius).toDouble(&radiusOk);
            const double temperature =
                planetComboBox_->currentData(kRolePrimaryStarTemperature).toDouble(&temperatureOk);
            if (radiusOk && temperatureOk && radius > 0.0 && temperature > 0.0) {
                primary.radiusInSolarRadii = radius;
                primary.temperatureKelvin = temperature;
                hasPrimary = true;
            }
        }

        StellarParameters secondary{};
        bool hasSecondary = false;
        if (secondStarCheckBox_ && secondStarCheckBox_->isChecked()) {
            hasSecondary = readStellarParametersForRender(secondaryRadiusInput_,
                                                          secondaryTemperatureInput_,
                                                          secondary);
        }
        if (!hasSecondary && planetComboBox_->currentIndex() >= 0) {
            bool radiusOk = false;
            bool temperatureOk = false;
            const double radius =
                planetComboBox_->currentData(kRoleSecondaryStarRadius).toDouble(&radiusOk);
            const double temperature =
                planetComboBox_->currentData(kRoleSecondaryStarTemperature).toDouble(&temperatureOk);
            if (radiusOk && temperatureOk && radius > 0.0 && temperature > 0.0) {
                secondary.radiusInSolarRadii = radius;
                secondary.temperatureKelvin = temperature;
                hasSecondary = true;
            }
        }

        QPointF primaryPositionAu(0.0, 0.0);
        QPointF secondaryPositionAu(0.0, 0.0);
        const bool hasBinaryOrbit = planetComboBox_->currentData(kRoleHasBinaryOrbit).toBool();
        if (hasBinaryOrbit) {
            const double semiMajorAxisAu =
                planetComboBox_->currentData(kRoleBinarySemiMajorAxis).toDouble();
            const double periodDays =
                planetComboBox_->currentData(kRoleBinaryPeriodDays).toDouble();
            const double eccentricity =
                planetComboBox_->currentData(kRoleBinaryEccentricity).toDouble();
            const double inclinationDegrees =
                planetComboBox_->currentData(kRoleBinaryInclinationDegrees).toDouble();
            const double argumentPericenterDegreesA =
                planetComboBox_->currentData(kRoleBinaryArgumentPericenterA).toDouble();
            const double argumentPericenterDegreesB =
                planetComboBox_->currentData(kRoleBinaryArgumentPericenterB).toDouble();
            if (semiMajorAxisAu > 0.0 && periodDays > 0.0) {
                // Положение звёзд вокруг барицентра считаем по тем же элементам,
                // что используются в StarSystemTopView.
                const double elapsedDays = resolveStarSystemElapsedDays();
                const double meanAnomaly = 2.0 * M_PI * (elapsedDays / periodDays);
                const double trueAnomaly =
                    trueAnomalyFromMeanAnomalyRadians(meanAnomaly, eccentricity);
                const double inclinationRadians = qDegreesToRadians(inclinationDegrees);
                primaryPositionAu = orbitPointAuInclined(
                    semiMajorAxisAu,
                    eccentricity,
                    qDegreesToRadians(argumentPericenterDegreesA),
                    inclinationRadians,
                    trueAnomaly);
                secondaryPositionAu = orbitPointAuInclined(
                    semiMajorAxisAu,
                    eccentricity,
                    qDegreesToRadians(argumentPericenterDegreesB),
                    inclinationRadians,
                    trueAnomaly);
            }
        }

        // Поток рассчитываем по текущему расстоянию до каждой звезды:
        // F = L / (4πr^2), где r — актуальная дистанция «звезда–планета».
        if (hasPrimary) {
            breakdown.primaryDistanceAu =
                std::hypot(planetPositionAu.x() - primaryPositionAu.x(),
                           planetPositionAu.y() - primaryPositionAu.y());
            breakdown.primaryFluxWPerM2 = SolarCalculator::stellarFluxAtDistance(
                primary.radiusInSolarRadii,
                primary.temperatureKelvin,
                breakdown.primaryDistanceAu);
            breakdown.totalFluxWPerM2 += breakdown.primaryFluxWPerM2;
            breakdown.hasPrimary = true;
        }
        if (hasSecondary) {
            breakdown.secondaryDistanceAu =
                std::hypot(planetPositionAu.x() - secondaryPositionAu.x(),
                           planetPositionAu.y() - secondaryPositionAu.y());
            breakdown.secondaryFluxWPerM2 = SolarCalculator::stellarFluxAtDistance(
                secondary.radiusInSolarRadii,
                secondary.temperatureKelvin,
                breakdown.secondaryDistanceAu);
            breakdown.totalFluxWPerM2 += breakdown.secondaryFluxWPerM2;
            breakdown.hasSecondary = true;
        }

        return breakdown;
    }

    QString buildSurfaceFluxTooltip(const SurfaceFluxBreakdown &breakdown) const {
        if (breakdown.totalFluxWPerM2 <= 0.0) {
            return QStringLiteral("Солнечная постоянная недоступна.");
        }

        QStringList lines;
        lines << QStringLiteral("Суммарный поток: %1 Вт/м²")
                     .arg(breakdown.totalFluxWPerM2, 0, 'f', 1);
        if (breakdown.hasPrimary) {
            lines << QStringLiteral("Первая звезда: %1 Вт/м² (r=%2 а.е.)")
                         .arg(breakdown.primaryFluxWPerM2, 0, 'f', 1)
                         .arg(breakdown.primaryDistanceAu, 0, 'f', 3);
        }
        if (breakdown.hasSecondary) {
            lines << QStringLiteral("Вторая звезда: %1 Вт/м² (r=%2 а.е.)")
                         .arg(breakdown.secondaryFluxWPerM2, 0, 'f', 1)
                         .arg(breakdown.secondaryDistanceAu, 0, 'f', 3);
        }
        lines << QStringLiteral("Расчёт выполнен по текущим расстояниям.");
        return lines.join(QStringLiteral("\n"));
    }

    void updateSurfaceSolarConstantLabel(const SurfaceFluxBreakdown &breakdown) {
        if (!surfaceSolarConstantLabel_) {
            return;
        }
        if (breakdown.totalFluxWPerM2 <= 0.0) {
            surfaceSolarConstantLabel_->setText(QStringLiteral("S0: —"));
            surfaceSolarConstantLabel_->setToolTip(buildSurfaceFluxTooltip(breakdown));
            return;
        }
        surfaceSolarConstantLabel_->setText(
            QStringLiteral("S0: %1 Вт/м²").arg(breakdown.totalFluxWPerM2, 0, 'f', 1));
        surfaceSolarConstantLabel_->setToolTip(buildSurfaceFluxTooltip(breakdown));
    }

    void updateSurfaceOrbitLighting() {
        if (!surfaceGlobeWidget_ || !planetComboBox_ || planetComboBox_->currentIndex() < 0) {
            updateSurfaceSolarConstantLabel(SurfaceFluxBreakdown{});
            return;
        }

        syncSurfaceOrbitAnimationToTime();
        const double dayLengthDays = planetComboBox_->currentData(kRoleDayLength).toDouble();
        const double yearLengthDays =
            planetComboBox_->currentData(kRoleYearLengthDays).toDouble();
        const bool isRetrograde =
            isRetrogradeRotation(planetComboBox_->currentData(kRoleObliquity).toDouble());
        const double siderealDayLengthDays =
            resolveSiderealDayLengthDays(dayLengthDays, yearLengthDays, isRetrograde);
        const int solarStepsPerDay = qMax(1, surfaceTime_.stepsPerDay);
        const double elapsedTimeDays =
            (surfaceTime_.surfaceElapsedDays() + 0.5 / static_cast<double>(solarStepsPerDay)) *
            dayLengthDays;
        const double declinationDegrees = surfaceOrbitAnimation_.declinationDegrees();
        const double orbitalPhaseRadians = surfaceOrbitAnimation_.orbitalPhaseRadians();
        const RotationMode rotationMode =
            resolveRotationModeForPlanetIndex(planetComboBox_->currentIndex());
        const int spinOrbitP = planetComboBox_->currentData(kRoleSpinOrbitP).toInt();
        const int spinOrbitQ = planetComboBox_->currentData(kRoleSpinOrbitQ).toInt();

        const SurfaceFluxBreakdown breakdown = computeSurfaceFluxBreakdown();
        const double renderDistanceAu =
            breakdown.hasPrimary && breakdown.primaryDistanceAu > 0.0
                ? breakdown.primaryDistanceAu
                : surfaceOrbitAnimation_.distanceAU();
        // Направление на звезду совпадает с нормалью подсолнечной точки.
        updateSurfaceStarRendering(declinationDegrees,
                                   elapsedTimeDays,
                                   siderealDayLengthDays,
                                   rotationMode,
                                   orbitalPhaseRadians,
                                   spinOrbitP,
                                   spinOrbitQ,
                                   isRetrograde,
                                   renderDistanceAu);

        updateSurfaceSolarConstantLabel(breakdown);
    }

    void updateSurfaceWindField(const AtmosphereComposition &atmosphere,
                                double atmospherePressureAtm,
                                double dayLengthDays,
                                double surfaceGravity) {
        Q_UNUSED(surfaceGravity)
        if (surfaceGrid_.points().isEmpty()) {
            updateSurfaceWindLegend(false, 0.0, 0.0);
            return;
        }

        constexpr double kMinAtmospherePressureAtm = 1e-6;
        if (atmospherePressureAtm <= kMinAtmospherePressureAtm ||
            atmosphere.totalMassGigatons() <= 0.0) {
            // Без атмосферы нет ветра и связанного переноса энергии по поверхности.
            for (auto &point : surfaceGrid_.points()) {
                point.windEastMps = 0.0;
                point.windNorthMps = 0.0;
                point.windSpeedMps = 0.0;
            }
            updateSurfaceWindLegend(false, 0.0, 0.0);
            return;
        }

        QVector<double> pressuresAtm;
        QVector<double> temperatures;
        pressuresAtm.reserve(surfaceGrid_.points().size());
        temperatures.reserve(surfaceGrid_.points().size());
        for (auto &point : surfaceGrid_.points()) {
            pressuresAtm.push_back(qMax(0.0, point.pressureAtm));
            temperatures.push_back(point.temperatureK);
        }

        WindFieldModel windModel;
        const QVector<WindVector> wind =
            windModel.buildField(surfaceGrid_,
                                 pressuresAtm,
                                 temperatures,
                                 qMax(0.0, dayLengthDays) * 86400.0,
                                 2);
        if (wind.size() != surfaceGrid_.points().size()) {
            for (auto &point : surfaceGrid_.points()) {
                point.windEastMps = 0.0;
                point.windNorthMps = 0.0;
                point.windSpeedMps = 0.0;
            }
            updateSurfaceWindLegend(false, 0.0, 0.0);
            return;
        }

        double minWind = std::numeric_limits<double>::max();
        double maxWind = std::numeric_limits<double>::lowest();
        for (int i = 0; i < wind.size(); ++i) {
            auto &point = surfaceGrid_.points()[i];
            point.windEastMps = wind.at(i).eastMps;
            point.windNorthMps = wind.at(i).northMps;
            point.windSpeedMps = std::hypot(point.windEastMps, point.windNorthMps);
            minWind = qMin(minWind, point.windSpeedMps);
            maxWind = qMax(maxWind, point.windSpeedMps);
        }

        if (minWind <= maxWind) {
            applySurfaceWindRangeToViews(minWind, maxWind);
            updateSurfaceWindLegend(true, minWind, maxWind);
        } else {
            updateSurfaceWindLegend(false, 0.0, 0.0);
        }
    }

    bool isSurfaceTabActive() const {
        return rightTabs_ && surfaceMapContainer_ &&
            rightTabs_->currentWidget() == surfaceMapContainer_;
    }

    void setSurfaceGridCalculationRunning(bool running) {
        surfaceGridCalculationRunning_ = running;
        updateSurfaceCalculationIndicator();
    }

    void updateSurfaceCalculationIndicator() {
        const bool showIndicator = surfaceGridCalculationRunning_ && isSurfaceTabActive();
        if (surfaceCalculationLabel_) {
            surfaceCalculationLabel_->setVisible(showIndicator);
        }
        if (surfaceCalculationProgress_) {
            surfaceCalculationProgress_->setVisible(showIndicator);
        }
    }

    void scheduleSurfaceGridTemperatureUpdate(bool debounce) {
        if (debounce) {
            if (!surfaceGridDebounceTimer_) {
                surfaceGridDebounceTimer_ = new QTimer(this);
                surfaceGridDebounceTimer_->setSingleShot(true);
                connect(surfaceGridDebounceTimer_, &QTimer::timeout, this,
                        [this]() { startSurfaceGridTemperatureUpdate(); });
            }
            surfaceGridDebounceTimer_->start(150);
            return;
        }

        if (surfaceGridDebounceTimer_) {
            surfaceGridDebounceTimer_->stop();
        }
        startSurfaceGridTemperatureUpdate();
    }

    void cancelSurfaceGridTemperatureUpdate() {
        if (surfaceGridCancelFlag_) {
            surfaceGridCancelFlag_->store(true);
        }
        if (surfaceGridWatcher_) {
            auto future = surfaceGridWatcher_->future();
            if (future.isRunning()) {
                future.cancel();
            }
        }
    }

    static SurfaceGridComputationResult computeSurfaceGridTemperatures(
        const SurfaceGridComputationInput &input,
        const std::shared_ptr<std::atomic_bool> &cancelFlag,
        const std::shared_ptr<std::atomic<int>> &requestIdCounter,
        int requestId) {
        SurfaceGridComputationResult result;
        result.grid = input.grid;
        result.declinationDegrees = input.declinationDegrees;
        result.orbitalPhaseRadians = input.orbitalPhaseRadians;
        result.stepsPerDay = input.stepsPerDay;
        result.rotationMode = input.rotationMode;
        result.spinOrbitP = input.spinOrbitP;
        result.spinOrbitQ = input.spinOrbitQ;
        result.isRetrograde = input.isRetrograde;
        result.dayLengthDays = input.dayLengthDays;
        result.yearLengthDays = input.yearLengthDays;
        result.elapsedDays = input.elapsedDays;

        const auto shouldCancel = [&cancelFlag, requestIdCounter, requestId]() {
            if (cancelFlag && cancelFlag->load()) {
                return true;
            }
            return requestIdCounter && requestIdCounter->load() != requestId;
        };

        if (result.grid.points().isEmpty()) {
            return result;
        }

        double maxHeightKm = 0.0;
        for (const auto &point : result.grid.points()) {
            maxHeightKm = qMax(maxHeightKm, point.heightKm);
        }

        for (int i = 0; i < result.grid.points().size(); ++i) {
            if ((i % 64) == 0 && shouldCancel()) {
                result.cancelled = true;
                return result;
            }
            auto &point = result.grid.points()[i];
            if (maxHeightKm > 0.0 && point.heightKm >= 0.75 * maxHeightKm) {
                // Высотный порог для скальных пород: доля от максимума дает масштабируемый
                // критерий для разных планетных рельефов.
                point.materialId = QStringLiteral("rocky");
            } else {
                point.materialId = input.baseMaterialId;
            }
        }

        const double dayLengthDays = input.dayLengthDays;
        double atmospherePressureAtm = 0.0;
        if (input.massEarths > 0.0 && input.radiusKm > 0.0) {
            atmospherePressureAtm = input.atmosphere.totalPressureAtm(input.massEarths,
                                                                      input.radiusKm);
        }
        const bool useAtmosphericModel = input.atmosphere.totalMassGigatons() > 0.0;
        const double localSeaLevelPressureAtm =
            calculateCellPressureAtmFromKg(input.atmosphere.totalMassKg(),
                                           input.massEarths,
                                           input.radiusKm,
                                           result.grid.pointAreaKm2());
        double surfaceGravity = 0.0;
        if (input.massEarths > 0.0 && input.radiusKm > 0.0) {
            const double radiusMeters = input.radiusKm * 1000.0;
            const double planetMassKg = input.massEarths * kEarthMassKg;
            surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);
        }
        const double gravity = (surfaceGravity > 0.0) ? surfaceGravity : 9.80665;
        const double atmosphereMassGt = input.atmosphere.totalMassGigatons();
        const double co2MassGt = input.atmosphere.massGigatons(QStringLiteral("co2"));
        const double co2Share = atmosphereMassGt > 0.0 ? co2MassGt / atmosphereMassGt : 0.0;

        SurfaceTileTemperatureDefaults tileDefaults;
        tileDefaults.minTemperatureKelvin = input.stateDefaults.minTemperatureKelvin;
        tileDefaults.greenhouseOpacity = input.stateDefaults.greenhouseOpacity;
        tileDefaults.radiationModelType = input.stateDefaults.radiationModelType;
        tileDefaults.diurnalCoolingBiasK = input.stateDefaults.diurnalCoolingBiasK;
        tileDefaults.subsurfaceSettings = input.stateDefaults.subsurfaceSettings;

        SurfaceTileTemperatureSettings tileSettings;
        tileSettings.segmentSolarConstant = input.segmentSolarConstant;
        tileSettings.cloudAlbedo = input.cloudAlbedo;
        tileSettings.hasSolarConstant = input.hasSolarConstant;

        SurfaceTileTemperatureCalculator tileCalculator;
        SurfaceTileTemperatureResult tileResult =
            tileCalculator.initializeSurface(result.grid,
                                             tileSettings,
                                             tileDefaults,
                                             input.materialsById,
                                             input.stateDefaults.material);
        result.resolvedSubsurfaceSettings = tileResult.resolvedSubsurfaceSettings;

        result.hasTemperatureRange = tileResult.hasTemperatureRange;
        result.minTemperatureK = tileResult.minSurfaceTemperatureK;
        result.maxTemperatureK = tileResult.maxSurfaceTemperatureK;
        double minAirTemperature = std::numeric_limits<double>::max();
        double maxAirTemperature = std::numeric_limits<double>::lowest();
        QVector<double> localInsolations = tileResult.blendedInsolations;
        QVector<double> baselineAirTemperatures = tileResult.baselineAirTemperatures;
        const double timeStepSeconds = 3600.0;
        // Коэффициент прохождения коротковолнового излучения через облака.
        double cloudShortwaveTransmission = 1.0 - input.cloudAlbedo;
        if (input.cloudAlbedo > 0.7) {
            // Для плотных сернокислотных облаков дополнительно ослабляем поток к поверхности.
            cloudShortwaveTransmission *= 0.2;
        }
        cloudShortwaveTransmission = qBound(0.0, cloudShortwaveTransmission, 1.0);
        if (localInsolations.size() != result.grid.points().size()) {
            localInsolations.clear();
            localInsolations.resize(result.grid.points().size());
        }
        if (baselineAirTemperatures.size() != result.grid.points().size()) {
            baselineAirTemperatures.clear();
            baselineAirTemperatures.resize(result.grid.points().size());
        }

        // Радиационный профиль должен опираться на средний ТОА-поток, а не на локальную
        // инсоляцию: вертикальная структура атмосферы задаётся глобальным балансом энергии
        // (суточно-усреднённый приток), а локальная инсоляция влияет только на SW-поток внизу.
        const double planetaryAlbedo = input.cloudAlbedo;
        const double meanToaFlux =
            input.segmentSolarConstant * qMax(0.0, 1.0 - planetaryAlbedo) / 4.0;
        const double effectiveTemperatureKelvin =
            std::pow(qMax(0.0, meanToaFlux) / kStefanBoltzmannConstant, 0.25);

        double minPressureAtm = std::numeric_limits<double>::max();
        double maxPressureAtm = std::numeric_limits<double>::lowest();
        const int logPointIndex = qBound(0, input.logPointIndex, result.grid.points().size() - 1);
        for (int i = 0; i < result.grid.points().size(); ++i) {
            if ((i % 64) == 0 && shouldCancel()) {
                result.cancelled = true;
                return result;
            }
            auto &point = result.grid.points()[i];
            const double baselineAirTemperature =
                (i < baselineAirTemperatures.size()) ? baselineAirTemperatures.at(i) : 0.0;
            const double surfaceAirTemperatureK =
                (point.airTemperatureK > 0.0)
                    ? point.airTemperatureK
                    : ((baselineAirTemperature > 0.0) ? baselineAirTemperature : point.temperatureK);
            double airColumnTemperatureK = surfaceAirTemperatureK;
            bool usedProfileTemperature = false;
            if (airColumnTemperatureK > 0.0) {
                AtmosphericProfileInitializer profileInitializer(input.atmosphere);
                AtmosphericProfileInitializer::Settings profileSettings;
                profileSettings.surfaceTemperatureKelvin = airColumnTemperatureK;
                profileSettings.surfacePressureAtm = localSeaLevelPressureAtm;
                profileSettings.gravityMps2 = gravity;
                profileSettings.layerCount = 12;
                profileSettings.useDryAdiabatic = true;
                profileSettings.minTopPressureFraction = 0.01;
                const QVector<AtmosphericProfileLayer> profileLayers =
                    profileInitializer.buildProfile(profileSettings);
                double columnMass = 0.0;
                double columnTemperatureSum = 0.0;
                for (const auto &layer : profileLayers) {
                    const double layerMass = layer.densityKgPerM3 * layer.thicknessMeters;
                    if (layerMass > 0.0 && layer.temperatureKelvin > 0.0) {
                        columnMass += layerMass;
                        columnTemperatureSum += layer.temperatureKelvin * layerMass;
                    }
                }
                if (columnMass > 0.0) {
                    airColumnTemperatureK = columnTemperatureSum / columnMass;
                    usedProfileTemperature = true;
                }
            }
            if (!usedProfileTemperature && airColumnTemperatureK > 0.0) {
                // Fallback: стандартный градиент 6.5 K/км, берём среднюю температуру
                // столба на половине высоты.
                constexpr double kStandardLapseRateKPerM = 0.0065;
                const double halfHeightMeters = point.heightKm * 1000.0 * 0.5;
                airColumnTemperatureK =
                    qMax(1.0, airColumnTemperatureK - kStandardLapseRateKPerM * halfHeightMeters);
            }
            // Инициализируем поверхностное давление (на уровне рельефа), чтобы дальше
            // переносить его ветром без пересчёта из всей атмосферы каждый тик.
            // Это нужно и в fallback-ветке, чтобы карта давления/ветра не оставалась с нулями.
            // Масса газового столба в ячейке пропорциональна её площади, поэтому P = (m_cell * g) / area_cell.
            // Для масштаба высоты нужен воздушный профиль, а не температура грунта.
            // Ожидаемый контроль: для Земли ~0.7 атм на 3 км при 280–290 K.
            const double pressureAtm =
                AtmosphericPressureModel::pressureAtHeightAtm(localSeaLevelPressureAtm,
                                                              point.heightKm * 1000.0,
                                                              airColumnTemperatureK,
                                                              input.atmosphere,
                                                              gravity);
            // Храним поверхностное давление в точке для последующего переноса по ветру.
            point.pressureAtm = qMax(0.0, pressureAtm);
            minPressureAtm = qMin(minPressureAtm, point.pressureAtm);
            maxPressureAtm = qMax(maxPressureAtm, point.pressureAtm);
            const double localInsolation =
                (i < localInsolations.size()) ? localInsolations.at(i) : 0.0;
            const auto materialIt = input.materialsById.constFind(point.materialId);
            const SurfaceMaterial material =
                materialIt != input.materialsById.cend()
                    ? *materialIt
                    : input.stateDefaults.material;
            const bool logDetails = i == logPointIndex;
            const double greenhouseBaseTemperature =
                (input.radiationModelType == RadiationModelType::Layered &&
                 point.airTemperatureK > 0.0)
                    ? point.airTemperatureK
                    : point.state.temperatureKelvin();
            const double localGreenhouseOpacity =
                computeLocalGreenhouseOpacity(input.atmosphere,
                                              material,
                                              input.cloudAlbedo,
                                              input.disableWaterAndClouds,
                                              point.pressureAtm,
                                              input.radiusKm,
                                              gravity,
                                              localInsolation,
                                              meanToaFlux,
                                              greenhouseBaseTemperature,
                                              input.manualGreenhouseOpacity,
                                              input.minDenseAtmosphereTemperatureK,
                                              useAtmosphericModel,
                                              input.radiationModelType,
                                              input.manualGreenhouseOnTop,
                                              logDetails);
            point.state.setGreenhouseOpacity(localGreenhouseOpacity);
            // Воздух интегрируется по времени, иначе он будет “сбрасываться” каждый тик.
            // Поэтому используем предыдущее значение, а базу только для первичной инициализации.
            const double initialAirTemperature =
                (point.airTemperatureK > 0.0)
                    ? point.airTemperatureK
                    : ((i < baselineAirTemperatures.size())
                           ? baselineAirTemperatures.at(i)
                           : point.temperatureK);
            double resolvedInitialAirTemperature = initialAirTemperature;
            if (point.pressureAtm >= 0.1) {
                // Для плотных атмосфер t_eff недостаточна (пример: Венера) —
                // верхняя граница излучает холодно, но нижние слои остаются
                // тёплыми из-за давления и состава. Поэтому используем базу,
                // зависящую от давления и оптической плотности.
                const double tMinDense = denseAtmosphereMinTemperatureKelvin(
                    effectiveTemperatureKelvin,
                    point.pressureAtm,
                    localGreenhouseOpacity,
                    co2Share,
                    input.minDenseAtmosphereTemperatureK);
                resolvedInitialAirTemperature =
                    qMax(resolvedInitialAirTemperature, tMinDense);
            }
            if (logDetails) {
                qCInfo(solarRadiationLog) << "Radiation inputs (init)"
                                          << "index=" << i
                                          << "localInsolation=" << localInsolation
                                          << "planetaryAlbedo=" << planetaryAlbedo
                                          << "meanToaFlux=" << meanToaFlux
                                          << "effectiveTemperatureKelvin="
                                          << effectiveTemperatureKelvin;
            }
            // Стабилизация для плотных атмосфер (см. computeLocalGreenhouseOpacity):
            // t_eff не даёт реалистичной нижней границы при больших давлениях (Венера),
            // поэтому минимальная база зависит от давления и состава.
            const double profileBaseTemperatureKelvin =
                (point.pressureAtm >= 0.1)
                    ? qMax(point.state.temperatureKelvin(),
                           denseAtmosphereMinTemperatureKelvin(effectiveTemperatureKelvin,
                                                               point.pressureAtm,
                                                               localGreenhouseOpacity,
                                                               co2Share,
                                                               input.minDenseAtmosphereTemperatureK))
                    : point.state.temperatureKelvin();
            const auto radiationModel =
                makeRadiationModel(input.atmosphere,
                                   point.pressureAtm,
                                   profileBaseTemperatureKelvin,
                                   effectiveTemperatureKelvin,
                                   gravity,
                                   input.radiationModelType);
            // Поток у поверхности: S_surf = S_local * T_cloud * T_atm.
            const double surfaceShortwaveFlux =
                localInsolation * cloudShortwaveTransmission *
                radiationModel->incomingTransmission();
            // Храним текущий солнечный поток у поверхности в точке для подсказки тайла.
            point.solarFluxWPerM2 = surfaceShortwaveFlux;
            point.shortwaveSurfaceWPerM2 = surfaceShortwaveFlux;
            point.airTemperatureK =
                resolveAirTemperatureKelvin(point.state,
                                            point.pressureAtm,
                                            gravity,
                                            resolvedInitialAirTemperature,
                                            localInsolation,
                                            cloudShortwaveTransmission,
                                            input.radiationModelType,
                                            *radiationModel,
                                            useAtmosphericModel,
                                            input.manualGreenhouseOnTop,
                                            timeStepSeconds);
            const double emittedFlux = point.state.surfaceEmittedFlux();
            const double surfaceLongwaveUp =
                resolveSurfaceLongwaveUpFlux(*radiationModel,
                                             input.radiationModelType,
                                             emittedFlux);
            const double airTemperatureForFlux =
                (point.airTemperatureK > 0.0) ? point.airTemperatureK : point.temperatureK;
            const LongwaveFluxes longwaveFluxes =
                computeLongwaveFluxes(*radiationModel,
                                      input.radiationModelType,
                                      emittedFlux,
                                      airTemperatureForFlux,
                                      point.state.greenhouseOpacity(),
                                      useAtmosphericModel,
                                      input.manualGreenhouseOnTop);
            point.longwaveUpWPerM2 = surfaceLongwaveUp;
            point.longwaveDownWPerM2 = longwaveFluxes.downwardFluxWPerM2;
            point.surfaceAirFluxWPerM2 = 0.0;
            point.subsurfaceFluxWPerM2 = 0.0;
            if (logDetails) {
                qCInfo(solarRadiationLog) << "Resolved air temperature (init)"
                                          << "index=" << i
                                          << "airTemperatureK=" << point.airTemperatureK;
            }
            minAirTemperature = qMin(minAirTemperature, point.airTemperatureK);
            maxAirTemperature = qMax(maxAirTemperature, point.airTemperatureK);
        }

        result.hasAirTemperatureRange = minAirTemperature <= maxAirTemperature;
        result.minAirTemperatureK = minAirTemperature;
        result.maxAirTemperatureK = maxAirTemperature;
        result.hasPressureRange = minPressureAtm <= maxPressureAtm;
        result.minPressureAtm = minPressureAtm;
        result.maxPressureAtm = maxPressureAtm;

        double meanSurfaceTemperatureKelvin = 0.0;
        int meanSurfaceSamples = 0;
        for (const auto &point : result.grid.points()) {
            if (point.temperatureK > 0.0) {
                meanSurfaceTemperatureKelvin += point.temperatureK;
                ++meanSurfaceSamples;
            }
        }
        if (meanSurfaceSamples > 0) {
            meanSurfaceTemperatureKelvin /= static_cast<double>(meanSurfaceSamples);
        }

        double baselineMeanTemperature = 0.0;
        int baselineSamples = 0;
        for (const auto &temperature : baselineAirTemperatures) {
            if (temperature > 0.0) {
                baselineMeanTemperature += temperature;
                ++baselineSamples;
            }
        }
        if (baselineSamples > 0) {
            baselineMeanTemperature /= static_cast<double>(baselineSamples);
        }

        // Базовый профиль атмосферы синхронизируем с рассчитанной поверхностью,
        // чтобы слои получали ненулевую температуру и давление на старте.
        double baseTemperatureKelvin =
            (effectiveTemperatureKelvin > 0.0) ? effectiveTemperatureKelvin
                                                : meanSurfaceTemperatureKelvin;
        if (baseTemperatureKelvin <= 0.0) {
            // Если и t_eff, и средняя поверхность равны нулю, используем разумный fallback,
            // чтобы профиль атмосферы не был полностью нулевым.
            baseTemperatureKelvin =
                (baselineMeanTemperature > 0.0) ? baselineMeanTemperature
                                                : qMax(1.0, input.stateDefaults.minTemperatureKelvin);
        }
        if (meanSurfaceTemperatureKelvin > 0.0) {
            result.grid.initializeAtmosphericGrid(input.atmosphere,
                                                  input.massEarths,
                                                  baseTemperatureKelvin,
                                                  0,
                                                  input.minBottomLayerThicknessMeters,
                                                  input.minTopPressureAtm);
        }

        if (shouldCancel()) {
            result.cancelled = true;
            return result;
        }

        constexpr double kMinAtmospherePressureAtm = 1e-6;
        if (result.grid.points().isEmpty() ||
            atmospherePressureAtm <= kMinAtmospherePressureAtm ||
            input.atmosphere.totalMassGigatons() <= 0.0) {
            for (auto &point : result.grid.points()) {
                point.windEastMps = 0.0;
                point.windNorthMps = 0.0;
                point.windSpeedMps = 0.0;
            }
            result.hasWindRange = false;
            return result;
        }

        QVector<double> pressuresAtm;
        QVector<double> temperatures;
        pressuresAtm.reserve(result.grid.points().size());
        temperatures.reserve(result.grid.points().size());
        for (const auto &point : result.grid.points()) {
            pressuresAtm.push_back(qMax(0.0, point.pressureAtm));
            temperatures.push_back(point.temperatureK);
        }

        WindFieldModel windModel;
        const QVector<WindVector> wind =
            windModel.buildField(result.grid,
                                 pressuresAtm,
                                 temperatures,
                                 qMax(0.0, dayLengthDays) * 86400.0,
                                 2);
        if (wind.size() != result.grid.points().size()) {
            for (auto &point : result.grid.points()) {
                point.windEastMps = 0.0;
                point.windNorthMps = 0.0;
                point.windSpeedMps = 0.0;
            }
            result.hasWindRange = false;
            return result;
        }

        double minWind = std::numeric_limits<double>::max();
        double maxWind = std::numeric_limits<double>::lowest();
        for (int i = 0; i < wind.size(); ++i) {
            if ((i % 64) == 0 && shouldCancel()) {
                result.cancelled = true;
                return result;
            }
            auto &point = result.grid.points()[i];
            point.windEastMps = wind.at(i).eastMps;
            point.windNorthMps = wind.at(i).northMps;
            point.windSpeedMps = std::hypot(point.windEastMps, point.windNorthMps);
            minWind = qMin(minWind, point.windSpeedMps);
            maxWind = qMax(maxWind, point.windSpeedMps);
        }

        if (minWind <= maxWind) {
            result.hasWindRange = true;
            result.minWindSpeedMps = minWind;
            result.maxWindSpeedMps = maxWind;
        } else {
            result.hasWindRange = false;
        }

        return result;
    }

    void handleSurfaceGridComputationFinished(int requestId) {
        if (!surfaceGridWatcher_) {
            return;
        }
        if (!surfaceGridRequestId_ || requestId != surfaceGridRequestId_->load()) {
            return;
        }

        SurfaceGridComputationResult result = surfaceGridWatcher_->future().result();
        if (result.cancelled) {
            setSurfaceGridCalculationRunning(false);
            return;
        }

        surfaceGrid_ = std::move(result.grid);
        updateSubsurfaceSettingsUi(result.resolvedSubsurfaceSettings);
        applySurfaceGridToViews();
        updateSurfaceHeightLegendFromGrid();
        const double distanceAu = surfaceOrbitAnimation_.distanceAU();
        const double siderealDayLengthDays =
            resolveSiderealDayLengthDays(result.dayLengthDays, result.yearLengthDays,
                                         result.isRetrograde);
        const double elapsedTimeDays = result.elapsedDays * result.dayLengthDays;
        updateSurfaceStarRendering(result.declinationDegrees,
                                   elapsedTimeDays,
                                   siderealDayLengthDays,
                                   result.rotationMode,
                                   result.orbitalPhaseRadians,
                                   result.spinOrbitP,
                                   result.spinOrbitQ,
                                   result.isRetrograde,
                                   distanceAu);
        if (result.hasTemperatureRange && result.minTemperatureK <= result.maxTemperatureK) {
            applySurfaceTemperatureRangeToViews(result.minTemperatureK, result.maxTemperatureK);
            updateSurfaceTemperatureLegend(true, result.minTemperatureK, result.maxTemperatureK);
        } else {
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
        }
        if (result.hasAirTemperatureRange &&
            result.minAirTemperatureK <= result.maxAirTemperatureK) {
            if (surfaceMapMode_ == SurfaceMapMode::AirTemperature) {
                applySurfaceTemperatureRangeToViews(result.minAirTemperatureK,
                                                    result.maxAirTemperatureK);
            }
            updateSurfaceAirTemperatureLegend(true,
                                              result.minAirTemperatureK,
                                              result.maxAirTemperatureK);
        } else {
            updateSurfaceAirTemperatureLegend(false, 0.0, 0.0);
        }
        if (result.hasPressureRange && result.minPressureAtm <= result.maxPressureAtm) {
            applySurfacePressureRangeToViews(result.minPressureAtm, result.maxPressureAtm);
            updateSurfacePressureLegend(true, result.minPressureAtm, result.maxPressureAtm);
        } else {
            updateSurfacePressureLegend(false, 0.0, 0.0);
        }
        if (result.hasWindRange && result.minWindSpeedMps <= result.maxWindSpeedMps) {
            applySurfaceWindRangeToViews(result.minWindSpeedMps, result.maxWindSpeedMps);
            updateSurfaceWindLegend(true, result.minWindSpeedMps, result.maxWindSpeedMps);
        } else {
            updateSurfaceWindLegend(false, 0.0, 0.0);
        }

        setSurfaceGridCalculationRunning(false);
        resumeSurfaceSimulationAfterGridUpdate();
    }

    void startSurfaceGridTemperatureUpdate() {
        cancelSurfaceGridTemperatureUpdate();
        const bool wasRunning = surfaceSimRunning_;
        resetSurfaceSimulation(false);
        surfaceSimResumeAfterGridUpdate_ = wasRunning;
        rebuildSurfaceGrid();
        surfaceGridDirty_ = false;
        updateSurfaceHeightLegendFromGrid();

        if (surfaceGrid_.points().isEmpty()) {
            applySurfaceGridToViews();
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceAirTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceWindLegend(false, 0.0, 0.0);
            updateSurfacePressureLegend(false, 0.0, 0.0);
            setSurfaceGridCalculationRunning(false);
            resumeSurfaceSimulationAfterGridUpdate();
            return;
        }

        auto stateDefaults = buildSurfacePointStateDefaults();
        if (!stateDefaults) {
            applySurfaceGridToViews();
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceAirTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceWindLegend(false, 0.0, 0.0);
            updateSurfacePressureLegend(false, 0.0, 0.0);
            setSurfaceGridCalculationRunning(false);
            resumeSurfaceSimulationAfterGridUpdate();
            return;
        }

        const AtmosphereComposition atmosphere = currentAtmosphereForCalculations();
        const double massEarths = planetComboBox_->currentData(kRoleMassEarths).toDouble();
        const double radiusKm = planetComboBox_->currentData(kRoleRadiusKm).toDouble();
        const double cloudAlbedo =
            qBound(0.0, planetComboBox_->currentData(kRoleCloudAlbedo).toDouble(), 1.0);
        const bool disableWaterAndClouds =
            planetComboBox_->currentData(kRoleDisableWaterAndClouds).toBool();
        const QString baseMaterialId =
            resolveSurfaceMaterialIdForPoint(atmosphere, massEarths, radiusKm);

        QHash<QString, SurfaceMaterial> materialsById;
        const auto materials = surfaceMaterials();
        materialsById.reserve(materials.size());
        for (const auto &material : materials) {
            materialsById.insert(material.id, material);
        }

        const double dayLengthDays = planetComboBox_->currentData(kRoleDayLength).toDouble();
        const double yearLengthDays =
            planetComboBox_->currentData(kRoleYearLengthDays).toDouble();
        const bool isRetrograde =
            isRetrogradeRotation(planetComboBox_->currentData(kRoleObliquity).toDouble());
        const double siderealDayLengthDays =
            resolveSiderealDayLengthDays(dayLengthDays, yearLengthDays, isRetrograde);
        const int stepsPerDay = qMax(1, qRound(siderealDayLengthDays * 24.0));
        const int spinOrbitP = planetComboBox_->currentData(kRoleSpinOrbitP).toInt();
        const int spinOrbitQ = planetComboBox_->currentData(kRoleSpinOrbitQ).toInt();
        const bool manualGreenhouseOnTop =
            planetComboBox_->currentData(kRoleManualGreenhouseOnTopOfAtmosphere).toBool();
        const auto radiationModelType = static_cast<RadiationModelType>(
            planetComboBox_->currentData(kRoleAdvancedRadiationModel).toInt());
        stateDefaults->radiationModelType = radiationModelType;
        const RotationMode rotationMode =
            resolveRotationModeForPlanetIndex(planetComboBox_->currentIndex());
        ensureSurfaceOrbitAnimationReady();
        const double declinationDegrees = surfaceOrbitAnimation_.declinationDegrees();
        const double orbitalPhaseRadians = surfaceOrbitAnimation_.orbitalPhaseRadians();
        double segmentSolarConstant = hasSolarConstant_ ? lastSolarConstant_ : 0.0;
        if (hasSolarConstant_ && lastSolarConstantDistanceAU_ > 0.0) {
            const double distanceAU = surfaceOrbitAnimation_.distanceAU();
            segmentSolarConstant =
                lastSolarConstant_ *
                std::pow(lastSolarConstantDistanceAU_ / distanceAU, 2.0);
            updateSurfaceSolarConstantLabel(computeSurfaceFluxBreakdown());
        } else {
            const SurfaceFluxBreakdown breakdown = computeSurfaceFluxBreakdown();
            if (!hasSolarConstant_ && breakdown.totalFluxWPerM2 > 0.0) {
                // Допускаем оценку по пресету для первичной визуализации тайлов,
                // чтобы не держать карту на фоне 3 K до ручного расчёта.
                segmentSolarConstant = breakdown.totalFluxWPerM2;
            }
            updateSurfaceSolarConstantLabel(breakdown.totalFluxWPerM2 > 0.0
                                                ? breakdown
                                                : SurfaceFluxBreakdown{});
        }
        const bool skipSpinUp = surfaceSkipSpinUp_;
        surfaceSkipSpinUp_ = false;

        double atmospherePressureAtm = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            atmospherePressureAtm = atmosphere.totalPressureAtm(massEarths, radiusKm);
        }
        const double atmosphereMassGt = atmosphere.totalMassGigatons();
        const double co2MassGt = atmosphere.massGigatons(QStringLiteral("co2"));
        const double co2Share = atmosphereMassGt > 0.0 ? co2MassGt / atmosphereMassGt : 0.0;
        const bool useAtmosphericModel = atmosphereMassGt > 0.0;
        const bool hasAtmospherePressure = atmospherePressureAtm > 0.0;
        const bool atmosphereEnabled = useAtmosphericModel || hasAtmospherePressure;

        double surfaceGravity = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            const double radiusMeters = radiusKm * 1000.0;
            const double planetMassKg = massEarths * kEarthMassKg;
            surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);
        }
        const double gravity = (surfaceGravity > 0.0) ? surfaceGravity : 9.80665;

        const double meanToaFlux =
            segmentSolarConstant * qMax(0.0, 1.0 - cloudAlbedo) / 4.0;
        const double effectiveTemperatureKelvin =
            std::pow(qMax(0.0, meanToaFlux) / kStefanBoltzmannConstant, 0.25);
        const double minDenseAtmosphereTemperatureK =
            currentMinDenseAtmosphereTemperatureKelvin();
        if (atmosphereEnabled && meanToaFlux > 0.0) {
            // Учитываем атмосферный парник уже в инициализации:
            // без этого плотные атмосферы могут «схлопнуться» в холодный режим на старте.
            stateDefaults->greenhouseOpacity =
                computeLocalGreenhouseOpacity(atmosphere,
                                              stateDefaults->material,
                                              cloudAlbedo,
                                              disableWaterAndClouds,
                                              qMax(0.0, atmospherePressureAtm),
                                              radiusKm,
                                              gravity,
                                              meanToaFlux,
                                              meanToaFlux,
                                              effectiveTemperatureKelvin,
                                              stateDefaults->manualGreenhouseOpacity,
                                              minDenseAtmosphereTemperatureK,
                                              useAtmosphericModel,
                                              radiationModelType,
                                              manualGreenhouseOnTop,
                                              false);
        } else {
            stateDefaults->greenhouseOpacity = stateDefaults->manualGreenhouseOpacity;
        }

        double meanSurfaceTemperatureKelvin = 0.0;
        int meanSurfaceSamples = 0;
        for (const auto &point : surfaceGrid_.points()) {
            if (point.temperatureK > 0.0) {
                meanSurfaceTemperatureKelvin += point.temperatureK;
                ++meanSurfaceSamples;
            }
        }
        if (meanSurfaceSamples > 0) {
            meanSurfaceTemperatureKelvin /= static_cast<double>(meanSurfaceSamples);
        }

        // Базу берём из t_eff (или текущей средней поверхности), чтобы старт атмосферы
        // был синхронизирован с картой поверхности и не давал скачка при первом тике.
        double baseTemperatureKelvin =
            (effectiveTemperatureKelvin > 0.0) ? effectiveTemperatureKelvin
                                               : meanSurfaceTemperatureKelvin;
        if (baseTemperatureKelvin <= 0.0) {
            // Fallback нужен, чтобы при нулевых входных потоках профиль атмосферы
            // не оставался в 0 K: это обеспечивает ненулевые T/P в слоях (например, Земля).
            baseTemperatureKelvin = qMax(1.0, stateDefaults->minTemperatureKelvin);
        }
        if (atmospherePressureAtm >= 0.1) {
            // Для плотных атмосфер профиль должен стартовать не от t_eff (ТОА),
            // а от нижней границы парникового прогрева, задаваемой давлением и составом.
            const double denseAtmosphereBase =
                denseAtmosphereMinTemperatureKelvin(effectiveTemperatureKelvin,
                                                    atmospherePressureAtm,
                                                    stateDefaults->greenhouseOpacity,
                                                    co2Share,
                                                    minDenseAtmosphereTemperatureK);
            baseTemperatureKelvin = qMax(baseTemperatureKelvin, denseAtmosphereBase);
        }
        surfaceGrid_.initializeAtmosphericGrid(atmosphere,
                                               massEarths,
                                               baseTemperatureKelvin,
                                               0,
                                               currentAtmosphereBottomLayerThicknessMeters(),
                                               currentMinTopPressureAtm());

        SurfaceGridComputationInput input;
        input.grid = surfaceGrid_;
        input.atmosphere = atmosphere;
        input.stateDefaults = *stateDefaults;
        input.materialsById = materialsById;
        input.baseMaterialId = baseMaterialId;
        input.massEarths = massEarths;
        input.radiusKm = radiusKm;
        input.cloudAlbedo = cloudAlbedo;
        input.disableWaterAndClouds = disableWaterAndClouds;
        input.dayLengthDays = dayLengthDays;
        input.yearLengthDays = yearLengthDays;
        input.elapsedDays = surfaceTime_.surfaceElapsedDays();
        input.stepsPerDay = stepsPerDay;
        input.spinOrbitP = spinOrbitP;
        input.spinOrbitQ = spinOrbitQ;
        input.isRetrograde = isRetrograde;
        input.manualGreenhouseOpacity = stateDefaults->manualGreenhouseOpacity;
        input.manualGreenhouseOnTop = manualGreenhouseOnTop;
        input.radiationModelType = radiationModelType;
        input.minDenseAtmosphereTemperatureK = currentMinDenseAtmosphereTemperatureKelvin();
        input.rotationMode = rotationMode;
        input.declinationDegrees = declinationDegrees;
        input.orbitalPhaseRadians = orbitalPhaseRadians;
        input.segmentSolarConstant = segmentSolarConstant;
        input.hasSolarConstant = segmentSolarConstant > 0.0;
        input.currentHourIndex = surfaceTime_.hourIndex;
        input.currentAbsoluteHourIndex = surfaceTime_.surfaceTotalHourIndex();
        input.verticalWindMixingCoefficientKz = currentVerticalWindMixingCoefficientKz();
        input.logPointIndex = selectedSurfacePointIndex_;
        input.skipSpinUp = skipSpinUp;
        input.minBottomLayerThicknessMeters = currentAtmosphereBottomLayerThicknessMeters();
        input.minTopPressureAtm = currentMinTopPressureAtm();

        surfaceGridCancelFlag_ = std::make_shared<std::atomic_bool>(false);
        if (!surfaceGridRequestId_) {
            surfaceGridRequestId_ = std::make_shared<std::atomic<int>>(0);
        }
        // Храним счётчик запросов в shared_ptr, чтобы фоновые задачи не
        // держали указатель на уже уничтоженный виджет.
        const int requestId = surfaceGridRequestId_->fetch_add(1) + 1;

        if (!surfaceGridWatcher_) {
            surfaceGridWatcher_ = new QFutureWatcher<SurfaceGridComputationResult>(this);
        }
        surfaceGridWatcher_->disconnect(this);
        connect(surfaceGridWatcher_, &QFutureWatcher<SurfaceGridComputationResult>::finished, this,
                [this, requestId]() { handleSurfaceGridComputationFinished(requestId); });

        setSurfaceGridCalculationRunning(true);

        // Запуск в фоне позволяет отменять устаревшие расчёты при быстрой смене параметров.
        QFuture<SurfaceGridComputationResult> future = QtConcurrent::run(
            [input,
             cancelFlag = surfaceGridCancelFlag_,
             requestIdCounter = surfaceGridRequestId_,
             requestId]() {
                return SolarCalculatorWidget::computeSurfaceGridTemperatures(input,
                                                                             cancelFlag,
                                                                             requestIdCounter,
                                                                             requestId);
            });
        surfaceGridWatcher_->setFuture(future);
    }

    void updateSurfaceGridTemperatures() {
        scheduleSurfaceGridTemperatureUpdate(false);
    }

    void toggleSurfaceSimulation() {
        if (surfaceSimRunning_) {
            surfaceSimRunning_ = false;
            if (surfaceSimTimer_) {
                surfaceSimTimer_->stop();
            }
            updateSurfaceSimulationUi();
            updateStarSystemOrbits();
            return;
        }

        if (!hasSolarConstant_ || planetComboBox_->currentIndex() < 0) {
            return;
        }

        if (surfaceGrid_.points().isEmpty()) {
            updateSurfaceGridTemperatures();
        }

        if (!surfaceSimTimer_) {
            surfaceSimTimer_ = new QTimer(this);
            connect(surfaceSimTimer_, &QTimer::timeout, this,
                    [this]() { advanceSurfaceSimulation(); });
        }
        updateSurfaceSimulationTimerInterval();

        surfaceSimRunning_ = true;
        updateSurfaceSimulationUi();
        if (!temperaturePauseFlag_.load()) {
            surfaceSimTimer_->start();
        }
        updateStarSystemOrbits();
    }

    void advanceSurfaceSimulation() {
        if (!surfaceSimRunning_ || surfaceGrid_.points().isEmpty()) {
            return;
        }

        const auto defaultMaterial = currentMaterial();
        if (!defaultMaterial) {
            return;
        }

        const double dayLengthDays = planetComboBox_->currentData(kRoleDayLength).toDouble();
        const RotationMode rotationMode =
            resolveRotationModeForPlanetIndex(planetComboBox_->currentIndex());
        const double yearLengthDays =
            planetComboBox_->currentData(kRoleYearLengthDays).toDouble();
        const bool isRetrograde =
            isRetrogradeRotation(planetComboBox_->currentData(kRoleObliquity).toDouble());
        const int spinOrbitP = planetComboBox_->currentData(kRoleSpinOrbitP).toInt();
        const int spinOrbitQ = planetComboBox_->currentData(kRoleSpinOrbitQ).toInt();
        const AtmosphereComposition atmosphere = currentAtmosphereForCalculations();
        const double massEarths = planetComboBox_->currentData(kRoleMassEarths).toDouble();
        const double radiusKm = planetComboBox_->currentData(kRoleRadiusKm).toDouble();
        double atmospherePressureAtm = 0.0;
        double surfaceGravity = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            atmospherePressureAtm = atmosphere.totalPressureAtm(massEarths, radiusKm);
            const double radiusMeters = radiusKm * 1000.0;
            const double planetMassKg = massEarths * kEarthMassKg;
            surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);
        }
        const double gravity = (surfaceGravity > 0.0) ? surfaceGravity : 9.80665;
        const double atmosphereMassGt = atmosphere.totalMassGigatons();
        const double co2MassGt = atmosphere.massGigatons(QStringLiteral("co2"));
        const double co2Share = atmosphereMassGt > 0.0 ? co2MassGt / atmosphereMassGt : 0.0;
        const double localSeaLevelPressureAtm =
            calculateCellPressureAtmFromKg(atmosphere.totalMassKg(),
                                           massEarths,
                                           radiusKm,
                                           surfaceGrid_.pointAreaKm2());
        const bool useAtmosphericModel = atmosphere.totalMassGigatons() > 0.0;
        const double manualGreenhouseOpacity =
            planetComboBox_->currentData(kRoleGreenhouseOpacity).toDouble();
        const bool manualGreenhouseOnTop =
            planetComboBox_->currentData(kRoleManualGreenhouseOnTopOfAtmosphere).toBool();
        const auto radiationModelType = static_cast<RadiationModelType>(
            planetComboBox_->currentData(kRoleAdvancedRadiationModel).toInt());
        const double cloudAlbedo =
            qBound(0.0, planetComboBox_->currentData(kRoleCloudAlbedo).toDouble(), 1.0);
        const bool disableWaterAndClouds =
            planetComboBox_->currentData(kRoleDisableWaterAndClouds).toBool();

        QHash<QString, SurfaceMaterial> materialsById;
        const auto materials = surfaceMaterials();
        materialsById.reserve(materials.size());
        for (const auto &material : materials) {
            materialsById.insert(material.id, material);
        }
        const auto materialForPoint = [&materialsById, &defaultMaterial](const QString &materialId) {
            const auto it = materialsById.constFind(materialId);
            return it != materialsById.cend() ? *it : *defaultMaterial;
        };

        syncSurfaceOrbitAnimationToTime();
        const double siderealDayLengthDays =
            resolveSiderealDayLengthDays(dayLengthDays, yearLengthDays, isRetrograde);
        const int siderealStepsPerDay = qMax(1, qRound(siderealDayLengthDays * 24.0));
        // Сезонная деклинация: δ = asin(sin(наклон оси) * sin(истинная долгота звезды)).
        const double declinationDegrees = surfaceOrbitAnimation_.declinationDegrees();
        const double orbitalPhaseRadians = surfaceOrbitAnimation_.orbitalPhaseRadians();

        // Инсоляция меняется с расстоянием как 1 / r^2 относительно опорной дистанции.
        const double distanceAU = surfaceOrbitAnimation_.distanceAU();
        const double segmentSolarConstant =
            lastSolarConstant_ *
            std::pow(lastSolarConstantDistanceAU_ / distanceAU, 2.0);
        updateSurfaceSolarConstantLabel(computeSurfaceFluxBreakdown());
        // Один тик = 1 час планетарных суток, ускорение реализовано уменьшением интервала таймера.
        const double timeStepSeconds = 3600.0;
        // Коэффициент прохождения коротковолнового излучения через облака.
        double cloudShortwaveTransmission = 1.0 - cloudAlbedo;
        if (cloudAlbedo > 0.7) {
            // Для плотных сернокислотных облаков дополнительно ослабляем поток к поверхности.
            cloudShortwaveTransmission *= 0.2;
        }
        cloudShortwaveTransmission = qBound(0.0, cloudShortwaveTransmission, 1.0);
        auto &atmosphereGrid = surfaceGrid_.atmosphericGrid();
        const bool useLayeredAtmosphere =
            useAtmosphericModel && radiationModelType == RadiationModelType::Layered &&
            atmosphereGrid.columnCount() == surfaceGrid_.points().size() &&
            atmosphereGrid.layerCount() > 0;
        const int solarStepsPerDay = qMax(1, surfaceTime_.stepsPerDay);
        const double elapsedTimeDays =
            (surfaceTime_.surfaceElapsedDays() + 0.5 / static_cast<double>(solarStepsPerDay)) *
            dayLengthDays;
        const double substellarLongitudeRadians =
            resolveSubstellarLongitudeRadians(elapsedTimeDays,
                                              siderealDayLengthDays,
                                              rotationMode,
                                              orbitalPhaseRadians,
                                              spinOrbitP,
                                              spinOrbitQ,
                                              isRetrograde);
        const double declinationRadians = qDegreesToRadians(declinationDegrees);
        const double sinDeclination = std::sin(declinationRadians);
        const double cosDeclination = std::cos(declinationRadians);
        // Радиационный профиль должен опираться на средний ТОА-поток, а не на локальную
        // инсоляцию: вертикальная структура атмосферы задаётся глобальным балансом энергии,
        // а локальная инсоляция влияет только на SW-поток у поверхности и в воздухе.
        const double meanToaFlux =
            segmentSolarConstant * qMax(0.0, 1.0 - cloudAlbedo) / 4.0;
        const double effectiveTemperatureKelvin =
            std::pow(qMax(0.0, meanToaFlux) / kStefanBoltzmannConstant, 0.25);
        // Фиксируем меридиан для агрегации: либо пользовательскую точку,
        // либо ближайший к подсолнечному меридиану, чтобы брать единый срез по широтам.
        const double targetLongitudeRadians =
            surfaceTemperatureAggregation_.hasTargetLongitude
                ? surfaceTemperatureAggregation_.targetLongitudeRadians
                : selectSurfaceAggregationLongitude(substellarLongitudeRadians);
        ensureSurfaceTemperatureAggregationTargets(targetLongitudeRadians);

        QVector<double> localInsolations;
        localInsolations.reserve(surfaceGrid_.points().size());
        for (auto &point : surfaceGrid_.points()) {
            const double localHourAngle = point.longitudeRadians - substellarLongitudeRadians;
            const double cosZenith =
                point.sinLatitude * sinDeclination +
                point.cosLatitude * cosDeclination * std::cos(localHourAngle);
            // S_inst = S0 * cos(zenith) при освещении, иначе 0.
            const double localInsolation =
                segmentSolarConstant * qMax(0.0, cosZenith);
            localInsolations.push_back(localInsolation);
            // Стабилизация для плотных атмосфер (см. computeLocalGreenhouseOpacity):
            // t_eff недостаточна для сверхплотных атмосфер (например, Венера),
            // поэтому минимальная база зависит от давления и состава.
            const double profileBaseTemperatureKelvin =
                (point.pressureAtm >= 0.1)
                    ? qMax(point.state.temperatureKelvin(),
                           denseAtmosphereMinTemperatureKelvin(effectiveTemperatureKelvin,
                                                               point.pressureAtm,
                                                               point.state.greenhouseOpacity(),
                                                               co2Share,
                                                               currentMinDenseAtmosphereTemperatureKelvin()))
                        : point.state.temperatureKelvin();
            const auto radiationModel =
                makeRadiationModel(atmosphere,
                                   point.pressureAtm,
                                   profileBaseTemperatureKelvin,
                                   effectiveTemperatureKelvin,
                                   gravity,
                                   radiationModelType);
            // Поток у поверхности: S_surf = S_local * T_cloud * T_atm.
            const double surfaceShortwaveFlux =
                localInsolation * cloudShortwaveTransmission *
                radiationModel->incomingTransmission();
            // Запоминаем солнечный поток у поверхности для UI и логов.
            point.solarFluxWPerM2 = surfaceShortwaveFlux;
            point.shortwaveSurfaceWPerM2 = surfaceShortwaveFlux;

            const double emittedFlux = point.state.surfaceEmittedFlux();
            const double surfaceLongwaveUp =
                resolveSurfaceLongwaveUpFlux(*radiationModel, radiationModelType, emittedFlux);
            const double airTemperatureForFlux =
                (point.airTemperatureK > 0.0) ? point.airTemperatureK : point.temperatureK;
            const double greenhouseOpacity = point.state.greenhouseOpacity();
            const LongwaveFluxes longwaveFluxes =
                computeLongwaveFluxes(*radiationModel,
                                      radiationModelType,
                                      emittedFlux,
                                      airTemperatureForFlux,
                                      greenhouseOpacity,
                                      useAtmosphericModel,
                                      manualGreenhouseOnTop);
            point.longwaveUpWPerM2 = surfaceLongwaveUp;
            point.longwaveDownWPerM2 = longwaveFluxes.downwardFluxWPerM2;
            point.surfaceAirFluxWPerM2 = 0.0;
            // Баланс поверхности: учитываем нисходящее LW излучение атмосферы.
            const double absorbedFlux =
                point.state.absorbedFlux(surfaceShortwaveFlux) + longwaveFluxes.downwardFluxWPerM2;
            // Положительный поток означает перенос энергии в грунт (вниз).
            point.subsurfaceFluxWPerM2 =
                point.state.stabilizedRadiativeFlux(absorbedFlux, timeStepSeconds);
            // Применяем шаговый радиационный баланс для состояния точки.
            point.state.updateTemperature(absorbedFlux, emittedFlux, timeStepSeconds);
            // Обновляем поверхностную температуру до переноса в соседние точки.
            point.temperatureK = point.state.temperatureKelvin();
        }

        updateSurfaceWindField(atmosphere, atmospherePressureAtm, dayLengthDays, surfaceGravity);
        // Перенос тепла теперь полностью зависит от расчёта ветра/адвекции без радиационного смешивания.

        QVector<double> temperatures;
        QVector<double> pressuresAtm;
        QVector<double> windEast;
        QVector<double> windNorth;
        temperatures.reserve(surfaceGrid_.points().size());
        pressuresAtm.reserve(surfaceGrid_.points().size());
        windEast.reserve(surfaceGrid_.points().size());
        windNorth.reserve(surfaceGrid_.points().size());
        for (const auto &point : surfaceGrid_.points()) {
            temperatures.push_back(point.temperatureK);
            pressuresAtm.push_back(qMax(0.0, point.pressureAtm));
            windEast.push_back(point.windEastMps);
            windNorth.push_back(point.windNorthMps);
        }

        // Переносим поверхностное давление (уровень рельефа) по полю ветра.
        SurfacePressureTransportModel pressureTransportModel;
        QVector<double> advectedPressures =
            pressureTransportModel.advectPressure(surfaceGrid_,
                                                  pressuresAtm,
                                                  windEast,
                                                  windNorth,
                                                  timeStepSeconds,
                                                  1,
                                                  0.0);
        if (advectedPressures.size() != surfaceGrid_.points().size()) {
            advectedPressures = pressuresAtm;
        }
        constexpr double kPressureRelaxFactor = 0.1;
        double minPressureAtm = std::numeric_limits<double>::max();
        double maxPressureAtm = std::numeric_limits<double>::lowest();
        for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
            auto &point = surfaceGrid_.points()[i];
            const double advectedAtm = advectedPressures.at(i);
            double airColumnTemperatureK = point.airTemperatureK;
            if (airColumnTemperatureK <= 0.0 && radiationModelType == RadiationModelType::Layered) {
                // Стабилизация для плотных атмосфер (см. computeLocalGreenhouseOpacity):
                // t_eff недостаточна для сверхплотных атмосфер (например, Венера),
                // поэтому минимальная база зависит от давления и состава.
                const double profileBaseTemperatureKelvin =
                    (point.pressureAtm >= 0.1)
                        ? qMax(point.state.temperatureKelvin(),
                               denseAtmosphereMinTemperatureKelvin(effectiveTemperatureKelvin,
                                                                   point.pressureAtm,
                                                                   point.state.greenhouseOpacity(),
                                                                   co2Share,
                                                                   currentMinDenseAtmosphereTemperatureKelvin()))
                        : point.state.temperatureKelvin();
                const auto radiationModel =
                    makeRadiationModel(atmosphere,
                                       point.pressureAtm,
                                       profileBaseTemperatureKelvin,
                                       effectiveTemperatureKelvin,
                                       gravity,
                                       radiationModelType);
                const auto *layeredModel =
                    dynamic_cast<const LayeredRadiationModel *>(radiationModel.get());
                if (layeredModel) {
                    airColumnTemperatureK = layeredModel->bottomLayerTemperatureKelvin();
                }
            }
            if (airColumnTemperatureK <= 0.0) {
                airColumnTemperatureK = effectiveTemperatureKelvin;
            }
            // Барометрическая релаксация возвращает поле давления к
            // физически согласованному профилю P(z) при температуре столба воздуха
            // (а не поверхности), чтобы ночное охлаждение не схлопывало давление.
            const double barometricAtm =
                AtmosphericPressureModel::pressureAtHeightAtm(localSeaLevelPressureAtm,
                                                              point.heightKm * 1000.0,
                                                              airColumnTemperatureK,
                                                              atmosphere,
                                                              gravity);
            const double relaxedAtm =
                advectedAtm + (barometricAtm - advectedAtm) * kPressureRelaxFactor;
            // Фиксируем перенесённое поверхностное давление для UI и динамики.
            point.pressureAtm = qMax(0.0, relaxedAtm);
            minPressureAtm = qMin(minPressureAtm, point.pressureAtm);
            maxPressureAtm = qMax(maxPressureAtm, point.pressureAtm);
            const double localInsolation =
                (i < localInsolations.size()) ? localInsolations.at(i) : 0.0;
            const SurfaceMaterial material = materialForPoint(point.materialId);
            const bool logDetails = shouldLogRadiationForPoint(i);
            const double greenhouseBaseTemperature =
                (radiationModelType == RadiationModelType::Layered && point.airTemperatureK > 0.0)
                    ? point.airTemperatureK
                    : point.state.temperatureKelvin();
            const double localGreenhouseOpacity =
                computeLocalGreenhouseOpacity(atmosphere,
                                              material,
                                              cloudAlbedo,
                                              disableWaterAndClouds,
                                              point.pressureAtm,
                                              radiusKm,
                                              gravity,
                                              localInsolation,
                                              meanToaFlux,
                                              greenhouseBaseTemperature,
                                              manualGreenhouseOpacity,
                                              currentMinDenseAtmosphereTemperatureKelvin(),
                                              useAtmosphericModel,
                                              radiationModelType,
                                              manualGreenhouseOnTop,
                                              logDetails);
            point.state.setGreenhouseOpacity(localGreenhouseOpacity);
        }

        SurfaceAdvectionModel advectionModel;
        constexpr double kMinAtmospherePressureAtm = 1e-6;
        QVector<double> advectedTemperatures;
        if (atmospherePressureAtm <= kMinAtmospherePressureAtm ||
            atmosphere.totalMassGigatons() <= 0.0) {
            // На безатмосферных телах нет ветрового переноса тепла.
            advectedTemperatures = temperatures;
        } else {
            advectedTemperatures =
                advectionModel.advectTemperature(surfaceGrid_,
                                                 temperatures,
                                                 windEast,
                                                 windNorth,
                                                 timeStepSeconds,
                                                 1,
                                                 1.0);
        }
        if (advectedTemperatures.size() != surfaceGrid_.points().size()) {
            advectedTemperatures = temperatures;
        }

        for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
            auto &point = surfaceGrid_.points()[i];
            // Перенос — поверхностная величина и не должен переписывать профиль глубины.
            point.state.setSurfaceLayerTemperatureKelvin(advectedTemperatures.at(i));
            // Температура после переноса хранится как поверхностная величина.
            point.temperatureK = point.state.temperatureKelvin();
        }

        if (useLayeredAtmosphere) {
            const double dayLengthSeconds = qMax(0.0, dayLengthDays) * 86400.0;
            AtmosphereStepSolver stepSolver(atmosphere,
                                            gravity,
                                            timeStepSeconds,
                                            dayLengthSeconds,
                                            isRetrograde);
            AtmosphereStepSolver::LayeredStepInput stepInput{surfaceGrid_,
                                                            atmosphereGrid,
                                                            localInsolations,
                                                            materialsById,
                                                            *defaultMaterial,
                                                            cloudShortwaveTransmission,
                                                            kDefaultHeatTransferWPerM2K};
            stepInput.verticalWindMixingCoefficientKz = currentVerticalWindMixingCoefficientKz();
            stepInput.logPointIndex = selectedSurfacePointIndex_;
            stepSolver.runLayeredStep(stepInput);
        }

        double minTemperature = std::numeric_limits<double>::max();
        double maxTemperature = std::numeric_limits<double>::lowest();
        double minAirTemperature = std::numeric_limits<double>::max();
        double maxAirTemperature = std::numeric_limits<double>::lowest();
        for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
            auto &point = surfaceGrid_.points()[i];
            const double initialAirTemperature =
                (point.airTemperatureK > 0.0) ? point.airTemperatureK : point.temperatureK;
            const double localInsolation =
                (i < localInsolations.size()) ? localInsolations.at(i) : 0.0;
            const SurfaceMaterial material = materialForPoint(point.materialId);
            const bool logDetails = shouldLogRadiationForPoint(i);
            if (logDetails) {
                qCInfo(solarRadiationLog) << "Radiation inputs (tick)"
                                          << "index=" << i
                                          << "localInsolation=" << localInsolation
                                          << "planetaryAlbedo=" << cloudAlbedo
                                          << "meanToaFlux=" << meanToaFlux
                                          << "effectiveTemperatureKelvin=" << effectiveTemperatureKelvin;
            }
            if (!useLayeredAtmosphere) {
                // Стабилизация для плотных атмосфер (см. computeLocalGreenhouseOpacity):
                // t_eff недостаточна для сверхплотных атмосфер (например, Венера),
                // поэтому минимальная база зависит от давления и состава.
                const double profileBaseTemperatureKelvin =
                    (point.pressureAtm >= 0.1)
                        ? qMax(point.state.temperatureKelvin(),
                               denseAtmosphereMinTemperatureKelvin(effectiveTemperatureKelvin,
                                                                   point.pressureAtm,
                                                                   point.state.greenhouseOpacity(),
                                                                   co2Share,
                                                                   currentMinDenseAtmosphereTemperatureKelvin()))
                        : point.state.temperatureKelvin();
                const auto radiationModel =
                    makeRadiationModel(atmosphere,
                                       point.pressureAtm,
                                       profileBaseTemperatureKelvin,
                                       effectiveTemperatureKelvin,
                                       gravity,
                                       radiationModelType);
                point.airTemperatureK =
                    resolveAirTemperatureKelvin(point.state,
                                                point.pressureAtm,
                                                gravity,
                                                initialAirTemperature,
                                                localInsolation,
                                                cloudShortwaveTransmission,
                                                radiationModelType,
                                                *radiationModel,
                                                useAtmosphericModel,
                                                manualGreenhouseOnTop,
                                                timeStepSeconds);
            }
            if (logDetails) {
                qCInfo(solarRadiationLog) << "Resolved air temperature (tick)"
                                          << "index=" << i
                                          << "airTemperatureK=" << point.airTemperatureK;
            }
            if (!useLayeredAtmosphere) {
                const double roughnessLengthMeters = material.roughnessLengthMeters;
                const double safePressureAtm = qMax(0.0, point.pressureAtm);
                if (safePressureAtm > 0.0 && gravity > 0.0) {
                    const double columnMassKgPerM2 =
                        (safePressureAtm * kStandardPressurePa) / gravity;
                    const double airHeatCapacity =
                        columnMassKgPerM2 * kDryAirSpecificHeatJPerKgK;
                    if (airHeatCapacity > 0.0) {
                        AtmosphericCellState airState(point.airTemperatureK, airHeatCapacity);
                        const double windSpeed =
                            std::hypot(point.windEastMps, point.windNorthMps);
                        const double opticalDepth =
                            opticalDepthFromGreenhouseOpacity(point.state.greenhouseOpacity(),
                                                              radiationModelType);
                        const double longwaveEmissivity =
                            qBound(0.0, 1.0 - std::exp(-opticalDepth), 1.0);
                        SurfaceAtmosphereCoupler coupler(kDefaultHeatTransferWPerM2K);
                        double surfaceAirFluxWPerM2 = 0.0;
                        coupler.exchangeHeat(point.state,
                                             airState,
                                             windSpeed,
                                             longwaveEmissivity,
                                             roughnessLengthMeters,
                                             timeStepSeconds,
                                             &surfaceAirFluxWPerM2);
                        point.airTemperatureK = airState.airTemperatureKelvin();
                        point.temperatureK = point.state.temperatureKelvin();
                        point.surfaceAirFluxWPerM2 = surfaceAirFluxWPerM2;
                        // Поток в грунт учитывает радиацию и обмен с воздухом.
                        point.subsurfaceFluxWPerM2 += -surfaceAirFluxWPerM2;
                    }
                }
            }
            minTemperature = qMin(minTemperature, point.temperatureK);
            maxTemperature = qMax(maxTemperature, point.temperatureK);
            minAirTemperature = qMin(minAirTemperature, point.airTemperatureK);
            maxAirTemperature = qMax(maxAirTemperature, point.airTemperatureK);
        }

        // Обновляем карту после каждого тика таймера, чтобы сразу отражать новую температуру.
        // Также обновляем виджеты, чтобы в них попали обновлённые поверхностные
        // температура и давление после шага переноса.
        applySurfaceGridToViews();
        updateSurfaceStarRendering(declinationDegrees,
                                   elapsedTimeDays,
                                   siderealDayLengthDays,
                                   rotationMode,
                                   orbitalPhaseRadians,
                                   spinOrbitP,
                                   spinOrbitQ,
                                   isRetrograde,
                                   distanceAU);
        if (minTemperature <= maxTemperature) {
            applySurfaceTemperatureRangeToViews(minTemperature, maxTemperature);
            updateSurfaceTemperatureLegend(true, minTemperature, maxTemperature);
        } else {
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
        }
        if (minAirTemperature <= maxAirTemperature) {
            if (surfaceMapMode_ == SurfaceMapMode::AirTemperature) {
                applySurfaceTemperatureRangeToViews(minAirTemperature, maxAirTemperature);
            }
            updateSurfaceAirTemperatureLegend(true, minAirTemperature, maxAirTemperature);
        } else {
            updateSurfaceAirTemperatureLegend(false, 0.0, 0.0);
        }
        if (minPressureAtm <= maxPressureAtm) {
            applySurfacePressureRangeToViews(minPressureAtm, maxPressureAtm);
            updateSurfacePressureLegend(true, minPressureAtm, maxPressureAtm);
        } else {
            updateSurfacePressureLegend(false, 0.0, 0.0);
        }

        const bool publishDaily = surfaceTime_.hourIndex + 1 >= solarStepsPerDay;
        updateSurfaceTemperatureAggregation(publishDaily, rotationMode, useAtmosphericModel);

        advanceSurfaceTimeFlow();
    }

    void updateSurfaceTemperatureLegend(bool hasRange, double minTemperature, double maxTemperature) {
        hasSurfaceTemperatureRange_ = hasRange;
        if (hasRange) {
            surfaceMinTemperatureK_ = minTemperature;
            surfaceMaxTemperatureK_ = maxTemperature;
        }
        if (surfaceMapMode_ == SurfaceMapMode::Temperature) {
            refreshSurfaceLegend();
        }
    }

    void updateSurfaceAirTemperatureLegend(bool hasRange,
                                           double minTemperature,
                                           double maxTemperature) {
        hasSurfaceAirTemperatureRange_ = hasRange;
        if (hasRange) {
            surfaceMinAirTemperatureK_ = minTemperature;
            surfaceMaxAirTemperatureK_ = maxTemperature;
        }
        if (surfaceMapMode_ == SurfaceMapMode::AirTemperature) {
            refreshSurfaceLegend();
        }
    }

    void updateSurfaceWindLegend(bool hasRange, double minWindSpeed, double maxWindSpeed) {
        hasSurfaceWindRange_ = hasRange;
        if (hasRange) {
            surfaceMinWindSpeedMps_ = minWindSpeed;
            surfaceMaxWindSpeedMps_ = maxWindSpeed;
        }
        if (surfaceMapMode_ == SurfaceMapMode::Wind) {
            refreshSurfaceLegend();
        }
    }

    void updateSurfacePressureLegend(bool hasRange, double minPressureAtm, double maxPressureAtm) {
        hasSurfacePressureRange_ = hasRange;
        if (hasRange) {
            surfaceMinPressureAtm_ = minPressureAtm;
            surfaceMaxPressureAtm_ = maxPressureAtm;
        }
        if (surfaceMapMode_ == SurfaceMapMode::Pressure) {
            refreshSurfaceLegend();
        }
    }

    void updateSurfaceHeightLegendFromGrid() {
        if (surfaceGrid_.points().isEmpty()) {
            hasSurfaceHeightRange_ = false;
            surfaceMinHeightKm_ = 0.0;
            surfaceMaxHeightKm_ = 0.0;
            if (surfaceMapMode_ == SurfaceMapMode::Height) {
                refreshSurfaceLegend();
            }
            return;
        }

        double minHeight = surfaceGrid_.points().first().heightKm;
        double maxHeight = minHeight;
        for (const auto &point : surfaceGrid_.points()) {
            minHeight = qMin(minHeight, point.heightKm);
            maxHeight = qMax(maxHeight, point.heightKm);
        }

        hasSurfaceHeightRange_ = minHeight <= maxHeight;
        surfaceMinHeightKm_ = minHeight;
        surfaceMaxHeightKm_ = maxHeight;
        if (surfaceMapMode_ == SurfaceMapMode::Height) {
            refreshSurfaceLegend();
        }
    }

    void refreshSurfaceLegend() {
        if (!surfaceMinTemperatureLabel_ || !surfaceMaxTemperatureLabel_) {
            return;
        }

        const QLocale locale;
        if (surfaceMapMode_ == SurfaceMapMode::Temperature) {
            if (surfaceLegendScaleStack_) {
                surfaceLegendScaleStack_->setCurrentWidget(temperatureScaleWidget_);
            }
            if (!hasSurfaceTemperatureRange_) {
                surfaceMinTemperatureLabel_->setText(QStringLiteral("Мин: —"));
                surfaceMaxTemperatureLabel_->setText(QStringLiteral("Макс: —"));
                if (temperatureScaleWidget_) {
                    temperatureScaleWidget_->clearRange();
                }
                return;
            }

            surfaceMinTemperatureLabel_->setText(
                QStringLiteral("Мин: %1 K").arg(locale.toString(surfaceMinTemperatureK_, 'f', 1)));
            surfaceMaxTemperatureLabel_->setText(
                QStringLiteral("Макс: %1 K").arg(locale.toString(surfaceMaxTemperatureK_, 'f', 1)));
            if (temperatureScaleWidget_) {
                temperatureScaleWidget_->setTemperatureRange(surfaceMinTemperatureK_,
                                                             surfaceMaxTemperatureK_);
            }
            return;
        }

        if (surfaceMapMode_ == SurfaceMapMode::AirTemperature) {
            if (surfaceLegendScaleStack_) {
                surfaceLegendScaleStack_->setCurrentWidget(temperatureScaleWidget_);
            }
            if (!hasSurfaceAirTemperatureRange_) {
                surfaceMinTemperatureLabel_->setText(QStringLiteral("Мин: —"));
                surfaceMaxTemperatureLabel_->setText(QStringLiteral("Макс: —"));
                if (temperatureScaleWidget_) {
                    temperatureScaleWidget_->clearRange();
                }
                return;
            }

            surfaceMinTemperatureLabel_->setText(
                QStringLiteral("Мин: %1 K").arg(locale.toString(surfaceMinAirTemperatureK_, 'f', 1)));
            surfaceMaxTemperatureLabel_->setText(
                QStringLiteral("Макс: %1 K").arg(locale.toString(surfaceMaxAirTemperatureK_, 'f', 1)));
            if (temperatureScaleWidget_) {
                temperatureScaleWidget_->setTemperatureRange(surfaceMinAirTemperatureK_,
                                                             surfaceMaxAirTemperatureK_);
            }
            return;
        }

        if (surfaceMapMode_ == SurfaceMapMode::Height) {
            if (surfaceLegendScaleStack_) {
                surfaceLegendScaleStack_->setCurrentWidget(heightScaleWidget_);
            }
            if (!hasSurfaceHeightRange_) {
                surfaceMinTemperatureLabel_->setText(QStringLiteral("Мин: —"));
                surfaceMaxTemperatureLabel_->setText(QStringLiteral("Макс: —"));
                if (heightScaleWidget_) {
                    heightScaleWidget_->clearRange();
                }
                return;
            }

            surfaceMinTemperatureLabel_->setText(
                QStringLiteral("Мин: %1 км").arg(locale.toString(surfaceMinHeightKm_, 'f', 1)));
            surfaceMaxTemperatureLabel_->setText(
                QStringLiteral("Макс: %1 км").arg(locale.toString(surfaceMaxHeightKm_, 'f', 1)));
            if (heightScaleWidget_) {
                heightScaleWidget_->setHeightRange(surfaceMinHeightKm_, surfaceMaxHeightKm_);
            }
            return;
        }

        if (surfaceMapMode_ == SurfaceMapMode::Wind) {
            if (surfaceLegendScaleStack_) {
                surfaceLegendScaleStack_->setCurrentWidget(windScaleWidget_);
            }
            if (!hasSurfaceWindRange_) {
                surfaceMinTemperatureLabel_->setText(QStringLiteral("Мин: —"));
                surfaceMaxTemperatureLabel_->setText(QStringLiteral("Макс: —"));
                if (windScaleWidget_) {
                    windScaleWidget_->clearRange();
                }
                return;
            }

            surfaceMinTemperatureLabel_->setText(
                QStringLiteral("Мин: %1 м/с").arg(locale.toString(surfaceMinWindSpeedMps_, 'f', 1)));
            surfaceMaxTemperatureLabel_->setText(
                QStringLiteral("Макс: %1 м/с").arg(locale.toString(surfaceMaxWindSpeedMps_, 'f', 1)));
            if (windScaleWidget_) {
                windScaleWidget_->setWindRange(surfaceMinWindSpeedMps_, surfaceMaxWindSpeedMps_);
            }
            return;
        }

        if (surfaceMapMode_ == SurfaceMapMode::Realistic) {
            if (surfaceLegendScaleStack_) {
                surfaceLegendScaleStack_->setCurrentWidget(realisticScaleWidget_);
            }
            surfaceMinTemperatureLabel_->setText(QStringLiteral("Мин: —"));
            surfaceMaxTemperatureLabel_->setText(QStringLiteral("Макс: —"));
            return;
        }

        if (surfaceLegendScaleStack_) {
            surfaceLegendScaleStack_->setCurrentWidget(pressureScaleWidget_);
        }
        if (!hasSurfacePressureRange_) {
            surfaceMinTemperatureLabel_->setText(QStringLiteral("Мин: —"));
            surfaceMaxTemperatureLabel_->setText(QStringLiteral("Макс: —"));
            if (pressureScaleWidget_) {
                pressureScaleWidget_->clearRange();
            }
            return;
        }

        surfaceMinTemperatureLabel_->setText(
            QStringLiteral("Мин: %1 атм").arg(locale.toString(surfaceMinPressureAtm_, 'f', 3)));
        surfaceMaxTemperatureLabel_->setText(
            QStringLiteral("Макс: %1 атм").arg(locale.toString(surfaceMaxPressureAtm_, 'f', 3)));
        if (pressureScaleWidget_) {
            pressureScaleWidget_->setPressureRange(surfaceMinPressureAtm_, surfaceMaxPressureAtm_);
        }
    }

    QProgressDialog *ensureTemperatureProgressDialog() {
        if (temperatureProgressDialog_) {
            return temperatureProgressDialog_;
        }

        temperatureProgressDialog_ = new QProgressDialog(this);
        temperatureProgressDialog_->setWindowTitle(QStringLiteral("Расчет температур"));
        temperatureProgressDialog_->setLabelText(QStringLiteral("Вычисление температурного профиля..."));
        temperatureProgressDialog_->setCancelButtonText(QStringLiteral("Отмена"));
        temperatureProgressDialog_->setWindowModality(Qt::WindowModal);
        temperatureProgressDialog_->setAutoClose(true);
        temperatureProgressDialog_->setAutoReset(true);
        // Диалог создается по требованию, чтобы не всплывать без запуска вычислений.
        temperatureProgressDialog_->hide();

        connect(temperatureProgressDialog_, &QProgressDialog::canceled, this, [this]() {
            cancelTemperatureCalculation();
        });

        auto *pauseButton = new QPushButton(QStringLiteral("Пауза"), temperatureProgressDialog_);
        pauseButton->setObjectName(QStringLiteral("temperaturePauseButton"));
        connect(pauseButton, &QPushButton::clicked, this, [this, pauseButton]() {
            const bool shouldPause = !temperaturePauseFlag_.load();
            temperaturePauseFlag_.store(shouldPause);
            updateTemperaturePauseUi(shouldPause);
            if (surfaceSimTimer_) {
                if (shouldPause) {
                    surfaceSimTimer_->stop();
                } else if (surfaceSimRunning_) {
                    surfaceSimTimer_->start();
                }
            }
        });
        if (auto *layout = temperatureProgressDialog_->layout()) {
            layout->addWidget(pauseButton);
        }

        return temperatureProgressDialog_;
    }

    void resetTemperatureUiState() {
        temperaturePauseFlag_.store(false);
        if (temperatureUiTimer_) {
            temperatureUiTimer_->stop();
        }
        updateTemperaturePauseUi(false);
    }

    void startTemperatureElapsedUi(int requestId, const QPointer<QProgressDialog> &dialogGuard) {
        if (temperatureUiTimer_ && temperatureUiTimer_->isActive() &&
            requestId == temperatureRequestId_) {
            return;
        }

        temperatureElapsed_.start();

        if (!temperatureUiTimer_) {
            temperatureUiTimer_ = new QTimer(this);
            temperatureUiTimer_->setInterval(1000);
        }
        temperatureUiTimer_->stop();
        temperatureUiTimer_->disconnect();
        auto updateElapsedLabel = [this, dialogGuard, requestId]() {
            if (requestId != temperatureRequestId_) {
                return;
            }
            const qint64 elapsedSeconds = temperatureElapsed_.elapsed() / 1000;
            const int minutes = static_cast<int>(elapsedSeconds / 60);
            const int seconds = static_cast<int>(elapsedSeconds % 60);
            const QString pauseSuffix =
                temperaturePauseFlag_.load() ? QStringLiteral(" (пауза)") : QString();
            if (dialogGuard) {
                dialogGuard->setLabelText(
                    QStringLiteral("Вычисление температурного профиля... %1:%2%3")
                        .arg(minutes, 2, 10, QLatin1Char('0'))
                        .arg(seconds, 2, 10, QLatin1Char('0'))
                        .arg(pauseSuffix));
            }
        };
        connect(temperatureUiTimer_, &QTimer::timeout, this, updateElapsedLabel);
        updateElapsedLabel();
        temperatureUiTimer_->start();
    }

    void updateTemperaturePlot() {
        cancelTemperatureCalculation();
        ensureTemperaturePlotSourceSurfaceGrid();

        // График температуры временно отключён: основная модель теперь "Поверхность" с ячейками.
        // Оставляем только очистку сегментов, чтобы UI не отображал устаревшие кривые.
        clearTemperatureSegments();
    }

    void cancelTemperatureCalculation() {
        if (temperatureCancelFlag_) {
            temperatureCancelFlag_->store(true);
        }
        resetTemperatureUiState();
        if (temperatureProgressDialog_) {
            temperatureProgressDialog_->hide();
        }
    }

    void resetSolarConstant() {
        hasSolarConstant_ = false;
        lastSolarConstant_ = 0.0;
        lastSolarConstantDistanceAU_ = 0.0;
        resultLabel_->setText(
            QStringLiteral("Введите параметры и нажмите \"Рассчитать\"."));
        updateTemperaturePlot();
        updateSurfaceSolarConstantLabel(SurfaceFluxBreakdown{});
    }
};
}  // namespace

ArgumentsParseResult parseParametersFromArguments(const QCoreApplication &app,
                                                  QTextStream &output,
                                                  BinarySystemParameters &parameters,
                                                  int &precision,
                                                  bool &enableRadiationLog) {
    QCommandLineParser parser;
    parser.setApplicationDescription(
        QStringLiteral("Вычисление солнечной постоянной по параметрам звезды."));
    parser.addHelpOption();

    QCommandLineOption radiusOption({QStringLiteral("r"), QStringLiteral("radius")},
                                    QStringLiteral("Радиус звезды в солнечных радиусах."),
                                    QStringLiteral("value"));
    QCommandLineOption temperatureOption({QStringLiteral("t"), QStringLiteral("temperature")},
                                         QStringLiteral("Температура поверхности в К."),
                                         QStringLiteral("value"));
    QCommandLineOption distanceOption({QStringLiteral("d"), QStringLiteral("distance")},
                                      QStringLiteral("Расстояние от барицентра до планеты в а.е."),
                                      QStringLiteral("value"));
    QCommandLineOption radius2Option({QStringLiteral("r2"), QStringLiteral("radius2")},
                                     QStringLiteral("Радиус второй звезды в солнечных радиусах."),
                                     QStringLiteral("value"));
    QCommandLineOption temperature2Option({QStringLiteral("t2"), QStringLiteral("temperature2")},
                                          QStringLiteral("Температура второй звезды в К."),
                                          QStringLiteral("value"));
    QCommandLineOption distance2Option({QStringLiteral("d2"), QStringLiteral("distance2")},
                                       QStringLiteral(
                                           "Расстояние от барицентра до планеты для второй звезды в а.е."),
                                       QStringLiteral("value"));
    QCommandLineOption precisionOption({QStringLiteral("p"), QStringLiteral("precision")},
                                       QStringLiteral("Количество значащих цифр в выводе."),
                                       QStringLiteral("digits"),
                                       QString::number(kDefaultPrecision));
    QCommandLineOption radiationLogOption(QStringLiteral("radiation-log"),
                                          QStringLiteral("Включить подробные логи радиационной модели."));

    parser.addOptions({radiusOption, temperatureOption, distanceOption, radius2Option, temperature2Option,
                       distance2Option, precisionOption, radiationLogOption});
    parser.process(app);
    if (parser.isSet(radiationLogOption)) {
        enableRadiationLog = true;
    }

    const auto parsePositive = [&](const QCommandLineOption &option, double &target,
                                   const QString &name) -> bool {
        const QString valueString = parser.value(option);
        bool ok = false;
        const double parsed = valueString.toDouble(&ok);
        if (!ok || !std::isfinite(parsed) || parsed <= 0.0) {
            output << "Значение для " << name << " должно быть положительным числом." << Qt::endl;
            return false;
        }
        target = parsed;
        return true;
    };

    const bool primaryParametersProvided = parser.isSet(radiusOption) || parser.isSet(temperatureOption) ||
                                           parser.isSet(distanceOption);
    const bool secondaryParametersProvided = parser.isSet(radius2Option) || parser.isSet(temperature2Option) ||
                                             parser.isSet(distance2Option);
    const bool precisionProvided = parser.isSet(precisionOption);

    if (precisionProvided) {
        bool ok = false;
        const int parsedPrecision = parser.value(precisionOption).toInt(&ok);
        if (!ok || parsedPrecision <= 0 || parsedPrecision > 15) {
            output << "Точность вывода должна быть целым числом от 1 до 15." << Qt::endl;
            return ArgumentsParseResult::Failure;
        }
        precision = parsedPrecision;
    }

    if (!primaryParametersProvided && !secondaryParametersProvided) {
        return ArgumentsParseResult::None;
    }

    if (!primaryParametersProvided) {
        output << "Сначала укажите параметры основной звезды: --radius, --temperature и --distance." << Qt::endl;
        return ArgumentsParseResult::Failure;
    }

    if (!parser.isSet(radiusOption) || !parser.isSet(temperatureOption) || !parser.isSet(distanceOption)) {
        output << "Для неинтерактивного режима укажите все параметры первой звезды: --radius, --temperature и --distance."
               << Qt::endl;
        return ArgumentsParseResult::Failure;
    }

    if (secondaryParametersProvided &&
        (!parser.isSet(radius2Option) || !parser.isSet(temperature2Option) || !parser.isSet(distance2Option))) {
        output << "Если указываете вторую звезду, задайте сразу три параметра: --radius2, --temperature2 и --distance2."
               << Qt::endl;
        return ArgumentsParseResult::Failure;
    }

    if (!parsePositive(radiusOption, parameters.primary.radiusInSolarRadii, QStringLiteral("radius")) ||
        !parsePositive(temperatureOption, parameters.primary.temperatureKelvin, QStringLiteral("temperature")) ||
        !parsePositive(distanceOption, parameters.primary.distanceInAU, QStringLiteral("distance"))) {
        return ArgumentsParseResult::Failure;
    }

    if (secondaryParametersProvided) {
        StellarParameters secondary{};
        if (!parsePositive(radius2Option, secondary.radiusInSolarRadii, QStringLiteral("radius2")) ||
            !parsePositive(temperature2Option, secondary.temperatureKelvin, QStringLiteral("temperature2")) ||
            !parsePositive(distance2Option, secondary.distanceInAU, QStringLiteral("distance2"))) {
            return ArgumentsParseResult::Failure;
        }
        parameters.secondary = secondary;
    }

    return ArgumentsParseResult::Success;
}

void promptAndComputeSolarConstant(QTextStream &input, QTextStream &output, int precision) {
    output.setRealNumberPrecision(precision);
    output.setRealNumberNotation(QTextStream::SmartNotation);

    BinarySystemParameters parameters{};

    output << "Введите радиус звезды (в солнечных радиусах):" << Qt::endl;
    output.flush();
    input >> parameters.primary.radiusInSolarRadii;
    if (input.status() != QTextStream::Ok || parameters.primary.radiusInSolarRadii <= 0.0) {
        output << "Радиус должен быть положительным числом." << Qt::endl;
        return;
    }

    output << "Введите температуру поверхности (в К):" << Qt::endl;
    output.flush();
    input >> parameters.primary.temperatureKelvin;
    if (input.status() != QTextStream::Ok || parameters.primary.temperatureKelvin <= 0.0) {
        output << "Температура должна быть положительным числом." << Qt::endl;
        return;
    }

    output << "Введите расстояние от барицентра до планеты (в а.е.):" << Qt::endl;
    output.flush();
    input >> parameters.primary.distanceInAU;
    if (input.status() != QTextStream::Ok || parameters.primary.distanceInAU <= 0.0) {
        output << "Расстояние от барицентра должно быть положительным числом." << Qt::endl;
        return;
    }

    output << "Есть вторая звезда в системе? (y/n):" << Qt::endl;
    output.flush();
    QString answer;
    input >> answer;

    const bool hasSecondary = answer.compare(QStringLiteral("y"), Qt::CaseInsensitive) == 0 ||
                              answer.compare(QStringLiteral("yes"), Qt::CaseInsensitive) == 0 ||
                              answer.compare(QStringLiteral("д"), Qt::CaseInsensitive) == 0;

    if (hasSecondary) {
        parameters.secondary = StellarParameters{};

        output << "Введите радиус второй звезды (в солнечных радиусах):" << Qt::endl;
        output.flush();
        input >> parameters.secondary->radiusInSolarRadii;
        if (input.status() != QTextStream::Ok || parameters.secondary->radiusInSolarRadii <= 0.0) {
            output << "Радиус второй звезды должен быть положительным числом." << Qt::endl;
            return;
        }

        output << "Введите температуру поверхности второй звезды (в К):" << Qt::endl;
        output.flush();
        input >> parameters.secondary->temperatureKelvin;
        if (input.status() != QTextStream::Ok || parameters.secondary->temperatureKelvin <= 0.0) {
            output << "Температура второй звезды должна быть положительным числом." << Qt::endl;
            return;
        }

        output << "Введите расстояние от барицентра до планеты для второй звезды (в а.е.):" << Qt::endl;
        output.flush();
        input >> parameters.secondary->distanceInAU;
        if (input.status() != QTextStream::Ok || parameters.secondary->distanceInAU <= 0.0) {
            output << "Расстояние от барицентра для второй звезды должно быть положительным числом." << Qt::endl;
            return;
        }
    }

    const double primaryFlux = SolarCalculator::solarConstant(parameters.primary);
    double totalFlux = primaryFlux;

    if (parameters.secondary) {
        const double secondaryFlux = SolarCalculator::solarConstant(*parameters.secondary);
        totalFlux += secondaryFlux;
        output << "Солнечная постоянная у планеты: " << totalFlux << " Вт/м²"
               << " (первая: " << primaryFlux << " Вт/м², вторая: "
               << secondaryFlux << " Вт/м²)" << Qt::endl;
        return;
    }

    output << "Солнечная постоянная у планеты: " << totalFlux << " Вт/м²" << Qt::endl;
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QTextStream output(stdout);
    QTextStream input(stdin);

    BinarySystemParameters parameters{};
    int precision = kDefaultPrecision;
    bool enableRadiationLog = isSolarRadiationLoggingEnabledFromEnvironment();
    const ArgumentsParseResult argsResult =
        parseParametersFromArguments(app, output, parameters, precision, enableRadiationLog);
    if (enableRadiationLog) {
        enableSolarRadiationLogging();
    }

    switch (argsResult) {
    case ArgumentsParseResult::Failure:
        return 1;
    case ArgumentsParseResult::Success: {
        output.setRealNumberPrecision(precision);
        output.setRealNumberNotation(QTextStream::SmartNotation);
        const double flux = SolarCalculator::solarConstant(parameters);
        output << "Солнечная постоянная у планеты: " << flux << " Вт/м²" << Qt::endl;
        return 0;
    }
    case ArgumentsParseResult::None:
        break;
    }

    SolarCalculatorWidget widget(precision, enableRadiationLog);
    widget.setWindowTitle(QStringLiteral("Калькулятор солнечной постоянной"));
    widget.show();

    return app.exec();
}
