#include "mode_illustration_widget.h"
#include "orbit_animation_model.h"
#include "planet_presets.h"
#include "segment_selector_widget.h"
#include "solar_calculator.h"
#include "solar_display.h"
#include "atmosphere_widget.h"
#include "surface_tile_temperature_calculator.h"
#include "surface_temperature_plot.h"
#include "surface_map_widget.h"
#include "surface_globe_widget.h"
#include "surface_point_status_dialog.h"
#include "surface_atmosphere_coupler.h"
#include "surface_temperature_scale_widget.h"
#include "surface_height_scale_widget.h"
#include "surface_wind_scale_widget.h"
#include "surface_pressure_scale_widget.h"
#include "surface_map_mode.h"
#include "surface_advection_model.h"
#include "surface_pressure_transport_model.h"
#include "wind_field_model.h"
#include "planet_surface_grid.h"
#include "subsurface_temperature_solver.h"
#include "radiation_model.h"
#include "layered_radiation_model.h"
#include "radiation_model_utils.h"
#include "atmospheric_pressure_model.h"
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
#include <QtConcurrent/QtConcurrent>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QProgressDialog>
#include <QtWidgets/QSpinBox>
#include <algorithm>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QStyle>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

#include <atomic>
#include <cmath>
#include <functional>
#include <limits>

namespace {
Q_LOGGING_CATEGORY(solarRadiationLog, "solar.radiation")

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
    QString rules = QString::fromLocal8Bit(qgetenv("QT_LOGGING_RULES"));
    if (!rules.isEmpty()) {
        rules.append('\n');
    }
    rules.append(QStringLiteral("solar.radiation.info=true\nsolar.radiation.debug=true"));
    QLoggingCategory::setFilterRules(rules);
}

constexpr int kRoleSemiMajorAxis = Qt::UserRole;
constexpr int kRoleIsCustom = Qt::UserRole + 1;
constexpr int kRolePlanetName = Qt::UserRole + 2;
constexpr int kRoleMaterialId = Qt::UserRole + 3;
constexpr int kRoleDayLength = Qt::UserRole + 4;
constexpr int kRoleEccentricity = Qt::UserRole + 5;
constexpr int kRoleObliquity = Qt::UserRole + 6;
constexpr int kRolePerihelionArgument = Qt::UserRole + 7;
constexpr int kRoleMassEarths = Qt::UserRole + 8;
constexpr int kRoleRadiusKm = Qt::UserRole + 9;
constexpr int kRoleRotationMode = Qt::UserRole + 10;
constexpr int kRoleAtmosphere = Qt::UserRole + 11;
constexpr int kRoleGreenhouseOpacity = Qt::UserRole + 12;
constexpr int kRoleCloudAlbedo = Qt::UserRole + 13;
constexpr int kRoleHeightSourceType = Qt::UserRole + 14;
constexpr int kRoleHeightmapPath = Qt::UserRole + 15;
constexpr int kRoleHeightmapScaleKm = Qt::UserRole + 16;
constexpr int kRoleHeightSeed = Qt::UserRole + 17;
constexpr int kRoleUseContinentsHeight = Qt::UserRole + 18;
constexpr int kRoleHasSeaLevel = Qt::UserRole + 19;
constexpr int kRoleFlatHeight = Qt::UserRole + 20;
constexpr int kRoleManualGreenhouseOnTopOfAtmosphere = Qt::UserRole + 21;
constexpr int kRoleAdvancedRadiationModel = Qt::UserRole + 22;
constexpr int kRoleGeothermalFlux = Qt::UserRole + 23;
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

double estimateSurfaceWaterGigatons(const SurfaceMaterial &material) {
    // В текущем UI нет явного управления гидросферой, поэтому применяем мягкую эвристику:
    // океаническую поверхность считаем водной, остальные материалы — сухими.
    if (material.id == QLatin1String("ocean")) {
        return kEarthWaterGigatons;
    }
    return 0.0;
}

double computeLocalGreenhouseOpacity(const AtmosphereComposition &atmosphere,
                                     const SurfaceMaterial &material,
                                     double pressureAtm,
                                     double planetRadiusKm,
                                     double surfaceGravity,
                                     double blendedInsolation,
                                     double manualGreenhouseOpacity,
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

    // Давление также усиливает облачность: давление выше порога повышает отражение.
    const double pressureClouds =
        pressureAtm > 0.05 ? 0.25 * (1.0 - std::exp(-pressureAtm)) : 0.0;
    const double albedo = qBound(0.0, material.albedo, 1.0);
    const double surfAlbedoPre =
        (1.0 - potentialCoverage) * albedo + potentialCoverage * 0.06;
    const double tEffPre =
        std::pow((blendedInsolation * (1.0 - qMax(surfAlbedoPre, pressureClouds))) /
                     (4.0 * kStefanBoltzmannConstant),
                 0.25);
    const auto preRadiationModel = makeRadiationModel(atmosphere,
                                                      pressureAtm,
                                                      tEffPre,
                                                      tEffPre,
                                                      surfaceGravity,
                                                      radiationModelType);
    const double baseLongwaveTransmission =
        qMax(1e-6, preRadiationModel->outgoingTransmission());
    const double tBasePre =
        tEffPre * std::pow(1.0 / baseLongwaveTransmission, 0.25);

    double evaporation = 0.0;
    if (potentialCoverage > 0.0 && tBasePre > 263.0) {
        evaporation = potentialCoverage * std::exp((tBasePre - 280.0) / 15.0);
    }
    const double waterTau = qMin(8.0, evaporation * 1.5);
    double extraTau = waterTau;
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

class SolarCalculatorWidget : public QWidget {
public:
    explicit SolarCalculatorWidget(int precision, QWidget *parent = nullptr)
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

        addPresetButton(QStringLiteral("Солнце"), [this, applyPrimary, applySecondary]() {
            applyPrimary(StellarParameters{1.0, 5772.0, 1.0});
            applySecondary(std::nullopt);
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

        addPresetButton(QStringLiteral("Сладкое Небо"), [this, applyPrimary, applySecondary]() {
            applyPrimary(StellarParameters{0.3761, 2576.0, 1.0});
            applySecondary(StellarParameters{0.3741, 2349.0, 1.0});
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
        heightSeedSpinBox_ = new QSpinBox(this);
        heightSeedSpinBox_->setRange(0, std::numeric_limits<int>::max());
        flatHeightButton_ = new QPushButton(QStringLiteral("Обнулить высоту"), this);
        flatHeightButton_->setCheckable(true);
        flatHeightButton_->setEnabled(false);
        updateFlatHeightButtonText(false);
        cloudAlbedoSpinBox_ = new QDoubleSpinBox(this);
        cloudAlbedoSpinBox_->setRange(0.0, 1.0);
        cloudAlbedoSpinBox_->setDecimals(2);
        cloudAlbedoSpinBox_->setSingleStep(0.05);
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
        modeIllustrationWidget_ = new ModeIllustrationWidget(this);
        modeIllustrationWidget_->setRotationMode(
            static_cast<RotationMode>(rotationModeComboBox_->currentData().toInt()));
        auto *latitudeStepLabel = new QLabel(QStringLiteral("Через 1° (быстрые)"), this);
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
        auto *latitudeStepWidget = new QWidget(this);
        auto *latitudeStepLayout = new QHBoxLayout();
        latitudeStepLayout->addWidget(latitudeStepLabel);
        latitudeStepLayout->addStretch();
        latitudeStepLayout->setContentsMargins(0, 0, 0, 0);
        latitudeStepWidget->setLayout(latitudeStepLayout);

        auto *planetColumnsLayout = new QHBoxLayout();
        planetColumnsLayout->addLayout(planetLeftFormLayout);
        planetColumnsLayout->addLayout(planetRightFormLayout);

        auto *planetControlsLayout = new QFormLayout();
        planetControlsLayout->addRow(QStringLiteral("Материал поверхности:"), materialComboBox_);
        planetControlsLayout->addRow(QStringLiteral("Режим вращения:"), rotationModeWidget);
        planetControlsLayout->addRow(QStringLiteral("Семя рельефа:"), heightSeedSpinBox_);
        planetControlsLayout->addRow(QStringLiteral("Высота поверхности:"), flatHeightButton_);
        planetControlsLayout->addRow(QStringLiteral("Альбедо облаков (0..1):"), cloudAlbedoSpinBox_);
        planetControlsLayout->addRow(QString(), manualGreenhouseOnTopCheckBox_);
        planetControlsLayout->addRow(QString(), advancedRadiationCheckBox_);
        planetControlsLayout->addRow(QStringLiteral("Шаг по широте:"), latitudeStepWidget);
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

        auto *starsPanelLayout = new QVBoxLayout();
        starsPanelLayout->addWidget(primaryGroupBox_);
        starsPanelLayout->addWidget(secondStarCheckBox_);
        starsPanelLayout->addWidget(secondaryGroupBox_);
        auto *starsPanel = new QGroupBox(QStringLiteral("Панель звезд"), this);
        starsPanel->setLayout(starsPanelLayout);

        temperaturePlot_ = new SurfaceTemperaturePlot(this);
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
                    surfacePointStatusDialog_->setPoint(*point);
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
        temperaturePauseButton_ = new QPushButton(QStringLiteral("Пауза"), this);
        surfaceSimToggleButton_ = new QPushButton(QStringLiteral("Старт"), this);
        surfaceSimSpeedComboBox_ = new QComboBox(this);
        surfaceSimSpeedComboBox_->addItem(QStringLiteral("1x"), 1.0);
        surfaceSimSpeedComboBox_->addItem(QStringLiteral("10x"), 10.0);
        surfaceSimTimeLabel_ = new QLabel(QStringLiteral("t = —"), this);
        temperatureElapsedLabel_ = new QLabel(QStringLiteral("Прошло: 00:00"), this);
        surfaceSeamlessCheckBox_ = new QCheckBox(QStringLiteral("Бесшовная карта"), this);
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
        surfaceLegendScaleStack_ = new QStackedWidget(this);
        surfaceLegendScaleStack_->addWidget(temperatureScaleWidget_);
        surfaceLegendScaleStack_->addWidget(heightScaleWidget_);
        surfaceLegendScaleStack_->addWidget(windScaleWidget_);
        surfaceLegendScaleStack_->addWidget(pressureScaleWidget_);
        surfaceLegendScaleStack_->setCurrentWidget(temperatureScaleWidget_);
        auto *surfaceLegendTopLayout = new QHBoxLayout();
        surfaceLegendTopLayout->addWidget(surfaceMinTemperatureLabel_);
        surfaceLegendTopLayout->addStretch();
        surfaceLegendTopLayout->addWidget(surfaceMaxTemperatureLabel_);
        auto *surfaceControlLayout = new QHBoxLayout();
        // Управление дублирует диалог, чтобы расчётом можно было управлять без модального окна.
        surfaceControlLayout->addWidget(temperaturePauseButton_);
        surfaceControlLayout->addWidget(surfaceSimToggleButton_);
        surfaceControlLayout->addWidget(surfaceSimSpeedComboBox_);
        surfaceControlLayout->addWidget(surfaceSimTimeLabel_);
        surfaceControlLayout->addWidget(temperatureElapsedLabel_);
        surfaceControlLayout->addWidget(surfaceSeamlessCheckBox_);
        surfaceControlLayout->addWidget(new QLabel(QStringLiteral("Карта:"), this));
        surfaceControlLayout->addWidget(surfaceMapModeComboBox_);
        surfaceControlLayout->addWidget(surfaceViewToggleButton_);
        surfaceControlLayout->addWidget(surfaceMarkupCheckBox_);
        surfaceControlLayout->addStretch();
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
        auto *surfaceMapLayout = new QVBoxLayout();
        surfaceMapLayout->addLayout(surfaceLegendTopLayout);
        surfaceMapLayout->addLayout(surfaceControlLayout);
        surfaceMapLayout->addWidget(subsurfaceGroupBox);
        surfaceMapLayout->addWidget(surfaceViewStack_, 1);
        surfaceMapLayout->addLayout(surfaceLegendBottomLayout);
        auto *surfaceMapContainer = new QWidget(this);
        surfaceMapContainer->setLayout(surfaceMapLayout);
        auto *plotGroupBox = new QGroupBox(QStringLiteral("Температурный профиль"), this);
        auto *plotLayout = new QVBoxLayout(plotGroupBox);
        auto *segmentLayout = new QHBoxLayout();
        segmentSelectorWidget_ = new SegmentSelectorWidget(plotGroupBox);
        segmentSelectorWidget_->setEnabled(false);
        segmentLayout->addWidget(new QLabel(QStringLiteral("Сегмент орбиты:"), plotGroupBox));
        segmentLayout->addWidget(segmentSelectorWidget_, 1);
        plotLayout->addLayout(segmentLayout);
        // auto *smoothingCheckBox = new QCheckBox(QStringLiteral("Сглаживать график"), plotGroupBox);
        // plotLayout->addWidget(smoothingCheckBox);
        plotLayout->addWidget(temperaturePlot_);
        plotLayout->addWidget(modeIllustrationWidget_, 0, Qt::AlignHCenter);
        plotGroupBox->setLayout(plotLayout);
        plotGroupBox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);

        auto *leftLayout = new QVBoxLayout();
        leftLayout->addLayout(presetsLayout);
        leftLayout->addWidget(starsPanel);
        leftLayout->addWidget(planetGroupBox);
        leftLayout->addWidget(calculateButton);
        leftLayout->addStretch();

        auto *rightTabs = new QTabWidget(this);
        rightTabs->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
        rightTabs->addTab(plotGroupBox, tr("Температура"));
        rightTabs->addTab(atmosphereWidget_, tr("Атмосфера"));
        rightTabs->addTab(surfaceMapContainer, tr("Поверхность"));
        connect(rightTabs, &QTabWidget::currentChanged, this,
                [this, rightTabs, surfaceMapContainer](int) {
                    if (rightTabs->currentWidget() != surfaceMapContainer) {
                        return;
                    }
                    // Пересчитываем поверхность только при первом открытии вкладки.
                    if (surfaceGrid_.points().isEmpty()) {
                        updateSurfaceGridTemperatures();
                    }
                });

        auto *rightLayout = new QVBoxLayout();
        rightLayout->addWidget(rightTabs, 1);

        auto *layout = new QHBoxLayout(this);
        layout->addLayout(leftLayout, 0);
        layout->addLayout(rightLayout, 1);

        setLayout(layout);
        resize(480, 360);
        // Оставляем один поток под UI, чтобы параллельные вычисления не блокировали интерфейс.
        QThreadPool::globalInstance()->setMaxThreadCount(
            qMax(1, QThread::idealThreadCount() - 1));

        connect(planetComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
                [this, rightTabs, surfaceMapContainer](int) {
            cancelTemperatureCalculation();
            updatePlanetSemiMajorAxisLabel();
            updatePlanetDayLengthLabel();
            updatePlanetMassLabel();
            updatePlanetRadiusLabel();
            updatePlanetDerivedLabels();
            updateAtmospherePlanetParameters();
            updateAtmosphereComposition();
            updatePlanetOrbitLabels();
            syncMaterialWithPlanet();
            if (rightTabs->currentWidget() == surfaceMapContainer) {
                updateSurfaceGridTemperatures();
            }
            syncRotationModeWithPlanet();
            syncHeightSeedWithPlanet();
            syncFlatHeightWithPlanet();
            syncCloudAlbedoWithPlanet();
            syncGeothermalFluxWithPlanet();
            syncManualGreenhouseOnTopWithPlanet();
            syncAdvancedRadiationWithPlanet();
            updatePlanetActions();
            if (autoCalculateEnabled_ && hasPrimaryInputs() &&
                (!secondStarCheckBox_->isChecked() || hasSecondaryInputs())) {
                onCalculateRequested();
            } else {
                updateTemperaturePlot();
            }
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
            updatePlanetMassLabel();
            updatePlanetRadiusLabel();
            updatePlanetDerivedLabels();
            updateAtmospherePlanetParameters();
            updateAtmosphereComposition();
            updatePlanetActions();
            updateTemperaturePlot();
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

        connect(heightSeedSpinBox_, QOverload<int>::of(&QSpinBox::valueChanged), this, [this](int) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetHeightSeedWithSelection();
            updateSurfaceGridTemperatures();
        });

        connect(flatHeightButton_, &QPushButton::toggled, this, [this](bool checked) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetFlatHeightWithSelection(checked);
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
                    surfaceSimSpeedMultiplier_ =
                        surfaceSimSpeedComboBox_->currentData().toDouble();
                    updateSurfaceSimulationTimerInterval();
                    updateSurfaceSimulationUi();
                });

        connect(surfaceSeamlessCheckBox_, &QCheckBox::toggled, this, [this](bool checked) {
            if (surfaceMapWidget_) {
                surfaceMapWidget_->setInterpolationEnabled(checked);
            }
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
                [onSubsurfaceChanged](int) { onSubsurfaceChanged(); });
        connect(subsurfaceTopThicknessSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [onSubsurfaceChanged](double) { onSubsurfaceChanged(); });
        connect(subsurfaceDepthSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
                [onSubsurfaceChanged](double) { onSubsurfaceChanged(); });
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
    }

private:
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
        updateTemperaturePlot();
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
    QLabel *planetMassLabel_ = nullptr;
    QLabel *planetRadiusLabel_ = nullptr;
    QLabel *planetSurfaceGravityLabel_ = nullptr;
    QLabel *planetSurfaceAreaLabel_ = nullptr;
    QLabel *planetEccentricityLabel_ = nullptr;
    QLabel *planetObliquityLabel_ = nullptr;
    QLabel *planetPerihelionArgumentLabel_ = nullptr;
    QComboBox *materialComboBox_ = nullptr;
    QComboBox *rotationModeComboBox_ = nullptr;
    QSpinBox *heightSeedSpinBox_ = nullptr;
    QPushButton *flatHeightButton_ = nullptr;
    QDoubleSpinBox *cloudAlbedoSpinBox_ = nullptr;
    QCheckBox *manualGreenhouseOnTopCheckBox_ = nullptr;
    QCheckBox *advancedRadiationCheckBox_ = nullptr;
    ModeIllustrationWidget *modeIllustrationWidget_ = nullptr;
    QPushButton *addPlanetButton_ = nullptr;
    QPushButton *deletePlanetButton_ = nullptr;
    AtmosphereWidget *atmosphereWidget_ = nullptr;

    QLabel *resultLabel_ = nullptr;
    SurfaceTemperaturePlot *temperaturePlot_ = nullptr;
    SurfaceMapWidget *surfaceMapWidget_ = nullptr;
    SurfaceGlobeWidget *surfaceGlobeWidget_ = nullptr;
    QPointer<SurfacePointStatusDialog> surfacePointStatusDialog_;
    int selectedSurfacePointIndex_ = -1;
    QStackedWidget *surfaceViewStack_ = nullptr;
    QLabel *surfaceMinTemperatureLabel_ = nullptr;
    QLabel *surfaceMaxTemperatureLabel_ = nullptr;
    QPushButton *temperaturePauseButton_ = nullptr;
    QPushButton *surfaceSimToggleButton_ = nullptr;
    QComboBox *surfaceSimSpeedComboBox_ = nullptr;
    QLabel *surfaceSimTimeLabel_ = nullptr;
    QLabel *temperatureElapsedLabel_ = nullptr;
    QCheckBox *surfaceSeamlessCheckBox_ = nullptr;
    QComboBox *surfaceMapModeComboBox_ = nullptr;
    QPushButton *surfaceViewToggleButton_ = nullptr;
    QCheckBox *surfaceMarkupCheckBox_ = nullptr;
    QSpinBox *subsurfaceLayersSpinBox_ = nullptr;
    QDoubleSpinBox *subsurfaceTopThicknessSpinBox_ = nullptr;
    QDoubleSpinBox *subsurfaceDepthSpinBox_ = nullptr;
    QComboBox *subsurfaceBoundaryComboBox_ = nullptr;
    QDoubleSpinBox *subsurfaceGeothermalFluxSpinBox_ = nullptr;
    SurfaceTemperatureScaleWidget *temperatureScaleWidget_ = nullptr;
    SurfaceHeightScaleWidget *heightScaleWidget_ = nullptr;
    SurfaceWindScaleWidget *windScaleWidget_ = nullptr;
    SurfacePressureScaleWidget *pressureScaleWidget_ = nullptr;
    QStackedWidget *surfaceLegendScaleStack_ = nullptr;
    SegmentSelectorWidget *segmentSelectorWidget_ = nullptr;
    QProgressDialog *temperatureProgressDialog_ = nullptr;
    QElapsedTimer temperatureElapsed_;
    QTimer *temperatureUiTimer_ = nullptr;
    QTimer *surfaceSimTimer_ = nullptr;
    PlanetSurfaceGrid surfaceGrid_;
    std::shared_ptr<std::atomic_bool> temperatureCancelFlag_;
    std::atomic_bool temperaturePauseFlag_{false};
    int temperatureRequestId_ = 0;
    int precision_ = kDefaultPrecision;
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
    std::optional<StellarCacheKey> lastStellarKey_;
    bool surfaceSimRunning_ = false;
    struct SurfaceSimulationState {
        int dayIndex = 0;
        int hourIndex = 0;
    } surfaceSimState_;
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

    void updateFlatHeightButtonText(bool useFlatHeight) {
        if (!flatHeightButton_) {
            return;
        }
        flatHeightButton_->setText(useFlatHeight ? QStringLiteral("Вернуть высоту")
                                                 : QStringLiteral("Обнулить высоту"));
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
            surfaceSimToggleButton_->setText(surfaceSimRunning_ ? QStringLiteral("Стоп")
                                                                : QStringLiteral("Старт"));
        }
        if (surfaceSimTimeLabel_) {
            const QString speedLabel =
                QStringLiteral("%1x").arg(surfaceSimSpeedMultiplier_, 0, 'g', 3);
            if (!surfaceSimRunning_ && surfaceSimState_.dayIndex == 0 &&
                surfaceSimState_.hourIndex == 0) {
                surfaceSimTimeLabel_->setText(
                    QStringLiteral("t = — (%1)").arg(speedLabel));
            } else {
                surfaceSimTimeLabel_->setText(
                    QStringLiteral("t = День %1, Час %2 (%3)")
                        .arg(surfaceSimState_.dayIndex + 1)
                        .arg(surfaceSimState_.hourIndex + 1)
                        .arg(speedLabel));
            }
        }
    }

    void updateSurfaceSimulationTimerInterval() {
        if (!surfaceSimTimer_) {
            return;
        }
        const int intervalMs = qMax(1, qRound(1000.0 / surfaceSimSpeedMultiplier_));
        surfaceSimTimer_->setInterval(intervalMs);
    }

    void resetSurfaceSimulation() {
        surfaceSimRunning_ = false;
        surfaceSimState_ = {};
        resetSurfaceTemperatureAggregation();
        resetSurfaceOrbitAnimation();
        if (surfaceSimTimer_) {
            surfaceSimTimer_->stop();
        }
        updateSurfaceSimulationUi();
    }

    void resetSurfaceOrbitAnimation() {
        if (!planetComboBox_) {
            surfaceOrbitAnimation_.reset(1.0, 0.0, 0.0, 0.0, kSurfaceOrbitSegmentsPerYear);
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
        surfaceOrbitAnimationInitialized_ = true;
    }

    void ensureSurfaceOrbitAnimationReady() {
        if (surfaceOrbitAnimationInitialized_) {
            return;
        }
        resetSurfaceOrbitAnimation();
    }

    void resetSurfaceTemperatureAggregation() {
        surfaceTemperatureAggregation_ = {};
        surfaceTemperatureAggregation_.historyDays = kSurfaceTemperatureHistoryDays;
        if (temperaturePlot_) {
            temperaturePlot_->clearSeries();
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
        surfacePointStatusDialog_->setPoint(*point);
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
        updatePlanetMassLabel();
        updatePlanetRadiusLabel();
        updatePlanetDerivedLabels();
        updateAtmospherePlanetParameters();
        updateAtmosphereComposition();
        updatePlanetOrbitLabels();
        syncMaterialWithPlanet();
        syncRotationModeWithPlanet();
        syncHeightSeedWithPlanet();
        syncFlatHeightWithPlanet();
        syncCloudAlbedoWithPlanet();
        syncManualGreenhouseOnTopWithPlanet();
        syncAdvancedRadiationWithPlanet();
        updatePlanetActions();
    }

    void clearPlanetPresets() {
        const QSignalBlocker blocker(planetComboBox_);
        planetComboBox_->clear();
        presetPlanetNames_.clear();
        planetSemiMajorAxisLabel_->setText(QStringLiteral("—"));
        planetDayLengthLabel_->setText(QStringLiteral("—"));
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
        if (flatHeightButton_) {
            const QSignalBlocker flatBlocker(flatHeightButton_);
            flatHeightButton_->setChecked(false);
            flatHeightButton_->setEnabled(false);
            updateFlatHeightButtonText(false);
        }
        if (cloudAlbedoSpinBox_) {
            const QSignalBlocker cloudBlocker(cloudAlbedoSpinBox_);
            cloudAlbedoSpinBox_->setValue(0.0);
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
        const QVariant compositionValue = planetComboBox_->currentData(kRoleAtmosphere);
        if (!compositionValue.isValid()) {
            atmosphereWidget_->setComposition(AtmosphereComposition{});
            return;
        }
        atmosphereWidget_->setComposition(compositionValue.value<AtmosphereComposition>());
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

    void addPlanetItem(const PlanetPreset &planet, bool isCustom) {
        planetComboBox_->addItem(formatPlanetName(planet), planet.semiMajorAxis);
        const int index = planetComboBox_->count() - 1;
        planetComboBox_->setItemData(index, planet.semiMajorAxis, kRoleSemiMajorAxis);
        planetComboBox_->setItemData(index, planet.dayLengthDays, kRoleDayLength);
        planetComboBox_->setItemData(index, planet.eccentricity, kRoleEccentricity);
        planetComboBox_->setItemData(index, planet.obliquityDegrees, kRoleObliquity);
        planetComboBox_->setItemData(index, planet.perihelionArgumentDegrees, kRolePerihelionArgument);
        planetComboBox_->setItemData(index, planet.massEarths, kRoleMassEarths);
        planetComboBox_->setItemData(index, planet.radiusKm, kRoleRadiusKm);
        const RotationMode rotationMode =
            planet.tidallyLocked ? RotationMode::TidalLocked : RotationMode::Normal;
        planetComboBox_->setItemData(index, static_cast<int>(rotationMode), kRoleRotationMode);
        planetComboBox_->setItemData(index, isCustom, kRoleIsCustom);
        planetComboBox_->setItemData(index, planet.name, kRolePlanetName);
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
        planetComboBox_->setItemData(index, planet.geothermalFluxWPerM2, kRoleGeothermalFlux);
        planetComboBox_->setItemData(index, static_cast<int>(planet.heightSourceType),
                                     kRoleHeightSourceType);
        planetComboBox_->setItemData(index, planet.heightmapPath, kRoleHeightmapPath);
        planetComboBox_->setItemData(index, planet.heightmapScaleKm, kRoleHeightmapScaleKm);
        planetComboBox_->setItemData(index, planet.heightSeed, kRoleHeightSeed);
        planetComboBox_->setItemData(index, planet.useContinentsHeight, kRoleUseContinentsHeight);
        planetComboBox_->setItemData(index, planet.hasSeaLevel, kRoleHasSeaLevel);
        planetComboBox_->setItemData(index, false, kRoleFlatHeight);
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
        auto *atmosphereInput = new AtmosphereWidget(&dialog, true);
        formLayout->addRow(QStringLiteral("Режим вращения:"), rotationModeInput);

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
                [&dialog, nameInput, axisInput, dayLengthInput, massInput, radiusInput,
                 eccentricityInput, obliquityInput, perihelionArgumentInput,
                 greenhouseOpacityInput, manualGreenhouseOnTopInput, cloudAlbedoInput,
                 heightSeedInput, materialInput,
                 rotationModeInput, atmosphereInput, this]() {
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
            const double geothermalFlux =
                subsurfaceGeothermalFluxSpinBox_ ? subsurfaceGeothermalFluxSpinBox_->value() : 0.0;

            const int existingIndex = findPlanetIndexByName(name);
            const QString materialId = materialInput->currentData().toString();
            const RotationMode rotationMode =
                static_cast<RotationMode>(rotationModeInput->currentData().toInt());
            const bool tidallyLocked = (rotationMode == RotationMode::TidalLocked);
            const quint32 heightSeed = static_cast<quint32>(heightSeedInput->value());
            const AtmosphereComposition composition = atmosphereInput->composition(false);
            const bool existingHasSeaLevel =
                existingIndex >= 0
                    ? planetComboBox_->itemData(existingIndex, kRoleHasSeaLevel).toBool()
                    : false;
            PlanetPreset preset{name, axis, dayLength, eccentricity, obliquity,
                                perihelionArgument, massEarths, radiusKm, materialId,
                                composition, greenhouseOpacity, manualGreenhouseOnTop,
                                cloudAlbedo, geothermalFlux, tidallyLocked};
            preset.heightSeed = heightSeed;
            preset.hasSeaLevel = existingHasSeaLevel;
            if (existingIndex >= 0) {
                if (!isCustomPlanetIndex(existingIndex)) {
                    showInputError(QStringLiteral("Нельзя заменить планету из пресета."));
                    return;
                }
                const QVariant atmosphereValue = QVariant::fromValue(composition);
                planetComboBox_->setItemText(existingIndex, formatPlanetName(preset));
                planetComboBox_->setItemData(existingIndex, axis, kRoleSemiMajorAxis);
                planetComboBox_->setItemData(existingIndex, dayLength, kRoleDayLength);
                planetComboBox_->setItemData(existingIndex, eccentricity, kRoleEccentricity);
                planetComboBox_->setItemData(existingIndex, obliquity, kRoleObliquity);
                planetComboBox_->setItemData(existingIndex, perihelionArgument,
                                             kRolePerihelionArgument);
                planetComboBox_->setItemData(existingIndex, massEarths, kRoleMassEarths);
                planetComboBox_->setItemData(existingIndex, radiusKm, kRoleRadiusKm);
                planetComboBox_->setItemData(existingIndex, static_cast<int>(rotationMode),
                                             kRoleRotationMode);
                planetComboBox_->setItemData(existingIndex, true, kRoleIsCustom);
                planetComboBox_->setItemData(existingIndex, name, kRolePlanetName);
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
                                             preset.geothermalFluxWPerM2,
                                             kRoleGeothermalFlux);
                planetComboBox_->setItemData(existingIndex,
                                             static_cast<int>(preset.heightSourceType),
                                             kRoleHeightSourceType);
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
                planetComboBox_->setCurrentIndex(existingIndex);
            } else {
                addPlanetItem(preset, true);
                planetComboBox_->setCurrentIndex(planetComboBox_->count() - 1);
            }

            updatePlanetSemiMajorAxisLabel();
            updatePlanetDayLengthLabel();
            updatePlanetMassLabel();
            updatePlanetRadiusLabel();
            updatePlanetDerivedLabels();
            updateAtmosphereComposition();
            updatePlanetOrbitLabels();
            syncMaterialWithPlanet();
            syncRotationModeWithPlanet();
            updatePlanetActions();
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

        const RotationMode rotationMode =
            static_cast<RotationMode>(planetComboBox_->itemData(index, kRoleRotationMode).toInt());
        const int modeIndex = rotationModeComboBox_->findData(static_cast<int>(rotationMode));
        if (modeIndex >= 0) {
            const QSignalBlocker blocker(rotationModeComboBox_);
            rotationModeComboBox_->setCurrentIndex(modeIndex);
        }
        updateRotationModeIllustration();
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

    void syncCloudAlbedoWithPlanet() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !cloudAlbedoSpinBox_) {
            return;
        }

        const double cloudAlbedo = planetComboBox_->itemData(index, kRoleCloudAlbedo).toDouble();
        const QSignalBlocker blocker(cloudAlbedoSpinBox_);
        cloudAlbedoSpinBox_->setValue(cloudAlbedo);
    }

    void syncGeothermalFluxWithPlanet() {
        if (!subsurfaceGeothermalFluxSpinBox_) {
            return;
        }
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }

        const QString planetName = planetComboBox_->itemData(index, kRolePlanetName).toString();
        if (planetName != QStringLiteral("Луна")) {
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
    }

    void syncPlanetHeightSeedWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !heightSeedSpinBox_) {
            return;
        }
        planetComboBox_->setItemData(index, heightSeedSpinBox_->value(), kRoleHeightSeed);
    }

    void syncPlanetCloudAlbedoWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0 || !cloudAlbedoSpinBox_) {
            return;
        }
        planetComboBox_->setItemData(index, cloudAlbedoSpinBox_->value(), kRoleCloudAlbedo);
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

    void clearTemperatureSegments() {
        segmentSelectorWidget_->setSegments({});
        segmentSelectorWidget_->setEnabled(false);
        temperaturePlot_->clearSeries();
    }

    void rebuildSurfaceGrid() {
        const double radiusKm = planetComboBox_->currentData(kRoleRadiusKm).toDouble();
        surfaceGrid_.setRadiusKm(radiusKm);
        const HeightSourceType heightSource =
            static_cast<HeightSourceType>(planetComboBox_->currentData(kRoleHeightSourceType)
                                              .toInt());
        const QString heightmapPath =
            planetComboBox_->currentData(kRoleHeightmapPath).toString();
        const double heightmapScaleKm =
            planetComboBox_->currentData(kRoleHeightmapScaleKm).toDouble();
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

        const int latitudePointCount = latitudePoints();
        const int pointsPerLatitude = 6;
        const int targetPointCount = qMax(1, latitudePointCount * pointsPerLatitude);
        // Число ячеек в геодезической сетке равно 20 * 4^n, подбираем n под желаемое
        // количество точек, чтобы сохранить приблизительную плотность сетки.
        const double ratio = qMax(1.0, static_cast<double>(targetPointCount) / 20.0);
        const int subdivisionLevel = qMax(0, static_cast<int>(qRound(qLn(ratio) / qLn(4.0))));
        surfaceGrid_.generateIcosahedronGrid(subdivisionLevel);

        if (useFlatHeight) {
            for (auto &point : surfaceGrid_.points()) {
                point.heightKm = 0.0;
            }
        }
    }

    struct SurfacePointStateDefaults {
        double albedo = 0.0;
        double greenhouseOpacity = 0.0;
        double minTemperatureKelvin = 3.0;
        SurfaceMaterial material;
        SubsurfaceModelSettings subsurfaceSettings;
    };

    std::optional<SurfacePointStateDefaults> buildSurfacePointStateDefaults() const {
        const auto material = currentMaterial();
        if (!material) {
            return std::nullopt;
        }

        const double greenhouseOpacity =
            planetComboBox_->currentData(kRoleGreenhouseOpacity).toDouble();
        const double albedo = qBound(0.0, material->albedo, 1.0);
        // Не задаем высокий нижний порог: карта поверхности должна стартовать от физического минимума,
        // чтобы при слабой инсоляции температура могла быть значительно ниже 200 K.
        const double minTemperatureKelvin = 3.0;

        SurfacePointStateDefaults defaults;
        defaults.albedo = albedo;
        defaults.greenhouseOpacity = greenhouseOpacity;
        defaults.minTemperatureKelvin = minTemperatureKelvin;
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
        SubsurfaceModelSettings settings;
        if (subsurfaceLayersSpinBox_) {
            settings.layerCount = subsurfaceLayersSpinBox_->value();
        }
        if (subsurfaceTopThicknessSpinBox_) {
            settings.topLayerThicknessMeters = subsurfaceTopThicknessSpinBox_->value();
        }
        if (subsurfaceDepthSpinBox_) {
            settings.bottomDepthMeters = subsurfaceDepthSpinBox_->value();
        }
        if (subsurfaceBoundaryComboBox_) {
            settings.bottomBoundary = static_cast<SubsurfaceBottomBoundaryCondition>(
                subsurfaceBoundaryComboBox_->currentData().toInt());
        }
        if (subsurfaceGeothermalFluxSpinBox_) {
            settings.geothermalFluxWPerM2 = subsurfaceGeothermalFluxSpinBox_->value();
        }
        return settings;
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

    double resolveAirTemperatureKelvin(const SurfacePointState &surfaceState,
                                       double pressureAtm,
                                       double gravity,
                                       double initialAirTemperature,
                                       double blendedInsolation,
                                       double cloudShortwaveTransmission,
                                       RadiationModelType radiationModelType,
                                       const RadiationModel &radiationModel,
                                       double timeStepSeconds) const {
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
                                            radiationModel,
                                            timeStepSeconds);
    }

    double estimateAirTemperatureKelvin(const SurfacePointState &surfaceState,
                                        double pressureAtm,
                                        double gravity,
                                        double initialAirTemperature,
                                        double blendedInsolation,
                                        double cloudShortwaveTransmission,
                                        const RadiationModel &radiationModel,
                                        double timeStepSeconds) const {
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
        const double couplingScale = qBound(0.0, pressureAtm, 1.0);
        const double heatTransfer = kDefaultHeatTransferWPerM2K * couplingScale;
        if (timeStepSeconds <= 0.0) {
            return airState.airTemperatureKelvin();
        }

        // Радиационный нагрев воздуха:
        // Q_sw_air = S_blend * T_cloud * (1 - T_atm), где T_atm = incomingTransmission().
        // Q_lw_air = F_surf * (1 - T_lw), где T_lw = outgoingTransmission().
        const double emittedFlux = surfaceState.emittedFlux();
        const double shortwaveAbsorbedByAir =
            blendedInsolation * cloudShortwaveTransmission *
            (1.0 - radiationModel.incomingTransmission());
        const double longwaveAbsorbedByAir =
            emittedFlux * (1.0 - radiationModel.outgoingTransmission());
        const double airRadiativeHeatingFlux = shortwaveAbsorbedByAir + longwaveAbsorbedByAir;
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

    void updateSurfaceWindField(const AtmosphereComposition &atmosphere,
                                double atmospherePressureAtm,
                                double dayLengthDays,
                                double surfaceGravity) {
        Q_UNUSED(atmosphere)
        Q_UNUSED(atmospherePressureAtm)
        Q_UNUSED(surfaceGravity)
        if (surfaceGrid_.points().isEmpty()) {
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

    void updateSurfaceGridTemperatures() {
        resetSurfaceSimulation();
        rebuildSurfaceGrid();
        updateSurfaceHeightLegendFromGrid();
        AtmosphereComposition atmosphere;
        const QVariant atmosphereValue = planetComboBox_->currentData(kRoleAtmosphere);
        if (atmosphereValue.isValid()) {
            atmosphere = atmosphereValue.value<AtmosphereComposition>();
        }
        const double massEarths = planetComboBox_->currentData(kRoleMassEarths).toDouble();
        const double radiusKm = planetComboBox_->currentData(kRoleRadiusKm).toDouble();
        const double cloudAlbedo =
            qBound(0.0, planetComboBox_->currentData(kRoleCloudAlbedo).toDouble(), 1.0);
        const QString baseMaterialId =
            resolveSurfaceMaterialIdForPoint(atmosphere, massEarths, radiusKm);
        for (auto &point : surfaceGrid_.points()) {
            if (surfaceMaxHeightKm_ > 0.0 &&
                point.heightKm >= 0.75 * surfaceMaxHeightKm_) {
                // Высотный порог для скальных пород: доля от максимума дает масштабируемый
                // критерий для разных планетных рельефов.
                point.materialId = QStringLiteral("rocky");
            } else {
                point.materialId = baseMaterialId;
            }
        }
        if (surfaceGrid_.points().isEmpty()) {
            applySurfaceGridToViews();
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceAirTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceWindLegend(false, 0.0, 0.0);
            updateSurfacePressureLegend(false, 0.0, 0.0);
            return;
        }

        const auto stateDefaults = buildSurfacePointStateDefaults();
        if (!stateDefaults) {
            applySurfaceGridToViews();
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceAirTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceWindLegend(false, 0.0, 0.0);
            updateSurfacePressureLegend(false, 0.0, 0.0);
            return;
        }

        QHash<QString, SurfaceMaterial> materialsById;
        const auto materials = surfaceMaterials();
        materialsById.reserve(materials.size());
        for (const auto &material : materials) {
            materialsById.insert(material.id, material);
        }
        const auto materialForPoint = [&materialsById, &stateDefaults](const QString &materialId) {
            const auto it = materialsById.constFind(materialId);
            return it != materialsById.cend() ? *it : stateDefaults->material;
        };

        const double dayLengthDays = planetComboBox_->currentData(kRoleDayLength).toDouble();
        double atmospherePressureAtm = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            atmospherePressureAtm = atmosphere.totalPressureAtm(massEarths, radiusKm);
        }
        const bool useAtmosphericModel = atmosphere.totalMassGigatons() > 0.0;
        const double localSeaLevelPressureAtm =
            calculateCellPressureAtmFromKg(atmosphere.totalMassKg(),
                                           massEarths,
                                           radiusKm,
                                           surfaceGrid_.pointAreaKm2());
        double surfaceGravity = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            const double radiusMeters = radiusKm * 1000.0;
            const double planetMassKg = massEarths * kEarthMassKg;
            surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);
        }
        const double gravity = (surfaceGravity > 0.0) ? surfaceGravity : 9.80665;
        const double manualGreenhouseOpacity = stateDefaults->greenhouseOpacity;
        const bool manualGreenhouseOnTop =
            planetComboBox_->currentData(kRoleManualGreenhouseOnTopOfAtmosphere).toBool();
        const auto radiationModelType = static_cast<RadiationModelType>(
            planetComboBox_->currentData(kRoleAdvancedRadiationModel).toInt());
        const RotationMode rotationMode =
            static_cast<RotationMode>(planetComboBox_->currentData(kRoleRotationMode).toInt());
        ensureSurfaceOrbitAnimationReady();
        const double declinationDegrees = surfaceOrbitAnimation_.declinationDegrees();
        double segmentSolarConstant = hasSolarConstant_ ? lastSolarConstant_ : 0.0;
        if (hasSolarConstant_ && lastSolarConstantDistanceAU_ > 0.0) {
            const double distanceAU = surfaceOrbitAnimation_.distanceAU();
            segmentSolarConstant =
                lastSolarConstant_ *
                std::pow(lastSolarConstantDistanceAU_ / distanceAU, 2.0);
        }

        SurfaceTileTemperatureDefaults tileDefaults;
        tileDefaults.minTemperatureKelvin = stateDefaults->minTemperatureKelvin;
        tileDefaults.greenhouseOpacity = stateDefaults->greenhouseOpacity;
        tileDefaults.subsurfaceSettings = stateDefaults->subsurfaceSettings;

        SurfaceTileTemperatureSettings tileSettings;
        tileSettings.segmentSolarConstant = segmentSolarConstant;
        tileSettings.dayLengthDays = dayLengthDays;
        tileSettings.rotationMode = rotationMode;
        tileSettings.declinationDegrees = declinationDegrees;
        tileSettings.atmospherePressureAtm = atmospherePressureAtm;
        tileSettings.cloudAlbedo = cloudAlbedo;
        tileSettings.currentHourIndex = surfaceSimState_.hourIndex;
        tileSettings.hasSolarConstant = hasSolarConstant_;

        SurfaceTileTemperatureCalculator tileCalculator;
        SurfaceTileTemperatureResult tileResult =
            tileCalculator.initializeSurface(surfaceGrid_,
                                             tileSettings,
                                             tileDefaults,
                                             materialsById,
                                             stateDefaults->material);

        bool hasTemperatureRange = tileResult.hasTemperatureRange;
        double minTemperature = tileResult.minSurfaceTemperatureK;
        double maxTemperature = tileResult.maxSurfaceTemperatureK;
        double minAirTemperature = std::numeric_limits<double>::max();
        double maxAirTemperature = std::numeric_limits<double>::lowest();
        QVector<double> blendedInsolations = tileResult.blendedInsolations;
        QVector<double> baselineAirTemperatures = tileResult.baselineAirTemperatures;
        const double timeStepSeconds = 3600.0;
        // Коэффициент прохождения коротковолнового излучения через облака.
        double cloudShortwaveTransmission = 1.0 - cloudAlbedo;
        if (cloudAlbedo > 0.7) {
            // Для плотных сернокислотных облаков дополнительно ослабляем поток к поверхности.
            cloudShortwaveTransmission *= 0.2;
        }
        cloudShortwaveTransmission = qBound(0.0, cloudShortwaveTransmission, 1.0);
        if (blendedInsolations.size() != surfaceGrid_.points().size()) {
            blendedInsolations.clear();
            blendedInsolations.resize(surfaceGrid_.points().size());
        }
        if (baselineAirTemperatures.size() != surfaceGrid_.points().size()) {
            baselineAirTemperatures.clear();
            baselineAirTemperatures.resize(surfaceGrid_.points().size());
        }

        double minPressureAtm = std::numeric_limits<double>::max();
        double maxPressureAtm = std::numeric_limits<double>::lowest();
        for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
            auto &point = surfaceGrid_.points()[i];
            // Инициализируем поверхностное давление (на уровне рельефа), чтобы дальше
            // переносить его ветром без пересчёта из всей атмосферы каждый тик.
            // Это нужно и в fallback-ветке, чтобы карта давления/ветра не оставалась с нулями.
            // Масса газового столба в ячейке пропорциональна её площади, поэтому P = (m_cell * g) / area_cell.
            const double pressureAtm =
                AtmosphericPressureModel::pressureAtHeightAtm(localSeaLevelPressureAtm,
                                                              point.heightKm * 1000.0,
                                                              point.temperatureK,
                                                              atmosphere,
                                                              gravity);
            // Храним поверхностное давление в точке для последующего переноса по ветру.
            point.pressureAtm = qMax(0.0, pressureAtm);
            minPressureAtm = qMin(minPressureAtm, point.pressureAtm);
            maxPressureAtm = qMax(maxPressureAtm, point.pressureAtm);
            const double blendedInsolation =
                (i < blendedInsolations.size()) ? blendedInsolations.at(i) : 0.0;
            const SurfaceMaterial material = materialForPoint(point.materialId);
            const bool logDetails = shouldLogRadiationForPoint(i);
            const double localGreenhouseOpacity =
                computeLocalGreenhouseOpacity(atmosphere,
                                              material,
                                              point.pressureAtm,
                                              radiusKm,
                                              gravity,
                                              blendedInsolation,
                                              manualGreenhouseOpacity,
                                              useAtmosphericModel,
                                              radiationModelType,
                                              manualGreenhouseOnTop,
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
            // Используем локальную оценку ТОА-потока через blendedInsolation и планетарное
            // альбедо (облака перекрывают поверхность), чтобы связать оптическую толщину
            // с физически корректным источником излучения.
            const double materialAlbedo = qBound(0.0, material.albedo, 1.0);
            const double planetaryAlbedo = qMax(materialAlbedo, cloudAlbedo);
            const double effectiveFlux =
                blendedInsolation * qMax(0.0, 1.0 - planetaryAlbedo);
            const double effectiveTemperatureKelvin =
                std::pow(qMax(0.0, effectiveFlux) / kStefanBoltzmannConstant, 0.25);
            if (logDetails) {
                qCInfo(solarRadiationLog) << "Radiation inputs (init)"
                                          << "index=" << i
                                          << "blendedInsolation=" << blendedInsolation
                                          << "planetaryAlbedo=" << planetaryAlbedo
                                          << "effectiveFlux=" << effectiveFlux
                                          << "effectiveTemperatureKelvin=" << effectiveTemperatureKelvin;
            }
            const auto radiationModel =
                makeRadiationModel(atmosphere,
                                   point.pressureAtm,
                                   point.state.temperatureKelvin(),
                                   effectiveTemperatureKelvin,
                                   gravity,
                                   radiationModelType);
            point.airTemperatureK =
                resolveAirTemperatureKelvin(point.state,
                                            point.pressureAtm,
                                            gravity,
                                            initialAirTemperature,
                                            blendedInsolation,
                                            cloudShortwaveTransmission,
                                            radiationModelType,
                                            *radiationModel,
                                            timeStepSeconds);
            if (logDetails) {
                qCInfo(solarRadiationLog) << "Resolved air temperature (init)"
                                          << "index=" << i
                                          << "airTemperatureK=" << point.airTemperatureK;
            }
            minAirTemperature = qMin(minAirTemperature, point.airTemperatureK);
            maxAirTemperature = qMax(maxAirTemperature, point.airTemperatureK);
        }
        updateSurfaceWindField(atmosphere, atmospherePressureAtm, dayLengthDays, surfaceGravity);

        applySurfaceGridToViews();
        if (hasTemperatureRange && minTemperature <= maxTemperature) {
            applySurfaceTemperatureRangeToViews(minTemperature, maxTemperature);
            updateSurfaceTemperatureLegend(true, minTemperature, maxTemperature);
        } else {
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
        }
        if (hasTemperatureRange && minAirTemperature <= maxAirTemperature) {
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
    }

    void toggleSurfaceSimulation() {
        if (surfaceSimRunning_) {
            surfaceSimRunning_ = false;
            if (surfaceSimTimer_) {
                surfaceSimTimer_->stop();
            }
            updateSurfaceSimulationUi();
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
            static_cast<RotationMode>(planetComboBox_->currentData(kRoleRotationMode).toInt());
        AtmosphereComposition atmosphere;
        const QVariant atmosphereValue = planetComboBox_->currentData(kRoleAtmosphere);
        if (atmosphereValue.isValid()) {
            atmosphere = atmosphereValue.value<AtmosphereComposition>();
        }
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

        const int stepsPerDay = qMax(1, qRound(dayLengthDays * 24.0));

        ensureSurfaceOrbitAnimationReady();
        // Сезонная деклинация: δ = asin(sin(наклон оси) * sin(истинная долгота звезды)).
        const double declinationDegrees = surfaceOrbitAnimation_.declinationDegrees();

        // Инсоляция меняется с расстоянием как 1 / r^2 относительно опорной дистанции.
        const double distanceAU = surfaceOrbitAnimation_.distanceAU();
        const double segmentSolarConstant =
            lastSolarConstant_ *
            std::pow(lastSolarConstantDistanceAU_ / distanceAU, 2.0);
        const double transport =
            (atmospherePressureAtm > 50.0)
                ? 0.99
                : (atmospherePressureAtm > 0.001
                       ? qMin(1.0, 0.15 * std::log(atmospherePressureAtm * 100.0 + 1.0))
                       : 0.0);
        const double rotBlock =
            (dayLengthDays < 2.0 && atmospherePressureAtm < 10.0) ? 0.65 : 1.0;
        const double meridionalTransport = transport * rotBlock;
        // Глобальный средний поток перед альбедо, как в SurfaceTemperatureCalculator.
        const double globalAverageInsolation = segmentSolarConstant / 4.0;

        // Один тик = 1 час планетарных суток, ускорение реализовано уменьшением интервала таймера.
        const double timeStepSeconds = 3600.0;
        // Коэффициент прохождения коротковолнового излучения через облака.
        double cloudShortwaveTransmission = 1.0 - cloudAlbedo;
        if (cloudAlbedo > 0.7) {
            // Для плотных сернокислотных облаков дополнительно ослабляем поток к поверхности.
            cloudShortwaveTransmission *= 0.2;
        }
        cloudShortwaveTransmission = qBound(0.0, cloudShortwaveTransmission, 1.0);
        const double phase =
            2.0 * M_PI *
            (static_cast<double>(surfaceSimState_.hourIndex) + 0.5) /
            static_cast<double>(stepsPerDay);
        const double hourAngle =
            (rotationMode == RotationMode::TidalLocked) ? 0.0 : (phase - M_PI);
        // Субзвёздная долгота задает меридиан с нулевым часовым углом.
        const double substellarLongitudeRadians =
            (rotationMode == RotationMode::TidalLocked) ? 0.0 : -hourAngle;
        const double declinationRadians = qDegreesToRadians(declinationDegrees);
        const double sinDeclination = std::sin(declinationRadians);
        const double cosDeclination = std::cos(declinationRadians);
        // Фиксируем меридиан для агрегации: либо пользовательскую точку,
        // либо ближайший к подсолнечному меридиану, чтобы брать единый срез по широтам.
        const double targetLongitudeRadians =
            surfaceTemperatureAggregation_.hasTargetLongitude
                ? surfaceTemperatureAggregation_.targetLongitudeRadians
                : selectSurfaceAggregationLongitude(substellarLongitudeRadians);
        ensureSurfaceTemperatureAggregationTargets(targetLongitudeRadians);

        QVector<double> blendedInsolations;
        blendedInsolations.reserve(surfaceGrid_.points().size());
        for (auto &point : surfaceGrid_.points()) {
            const double localHourAngle = point.longitudeRadians - substellarLongitudeRadians;
            const double cosZenith =
                point.sinLatitude * sinDeclination +
                point.cosLatitude * cosDeclination * std::cos(localHourAngle);
            // S_inst = S0 * cos(zenith) при освещении, иначе 0.
            const double localInsolation =
                segmentSolarConstant * qMax(0.0, cosZenith);
            const double blendedInsolation =
                localInsolation * (1.0 - meridionalTransport) +
                globalAverageInsolation * meridionalTransport;
            blendedInsolations.push_back(blendedInsolation);

            const double absorbedFlux = point.state.absorbedFlux(blendedInsolation);
            const double emittedFlux = point.state.emittedFlux();
            // Применяем шаговый радиационный баланс для состояния точки.
            point.state.updateTemperature(absorbedFlux, emittedFlux, timeStepSeconds);
            // Обновляем поверхностную температуру до переноса в соседние точки.
            point.temperatureK = point.state.temperatureKelvin();
        }

        updateSurfaceWindField(atmosphere, atmospherePressureAtm, dayLengthDays, surfaceGravity);

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
        const double gravity = (surfaceGravity > 0.0) ? surfaceGravity : 9.80665;
        constexpr double kPressureRelaxFactor = 0.1;
        double minPressureAtm = std::numeric_limits<double>::max();
        double maxPressureAtm = std::numeric_limits<double>::lowest();
        for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
            auto &point = surfaceGrid_.points()[i];
            const double advectedAtm = advectedPressures.at(i);
            // Барометрическая релаксация возвращает поле давления к
            // физически согласованному профилю P(z) при текущей температуре.
            const double barometricAtm =
                AtmosphericPressureModel::pressureAtHeightAtm(localSeaLevelPressureAtm,
                                                              point.heightKm * 1000.0,
                                                              point.temperatureK,
                                                              atmosphere,
                                                              gravity);
            const double relaxedAtm =
                advectedAtm + (barometricAtm - advectedAtm) * kPressureRelaxFactor;
            // Фиксируем перенесённое поверхностное давление для UI и динамики.
            point.pressureAtm = qMax(0.0, relaxedAtm);
            minPressureAtm = qMin(minPressureAtm, point.pressureAtm);
            maxPressureAtm = qMax(maxPressureAtm, point.pressureAtm);
            const double blendedInsolation =
                (i < blendedInsolations.size()) ? blendedInsolations.at(i) : 0.0;
            const SurfaceMaterial material = materialForPoint(point.materialId);
            const bool logDetails = shouldLogRadiationForPoint(i);
            const double localGreenhouseOpacity =
                computeLocalGreenhouseOpacity(atmosphere,
                                              material,
                                              point.pressureAtm,
                                              radiusKm,
                                              gravity,
                                              blendedInsolation,
                                              manualGreenhouseOpacity,
                                              useAtmosphericModel,
                                              radiationModelType,
                                              manualGreenhouseOnTop,
                                              logDetails);
            point.state.setGreenhouseOpacity(localGreenhouseOpacity);
        }

        SurfaceAdvectionModel advectionModel;
        QVector<double> advectedTemperatures =
            advectionModel.advectTemperature(surfaceGrid_,
                                             temperatures,
                                             windEast,
                                             windNorth,
                                             timeStepSeconds,
                                             1,
                                             1.0);
        if (advectedTemperatures.size() != surfaceGrid_.points().size()) {
            advectedTemperatures = temperatures;
        }

        double minTemperature = std::numeric_limits<double>::max();
        double maxTemperature = std::numeric_limits<double>::lowest();
        double minAirTemperature = std::numeric_limits<double>::max();
        double maxAirTemperature = std::numeric_limits<double>::lowest();
        for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
            auto &point = surfaceGrid_.points()[i];
            point.state.setTemperatureKelvin(advectedTemperatures.at(i));
            // Температура после переноса хранится как поверхностная величина.
            point.temperatureK = point.state.temperatureKelvin();
            minTemperature = qMin(minTemperature, point.temperatureK);
            maxTemperature = qMax(maxTemperature, point.temperatureK);
            const double initialAirTemperature =
                (point.airTemperatureK > 0.0) ? point.airTemperatureK : point.temperatureK;
            const double blendedInsolation =
                (i < blendedInsolations.size()) ? blendedInsolations.at(i) : 0.0;
            const SurfaceMaterial material = materialForPoint(point.materialId);
            const bool logDetails = shouldLogRadiationForPoint(i);
            const double materialAlbedo = qBound(0.0, material.albedo, 1.0);
            const double planetaryAlbedo = qMax(materialAlbedo, cloudAlbedo);
            const double effectiveFlux =
                blendedInsolation * qMax(0.0, 1.0 - planetaryAlbedo);
            const double effectiveTemperatureKelvin =
                std::pow(qMax(0.0, effectiveFlux) / kStefanBoltzmannConstant, 0.25);
            if (logDetails) {
                qCInfo(solarRadiationLog) << "Radiation inputs (tick)"
                                          << "index=" << i
                                          << "blendedInsolation=" << blendedInsolation
                                          << "planetaryAlbedo=" << planetaryAlbedo
                                          << "effectiveFlux=" << effectiveFlux
                                          << "effectiveTemperatureKelvin=" << effectiveTemperatureKelvin;
            }
            const auto radiationModel =
                makeRadiationModel(atmosphere,
                                   point.pressureAtm,
                                   point.state.temperatureKelvin(),
                                   effectiveTemperatureKelvin,
                                   gravity,
                                   radiationModelType);
            point.airTemperatureK =
                resolveAirTemperatureKelvin(point.state,
                                            point.pressureAtm,
                                            gravity,
                                            initialAirTemperature,
                                            blendedInsolation,
                                            cloudShortwaveTransmission,
                                            radiationModelType,
                                            *radiationModel,
                                            timeStepSeconds);
            if (logDetails) {
                qCInfo(solarRadiationLog) << "Resolved air temperature (tick)"
                                          << "index=" << i
                                          << "airTemperatureK=" << point.airTemperatureK;
            }
            minAirTemperature = qMin(minAirTemperature, point.airTemperatureK);
            maxAirTemperature = qMax(maxAirTemperature, point.airTemperatureK);
        }

        // Обновляем карту после каждого тика таймера, чтобы сразу отражать новую температуру.
        // Также обновляем виджеты, чтобы в них попали обновлённые поверхностные
        // температура и давление после шага переноса.
        applySurfaceGridToViews();
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

        const bool publishDaily = surfaceSimState_.hourIndex + 1 >= stepsPerDay;
        updateSurfaceTemperatureAggregation(publishDaily, rotationMode, useAtmosphericModel);

        ++surfaceSimState_.hourIndex;
        if (surfaceSimState_.hourIndex >= stepsPerDay) {
            surfaceSimState_.hourIndex = 0;
            ++surfaceSimState_.dayIndex;
            surfaceOrbitAnimation_.advanceSegment();
        }
        updateSurfaceSimulationUi();
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
        if (temperatureElapsedLabel_) {
            temperatureElapsedLabel_->setText(QStringLiteral("Прошло: 00:00"));
        }
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
            if (temperatureElapsedLabel_) {
                temperatureElapsedLabel_->setText(
                    QStringLiteral("Прошло: %1:%2%3")
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

    SolarCalculatorWidget widget(precision);
    widget.setWindowTitle(QStringLiteral("Калькулятор солнечной постоянной"));
    widget.show();

    return app.exec();
}
