#include "orbit_segment_calculator.h"
#include "mode_illustration_widget.h"
#include "orbit_animation_model.h"
#include "planet_presets.h"
#include "segment_selector_widget.h"
#include "solar_calculator.h"
#include "solar_display.h"
#include "atmosphere_widget.h"
#include "surface_temperature_calculator.h"
#include "surface_temperature_plot.h"
#include "surface_map_widget.h"
#include "surface_globe_widget.h"
#include "surface_point_status_dialog.h"
#include "surface_temperature_scale_widget.h"
#include "surface_height_scale_widget.h"
#include "surface_wind_scale_widget.h"
#include "surface_map_mode.h"
#include "surface_advection_model.h"
#include "surface_pressure_transport_model.h"
#include "wind_field_model.h"
#include "planet_surface_grid.h"
#include "subsurface_temperature_solver.h"
#include "atmospheric_pressure_model.h"

#include <QtCore/QCommandLineOption>
#include <QtCore/QCommandLineParser>
#include <QtCore/QElapsedTimer>
#include <QtCore/QLocale>
#include <QtCore/QPointer>
#include <QtCore/QSignalBlocker>
#include <QtCore/QSysInfo>
#include <QtCore/QThread>
#include <QtCore/QThreadPool>
#include <QtCore/QTimer>
// #include <QtCore/QOverload>
#include <QtGlobal>
#include <QtCore/QtMath>
#include <QtCore/QSet>
#include <QtCore/QFutureWatcher>
#include <QtCore/QHash>
#include <QtCore/QHashFunctions>
#include <QtCore/QStringList>
#include <QtEndian>
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
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QStyle>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <cstring>

namespace {
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
constexpr int kRoleHeightSourceType = Qt::UserRole + 13;
constexpr int kRoleHeightmapPath = Qt::UserRole + 14;
constexpr int kRoleHeightmapScaleKm = Qt::UserRole + 15;
constexpr int kRoleHeightSeed = Qt::UserRole + 16;
constexpr int kRoleUseContinentsHeight = Qt::UserRole + 17;
constexpr double kKelvinOffset = 273.15;
constexpr double kEarthRadiusKm = 6371.0;
constexpr double kEarthMassKg = 5.9722e24;
constexpr double kGravitationalConstant = 6.67430e-11;
constexpr int kSurfaceOrbitSegmentsPerYear = 360;

struct TemperatureCacheKey {
    double solarConstant = 0.0;
    QString materialId;
    QString atmosphereSignature;
    double atmospherePressureAtm = 0.0;
    double surfaceGravity = 0.0;
    double greenhouseOpacity = 0.0;
    double dayLength = 0.0;
    double referenceDistanceAU = 0.0;
    double semiMajorAxis = 0.0;
    double eccentricity = 0.0;
    double obliquity = 0.0;
    double perihelionArgument = 0.0;
    double planetRadiusKm = 0.0;
    HeightSourceType heightSourceType = HeightSourceType::Procedural;
    QString heightmapPath;
    double heightmapScaleKm = 0.0;
    quint32 heightSeed = 0;
    bool useContinentsHeight = false;
    int subsurfaceLayers = 0;
    double subsurfaceTopThicknessMeters = 0.0;
    double subsurfaceDepthMeters = 0.0;
    SubsurfaceBottomBoundaryCondition subsurfaceBoundary =
        SubsurfaceBottomBoundaryCondition::Insulating;
    int latitudePoints = 0;
    int segmentCount = 0;
    RotationMode rotationMode = RotationMode::Normal;

    bool operator==(const TemperatureCacheKey &other) const {
        return solarConstant == other.solarConstant &&
               materialId == other.materialId &&
               atmosphereSignature == other.atmosphereSignature &&
               atmospherePressureAtm == other.atmospherePressureAtm &&
               surfaceGravity == other.surfaceGravity &&
               greenhouseOpacity == other.greenhouseOpacity &&
               dayLength == other.dayLength &&
               referenceDistanceAU == other.referenceDistanceAU &&
               semiMajorAxis == other.semiMajorAxis &&
               eccentricity == other.eccentricity &&
               obliquity == other.obliquity &&
               perihelionArgument == other.perihelionArgument &&
               planetRadiusKm == other.planetRadiusKm &&
               heightSourceType == other.heightSourceType &&
               heightmapPath == other.heightmapPath &&
               heightmapScaleKm == other.heightmapScaleKm &&
               heightSeed == other.heightSeed &&
               useContinentsHeight == other.useContinentsHeight &&
               subsurfaceLayers == other.subsurfaceLayers &&
               subsurfaceTopThicknessMeters == other.subsurfaceTopThicknessMeters &&
               subsurfaceDepthMeters == other.subsurfaceDepthMeters &&
               subsurfaceBoundary == other.subsurfaceBoundary &&
               latitudePoints == other.latitudePoints &&
               segmentCount == other.segmentCount &&
               rotationMode == other.rotationMode;
    }
};

uint qHash(const TemperatureCacheKey &key, uint seed = 0) {
    using ::qHash;
    const auto hashDoubleBits = [](double value) -> quint64 {
        // Хэшируем битовое представление (в big-endian), чтобы не зависеть от qHash(double).
        quint64 bits = 0;
        static_assert(sizeof(bits) == sizeof(value), "Unexpected double size.");
        std::memcpy(&bits, &value, sizeof(bits));
        if (QSysInfo::ByteOrder == QSysInfo::LittleEndian) {
            bits = qbswap(bits);
        }
        return bits;
    };

    seed = qHash(hashDoubleBits(key.solarConstant), seed);
    seed = qHash(key.materialId, seed);
    seed = qHash(key.atmosphereSignature, seed);
    seed = qHash(hashDoubleBits(key.atmospherePressureAtm), seed);
    seed = qHash(hashDoubleBits(key.surfaceGravity), seed);
    seed = qHash(hashDoubleBits(key.greenhouseOpacity), seed);
    seed = qHash(hashDoubleBits(key.dayLength), seed);
    seed = qHash(hashDoubleBits(key.referenceDistanceAU), seed);
    seed = qHash(hashDoubleBits(key.semiMajorAxis), seed);
    seed = qHash(hashDoubleBits(key.eccentricity), seed);
    seed = qHash(hashDoubleBits(key.obliquity), seed);
    seed = qHash(hashDoubleBits(key.perihelionArgument), seed);
    seed = qHash(hashDoubleBits(key.planetRadiusKm), seed);
    seed = qHash(static_cast<int>(key.heightSourceType), seed);
    seed = qHash(key.heightmapPath, seed);
    seed = qHash(hashDoubleBits(key.heightmapScaleKm), seed);
    seed = qHash(static_cast<quint32>(key.heightSeed), seed);
    seed = qHash(key.useContinentsHeight, seed);
    seed = qHash(key.subsurfaceLayers, seed);
    seed = qHash(hashDoubleBits(key.subsurfaceTopThicknessMeters), seed);
    seed = qHash(hashDoubleBits(key.subsurfaceDepthMeters), seed);
    seed = qHash(static_cast<int>(key.subsurfaceBoundary), seed);
    seed = qHash(key.latitudePoints, seed);
    seed = qHash(key.segmentCount, seed);
    seed = qHash(static_cast<int>(key.rotationMode), seed);
    return seed;
}

static QString atmosphereSignature(const AtmosphereComposition &composition) {
    auto fractions = composition.fractions();
    std::sort(fractions.begin(), fractions.end(),
              [](const GasFraction &left, const GasFraction &right) {
                  return left.id < right.id;
              });
    QStringList parts;
    for (const auto &fraction : fractions) {
        if (fraction.massGigatons <= 0.0) {
            continue;
        }
        parts << QStringLiteral("%1:%2").arg(fraction.id,
                                             QString::number(fraction.massGigatons, 'g', 10));
    }
    return parts.isEmpty() ? QStringLiteral("none") : parts.join('|');
}

struct TemperatureCacheEntry {
    QVector<OrbitSegment> orbitSegments;
    QVector<QVector<TemperatureRangePoint>> temperatureSegments;
};

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
        modeIllustrationWidget_ = new ModeIllustrationWidget(this);
        modeIllustrationWidget_->setRotationMode(
            static_cast<RotationMode>(rotationModeComboBox_->currentData().toInt()));
        latitudeStepFastRadio_ = new QRadioButton(QStringLiteral("Через 1° (быстрые)"), this);
        latitudeStepSlowRadio_ = new QRadioButton(QStringLiteral("Через 10° (медленные)"), this);
        latitudeStepGroup_ = new QButtonGroup(this);
        latitudeStepGroup_->addButton(latitudeStepFastRadio_);
        latitudeStepGroup_->addButton(latitudeStepSlowRadio_);
        latitudeStepFastRadio_->setChecked(true);
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
        auto *latitudeStepLayout = new QHBoxLayout();
        latitudeStepLayout->addWidget(latitudeStepFastRadio_);
        latitudeStepLayout->addWidget(latitudeStepSlowRadio_);
        auto *latitudeStepWidget = new QWidget(this);
        latitudeStepWidget->setLayout(latitudeStepLayout);

        auto *planetColumnsLayout = new QHBoxLayout();
        planetColumnsLayout->addLayout(planetLeftFormLayout);
        planetColumnsLayout->addLayout(planetRightFormLayout);

        auto *planetControlsLayout = new QFormLayout();
        planetControlsLayout->addRow(QStringLiteral("Материал поверхности:"), materialComboBox_);
        planetControlsLayout->addRow(QStringLiteral("Режим вращения:"), rotationModeWidget);
        planetControlsLayout->addRow(QStringLiteral("Семя рельефа:"), heightSeedSpinBox_);
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
        surfaceViewStack_->setCurrentWidget(surfaceMapWidget_);
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
        surfaceMapModeComboBox_->addItem(QStringLiteral("Температура"),
                                         static_cast<int>(SurfaceMapMode::Temperature));
        surfaceMapModeComboBox_->addItem(QStringLiteral("Высота"),
                                         static_cast<int>(SurfaceMapMode::Height));
        surfaceMapModeComboBox_->addItem(QStringLiteral("Ветер"),
                                         static_cast<int>(SurfaceMapMode::Wind));
        surfaceViewToggleButton_ = new QPushButton(QStringLiteral("3D вид"), this);
        surfaceViewToggleButton_->setCheckable(true);
        subsurfaceLayersSpinBox_ = new QSpinBox(this);
        subsurfaceLayersSpinBox_->setRange(1, 200);
        subsurfaceLayersSpinBox_->setValue(24);
        subsurfaceTopThicknessSpinBox_ = new QDoubleSpinBox(this);
        subsurfaceTopThicknessSpinBox_->setRange(0.001, 10.0);
        subsurfaceTopThicknessSpinBox_->setDecimals(3);
        subsurfaceTopThicknessSpinBox_->setSingleStep(0.01);
        subsurfaceTopThicknessSpinBox_->setValue(0.02);
        subsurfaceTopThicknessSpinBox_->setSuffix(QStringLiteral(" м"));
        subsurfaceDepthSpinBox_ = new QDoubleSpinBox(this);
        subsurfaceDepthSpinBox_->setRange(0.1, 200.0);
        subsurfaceDepthSpinBox_->setDecimals(2);
        subsurfaceDepthSpinBox_->setSingleStep(0.1);
        subsurfaceDepthSpinBox_->setValue(2.0);
        subsurfaceDepthSpinBox_->setSuffix(QStringLiteral(" м"));
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
        surfaceLegendScaleStack_ = new QStackedWidget(this);
        surfaceLegendScaleStack_->addWidget(temperatureScaleWidget_);
        surfaceLegendScaleStack_->addWidget(heightScaleWidget_);
        surfaceLegendScaleStack_->addWidget(windScaleWidget_);
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

        connect(planetComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int) {
            cancelTemperatureCalculation();
            updatePlanetSemiMajorAxisLabel();
            updatePlanetDayLengthLabel();
            updatePlanetMassLabel();
            updatePlanetRadiusLabel();
            updatePlanetDerivedLabels();
            updateAtmospherePlanetParameters();
            updateAtmosphereComposition();
            updatePlanetOrbitLabels();
            updateLatitudePointsDefault();
            syncMaterialWithPlanet();
            syncRotationModeWithPlanet();
            syncHeightSeedWithPlanet();
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
            clearTemperatureCache();
            updateTemperaturePlot();
        });

        connect(materialComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int) {
            syncPlanetMaterialWithSelection();
            clearTemperatureCache();
            updateTemperaturePlot();
        });

        connect(rotationModeComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
                [this](int) {
            syncPlanetRotationModeWithSelection();
            updateRotationModeIllustration();
            clearTemperatureCache();
            updateTemperaturePlot();
        });

        connect(heightSeedSpinBox_, QOverload<int>::of(&QSpinBox::valueChanged), this, [this](int) {
            if (planetComboBox_->currentIndex() < 0) {
                return;
            }
            syncPlanetHeightSeedWithSelection();
            updateSurfaceGridTemperatures();
        });

        connect(latitudeStepFastRadio_, &QRadioButton::toggled, this, [this](bool checked) {
            if (!checked) {
                return;
            }
            latitudePointsManuallySet_ = true;
            clearTemperatureCache();
            updateTemperaturePlot();
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
            clearTemperatureCache();
            updateTemperaturePlot();
            updateSurfaceGridTemperatures();
        };
        connect(subsurfaceLayersSpinBox_, QOverload<int>::of(&QSpinBox::valueChanged), this,
                [onSubsurfaceChanged](int) { onSubsurfaceChanged(); });
        connect(subsurfaceTopThicknessSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [onSubsurfaceChanged](double) { onSubsurfaceChanged(); });
        connect(subsurfaceDepthSpinBox_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
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

        connect(latitudeStepSlowRadio_, &QRadioButton::toggled, this, [this](bool checked) {
            if (!checked) {
                return;
            }
            latitudePointsManuallySet_ = true;
            clearTemperatureCache();
            updateTemperaturePlot();
        });

        connect(segmentSelectorWidget_, &SegmentSelectorWidget::currentIndexChanged, this,
                [this](int) { updateTemperaturePlotForSelectedSegment(); });

        // Сглаживание влияет только на отображение кривых, а не на физический расчет.
        // connect(smoothingCheckBox, &QCheckBox::toggled, temperaturePlot_,
        //         &SurfaceTemperaturePlot::setSmoothingEnabled);

        connect(calculateButton, &QPushButton::clicked, this, [this]() {
            autoCalculateEnabled_ = true;
            onCalculateRequested();
        });

        applyPrimary(StellarParameters{1.0, 5772.0, 1.0});
        applySecondary(std::nullopt);
        setPlanetPresets(solarSystemPresets(), QStringLiteral("Земля"));
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
            clearTemperatureCache();
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
    ModeIllustrationWidget *modeIllustrationWidget_ = nullptr;
    QRadioButton *latitudeStepFastRadio_ = nullptr;
    QRadioButton *latitudeStepSlowRadio_ = nullptr;
    QButtonGroup *latitudeStepGroup_ = nullptr;
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
    QSpinBox *subsurfaceLayersSpinBox_ = nullptr;
    QDoubleSpinBox *subsurfaceTopThicknessSpinBox_ = nullptr;
    QDoubleSpinBox *subsurfaceDepthSpinBox_ = nullptr;
    QComboBox *subsurfaceBoundaryComboBox_ = nullptr;
    SurfaceTemperatureScaleWidget *temperatureScaleWidget_ = nullptr;
    SurfaceHeightScaleWidget *heightScaleWidget_ = nullptr;
    SurfaceWindScaleWidget *windScaleWidget_ = nullptr;
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
    QVector<OrbitSegment> lastOrbitSegments_;
    QVector<QVector<TemperatureRangePoint>> lastTemperatureSegments_;
    QVector<QVector<TemperatureRangePoint>> lastTemperatureSegmentsSurfaceOnly_;
    QVector<TemperatureSummaryPoint> temperatureSummary_;
    QVector<TemperatureSummaryPoint> temperatureSummarySurfaceOnly_;
    bool lastTemperatureUsesAtmosphere_ = false;
    double surfaceMinTemperatureK_ = 0.0;
    double surfaceMaxTemperatureK_ = 0.0;
    bool hasSurfaceTemperatureRange_ = false;
    double surfaceMinHeightKm_ = 0.0;
    double surfaceMaxHeightKm_ = 0.0;
    bool hasSurfaceHeightRange_ = false;
    double surfaceMinWindSpeedMps_ = 0.0;
    double surfaceMaxWindSpeedMps_ = 0.0;
    bool hasSurfaceWindRange_ = false;
    SurfaceMapMode surfaceMapMode_ = SurfaceMapMode::Temperature;
    bool latitudePointsManuallySet_ = false;
    bool autoCalculateEnabled_ = false;
    double surfaceSimSpeedMultiplier_ = 1.0;
    QHash<TemperatureCacheKey, TemperatureCacheEntry> temperatureCache_;
    QHash<TemperatureCacheKey, TemperatureCacheEntry> temperatureCacheSurfaceOnly_;
    std::optional<StellarCacheKey> lastStellarKey_;
    bool surfaceSimRunning_ = false;
    struct SurfaceSimulationState {
        int dayIndex = 0;
        int hourIndex = 0;
    } surfaceSimState_;
    OrbitAnimationModel surfaceOrbitAnimation_;
    bool surfaceOrbitAnimationInitialized_ = false;

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
        clearTemperatureCache();
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
        updateLatitudePointsDefault();
        syncMaterialWithPlanet();
        syncRotationModeWithPlanet();
        syncHeightSeedWithPlanet();
        updatePlanetActions();
    }

    void clearPlanetPresets() {
        const QSignalBlocker blocker(planetComboBox_);
        planetComboBox_->clear();
        presetPlanetNames_.clear();
        clearTemperatureCache();
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

    QString formatSegmentLabel(const OrbitSegment &segment) const {
        const double meanAnomalyDegrees = qRadiansToDegrees(segment.meanAnomalyRadians);
        return QStringLiteral("Сегмент %1 (M=%2°, r=%3 а.е.)")
            .arg(segment.index + 1)
            .arg(QLocale().toString(meanAnomalyDegrees, 'f', 0))
            .arg(formatDistance(segment.distanceAU));
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
        if (latitudeStepSlowRadio_ && latitudeStepSlowRadio_->isChecked()) {
            return 10;
        }
        return 1;
    }

    int latitudePoints() const {
        return 180 / latitudeStepDegrees() + 1;
    }

    void updateLatitudePointsDefault() {
        if (latitudePointsManuallySet_) {
            return;
        }

        const QVariant value = planetComboBox_->currentData(kRoleDayLength);
        if (!value.isValid()) {
            return;
        }

        const double dayLength = value.toDouble();
        // Для медленных планет увеличиваем шаг широты, чтобы профили быстро считались
        // и оставались читаемыми при малом числе характерных широт.
        const bool useSlowStep = (dayLength > 30.0);
        const QSignalBlocker fastBlocker(latitudeStepFastRadio_);
        const QSignalBlocker slowBlocker(latitudeStepSlowRadio_);
        if (useSlowStep) {
            latitudeStepSlowRadio_->setChecked(true);
        } else {
            latitudeStepFastRadio_->setChecked(true);
        }
    }

    void updatePlanetOrbitLabels() {
        const QVariant eccentricity = planetComboBox_->currentData(kRoleEccentricity);
        const QVariant obliquity = planetComboBox_->currentData(kRoleObliquity);
        const QVariant perihelionArgument = planetComboBox_->currentData(kRolePerihelionArgument);
        if (!eccentricity.isValid() || !obliquity.isValid() || !perihelionArgument.isValid()) {
            planetEccentricityLabel_->setText(QStringLiteral("—"));
            planetObliquityLabel_->setText(QStringLiteral("—"));
            planetPerihelionArgumentLabel_->setText(QStringLiteral("—"));
            return;
        }
        planetEccentricityLabel_->setText(formatEccentricity(eccentricity.toDouble()));
        planetObliquityLabel_->setText(formatAngle(obliquity.toDouble()));
        planetPerihelionArgumentLabel_->setText(formatAngle(perihelionArgument.toDouble()));
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
        planetComboBox_->setItemData(index, static_cast<int>(planet.heightSourceType),
                                     kRoleHeightSourceType);
        planetComboBox_->setItemData(index, planet.heightmapPath, kRoleHeightmapPath);
        planetComboBox_->setItemData(index, planet.heightmapScaleKm, kRoleHeightmapScaleKm);
        planetComboBox_->setItemData(index, planet.heightSeed, kRoleHeightSeed);
        planetComboBox_->setItemData(index, planet.useContinentsHeight, kRoleUseContinentsHeight);
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
                 greenhouseOpacityInput, heightSeedInput, materialInput, rotationModeInput,
                 atmosphereInput, this]() {
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

            const int existingIndex = findPlanetIndexByName(name);
            const QString materialId = materialInput->currentData().toString();
            const RotationMode rotationMode =
                static_cast<RotationMode>(rotationModeInput->currentData().toInt());
            const bool tidallyLocked = (rotationMode == RotationMode::TidalLocked);
            const quint32 heightSeed = static_cast<quint32>(heightSeedInput->value());
            const AtmosphereComposition composition = atmosphereInput->composition(false);
            PlanetPreset preset{name, axis, dayLength, eccentricity, obliquity,
                                perihelionArgument, massEarths, radiusKm, materialId,
                                composition, greenhouseOpacity, tidallyLocked};
            preset.heightSeed = heightSeed;
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
                                             static_cast<int>(preset.heightSourceType),
                                             kRoleHeightSourceType);
                planetComboBox_->setItemData(existingIndex, preset.heightmapPath, kRoleHeightmapPath);
                planetComboBox_->setItemData(existingIndex, preset.heightmapScaleKm,
                                             kRoleHeightmapScaleKm);
                planetComboBox_->setItemData(existingIndex, preset.heightSeed, kRoleHeightSeed);
                planetComboBox_->setItemData(existingIndex, preset.useContinentsHeight,
                                             kRoleUseContinentsHeight);
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
            clearTemperatureCache();
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

    void updateSegmentComboBox() {
        const QSignalBlocker blocker(segmentSelectorWidget_);
        const int previousIndex = segmentSelectorWidget_->currentIndex();

        segmentSelectorWidget_->setSegments(lastOrbitSegments_);

        if (lastOrbitSegments_.isEmpty()) {
            segmentSelectorWidget_->setEnabled(false);
            return;
        }

        if (previousIndex >= 0 && previousIndex < lastOrbitSegments_.size()) {
            segmentSelectorWidget_->setCurrentIndex(previousIndex);
        } else {
            segmentSelectorWidget_->setCurrentIndex(0);
        }
    }

    void rebuildTemperatureSummaryFromSegments(
        const QVector<QVector<TemperatureRangePoint>> &segments,
        QVector<TemperatureSummaryPoint> *summary) {
        summary->clear();
        if (segments.isEmpty()) {
            return;
        }

        const int latitudeCount = segments.front().size();
        summary->reserve(latitudeCount);

        for (int i = 0; i < latitudeCount; ++i) {
            double minimumOverall = std::numeric_limits<double>::max();
            double maximumOverall = 0.0;
            double meanAnnualSum = 0.0;
            double meanAnnualDaySum = 0.0;
            double meanAnnualNightSum = 0.0;
            int samples = 0;
            double latitudeDegrees = 0.0;

            for (const auto &segment : segments) {
                if (i >= segment.size()) {
                    continue;
                }
                const auto &point = segment.at(i);
                latitudeDegrees = point.latitudeDegrees;
                minimumOverall = qMin(minimumOverall, point.minimumKelvin);
                maximumOverall = qMax(maximumOverall, point.maximumKelvin);
                // Годовая средняя берется как среднее по сегментам равной длительности.
                meanAnnualSum += point.meanDailyKelvin;
                // Дневная/ночная годовая средняя получаются так же по сегментам.
                meanAnnualDaySum += point.meanDayKelvin;
                meanAnnualNightSum += point.meanNightKelvin;
                ++samples;
            }

            if (samples == 0) {
                continue;
            }

            TemperatureSummaryPoint summaryPoint;
            summaryPoint.latitudeDegrees = latitudeDegrees;
            summaryPoint.minimumKelvin = minimumOverall;
            summaryPoint.maximumKelvin = maximumOverall;
            summaryPoint.meanAnnualKelvin = meanAnnualSum / samples;
            summaryPoint.meanAnnualDayKelvin = meanAnnualDaySum / samples;
            summaryPoint.meanAnnualNightKelvin = meanAnnualNightSum / samples;
            summaryPoint.minimumCelsius = minimumOverall - kKelvinOffset;
            summaryPoint.maximumCelsius = maximumOverall - kKelvinOffset;
            summaryPoint.meanAnnualCelsius = summaryPoint.meanAnnualKelvin - kKelvinOffset;
            summaryPoint.meanAnnualDayCelsius = summaryPoint.meanAnnualDayKelvin - kKelvinOffset;
            summaryPoint.meanAnnualNightCelsius =
                summaryPoint.meanAnnualNightKelvin - kKelvinOffset;
            summary->push_back(summaryPoint);
        }
    }

    void rebuildTemperatureSummary() {
        rebuildTemperatureSummaryFromSegments(lastTemperatureSegments_, &temperatureSummary_);
    }

    void rebuildSurfaceOnlyTemperatureSummary() {
        rebuildTemperatureSummaryFromSegments(lastTemperatureSegmentsSurfaceOnly_,
                                              &temperatureSummarySurfaceOnly_);
    }

    void updateTemperaturePlotForSelectedSegment() {
        if (segmentSelectorWidget_->currentIndex() < 0 ||
            segmentSelectorWidget_->currentIndex() >= lastTemperatureSegments_.size()) {
            temperaturePlot_->clearSeries();
            return;
        }

        const int index = segmentSelectorWidget_->currentIndex();
        const QString label = formatSegmentLabel(lastOrbitSegments_.at(index));
        const RotationMode rotationMode =
            static_cast<RotationMode>(planetComboBox_->currentData(kRoleRotationMode).toInt());
        temperaturePlot_->setTemperatureSeries(lastTemperatureSegments_.at(index),
                                               temperatureSummary_,
                                               label,
                                               rotationMode,
                                               lastTemperatureUsesAtmosphere_);
    }

    void clearTemperatureSegments() {
        lastOrbitSegments_.clear();
        lastTemperatureSegments_.clear();
        lastTemperatureSegmentsSurfaceOnly_.clear();
        temperatureSummary_.clear();
        temperatureSummarySurfaceOnly_.clear();
        lastTemperatureUsesAtmosphere_ = false;
        segmentSelectorWidget_->setSegments({});
        segmentSelectorWidget_->setEnabled(false);
        temperaturePlot_->clearSeries();
    }

    void clearTemperatureCache() {
        temperatureCache_.clear();
        temperatureCacheSurfaceOnly_.clear();
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
        surfaceGrid_.setHeightSource(heightSource, heightmapPath, heightmapScaleKm,
                                     heightSeed, useContinentsHeight);

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

    void applySurfaceMapMode(SurfaceMapMode mode) {
        surfaceMapMode_ = mode;
        if (surfaceMapWidget_) {
            surfaceMapWidget_->setMapMode(mode);
        }
        if (surfaceGlobeWidget_) {
            surfaceGlobeWidget_->setMapMode(mode);
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
            updateSurfaceWindLegend(false, 0.0, 0.0);
            return;
        }

        const auto stateDefaults = buildSurfacePointStateDefaults();
        if (!stateDefaults) {
            applySurfaceGridToViews();
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceWindLegend(false, 0.0, 0.0);
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

        // Вкладка «Поверхность» всегда использует модель поверхности без атмосферных поправок.
        const QVector<TemperatureSummaryPoint> *summarySource = nullptr;
        if (!temperatureSummarySurfaceOnly_.isEmpty()) {
            summarySource = &temperatureSummarySurfaceOnly_;
        }

        const QVector<TemperatureRangePoint> *segmentSource = nullptr;
        if (!summarySource && !lastTemperatureSegmentsSurfaceOnly_.isEmpty()) {
            segmentSource = &lastTemperatureSegmentsSurfaceOnly_.first();
        }

        if (!summarySource && !segmentSource) {
            for (auto &point : surfaceGrid_.points()) {
                const SurfaceMaterial material = materialForPoint(point.materialId);
                // Альбедо и теплофизика теперь зависят от материала ячейки.
                const double albedo = qBound(0.0, material.albedo, 1.0);
                point.temperatureK = stateDefaults->minTemperatureKelvin;
                point.state = SurfacePointState(point.temperatureK,
                                                albedo,
                                                stateDefaults->greenhouseOpacity,
                                                stateDefaults->minTemperatureKelvin,
                                                material,
                                                stateDefaults->subsurfaceSettings);
            }
            applySurfaceGridToViews();
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
            updateSurfaceWindLegend(false, 0.0, 0.0);
            return;
        }

        const double dayLengthDays = planetComboBox_->currentData(kRoleDayLength).toDouble();
        const RotationMode rotationMode =
            static_cast<RotationMode>(planetComboBox_->currentData(kRoleRotationMode).toInt());
        const int stepsPerDay = qMax(1, qRound(dayLengthDays * 24.0));
        const double phase =
            2.0 * M_PI *
            (static_cast<double>(surfaceSimState_.hourIndex) + 0.5) /
            static_cast<double>(stepsPerDay);
        const double baseHourAngle =
            (rotationMode == RotationMode::TidalLocked) ? 0.0 : (phase - M_PI);
        // Субзвёздная долгота - долгота, где часовой угол равен нулю.
        // Для приливной синхронизации она фиксирована в системе долгот карты.
        const double substellarLongitudeRadians =
            (rotationMode == RotationMode::TidalLocked) ? 0.0 : -baseHourAngle;

        ensureSurfaceOrbitAnimationReady();
        const double declinationDegrees = surfaceOrbitAnimation_.declinationDegrees();
        const double declinationRadians = qDegreesToRadians(declinationDegrees);
        const double sinDeclination = std::sin(declinationRadians);
        const double cosDeclination = std::cos(declinationRadians);

        double atmospherePressureAtm = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            atmospherePressureAtm = atmosphere.totalPressureAtm(massEarths, radiusKm);
        }
        const double transport =
            (atmospherePressureAtm > 50.0)
                ? 0.99
                : (atmospherePressureAtm > 0.001
                       ? qMin(1.0, 0.15 * std::log(atmospherePressureAtm * 100.0 + 1.0))
                       : 0.0);
        const double rotBlock =
            (dayLengthDays < 2.0 && atmospherePressureAtm < 10.0) ? 0.65 : 1.0;
        const double meridionalTransport = transport * rotBlock;

        double segmentSolarConstant = lastSolarConstant_;
        if (lastSolarConstantDistanceAU_ > 0.0) {
            const double distanceAU = surfaceOrbitAnimation_.distanceAU();
            segmentSolarConstant =
                lastSolarConstant_ *
                std::pow(lastSolarConstantDistanceAU_ / distanceAU, 2.0);
        }
        // Глобальный средний поток перед альбедо, как в SurfaceTemperatureCalculator.
        const double globalAverageInsolation = segmentSolarConstant / 4.0;
        const double timeStepSeconds = 3600.0;

        auto interpolateBaselineTemperature = [summarySource, segmentSource](double latitudeDeg) {
            // Базовая температура зависит от широты и используется как старт для мгновенного шага.
            if (summarySource) {
                const auto &points = *summarySource;
                auto interpolate = [&points, latitudeDeg](auto valueForPoint) {
                    if (points.size() == 1) {
                        return valueForPoint(points.first());
                    }
                    if (latitudeDeg <= points.first().latitudeDegrees) {
                        return valueForPoint(points.first());
                    }
                    if (latitudeDeg >= points.last().latitudeDegrees) {
                        return valueForPoint(points.last());
                    }
                    for (int i = 1; i < points.size(); ++i) {
                        if (latitudeDeg <= points[i].latitudeDegrees) {
                            const auto &lower = points[i - 1];
                            const auto &upper = points[i];
                            const double span = upper.latitudeDegrees - lower.latitudeDegrees;
                            const double t =
                                span > 0.0 ? (latitudeDeg - lower.latitudeDegrees) / span : 0.0;
                            return valueForPoint(lower) + t * (valueForPoint(upper) - valueForPoint(lower));
                        }
                    }
                    return valueForPoint(points.last());
                };
                return interpolate([](const TemperatureSummaryPoint &point) {
                    return point.meanAnnualKelvin;
                });
            }

            const auto &points = *segmentSource;
            auto interpolate = [&points, latitudeDeg](auto valueForPoint) {
                if (points.size() == 1) {
                    return valueForPoint(points.first());
                }
                if (latitudeDeg <= points.first().latitudeDegrees) {
                    return valueForPoint(points.first());
                }
                if (latitudeDeg >= points.last().latitudeDegrees) {
                    return valueForPoint(points.last());
                }
                for (int i = 1; i < points.size(); ++i) {
                    if (latitudeDeg <= points[i].latitudeDegrees) {
                        const auto &lower = points[i - 1];
                        const auto &upper = points[i];
                        const double span = upper.latitudeDegrees - lower.latitudeDegrees;
                        const double t =
                            span > 0.0 ? (latitudeDeg - lower.latitudeDegrees) / span : 0.0;
                        return valueForPoint(lower) + t * (valueForPoint(upper) - valueForPoint(lower));
                    }
                }
                return valueForPoint(points.last());
            };
            return interpolate([](const TemperatureRangePoint &point) {
                return point.meanDailyKelvin;
            });
        };

        double minTemperature = std::numeric_limits<double>::max();
        double maxTemperature = std::numeric_limits<double>::lowest();
        for (auto &point : surfaceGrid_.points()) {
            const SurfaceMaterial material = materialForPoint(point.materialId);
            // Альбедо и теплофизика теперь зависят от материала ячейки.
            const double albedo = qBound(0.0, material.albedo, 1.0);
            const double baselineTemperature =
                interpolateBaselineTemperature(point.latitudeDeg);
            const double localHourAngle = point.longitudeRadians - substellarLongitudeRadians;
            const double cosZenith =
                point.sinLatitude * sinDeclination +
                point.cosLatitude * cosDeclination * std::cos(localHourAngle);
            // Синхронизация формул инициализации и шага нужна, чтобы при старте
            // не было скачков температуры и мерцания на освещенной стороне.
            const double localInsolation =
                segmentSolarConstant * qMax(0.0, cosZenith);
            const double blendedInsolation =
                localInsolation * (1.0 - meridionalTransport) +
                globalAverageInsolation * meridionalTransport;

            SurfacePointState state(baselineTemperature,
                                    albedo,
                                    stateDefaults->greenhouseOpacity,
                                    stateDefaults->minTemperatureKelvin,
                                    material,
                                    stateDefaults->subsurfaceSettings);
            const double absorbedFlux = state.absorbedFlux(blendedInsolation);
            const double emittedFlux = state.emittedFlux();
            state.updateTemperature(absorbedFlux, emittedFlux, timeStepSeconds);
            point.state = state;
            // Обновляем поверхностную температуру для карт и переноса.
            point.temperatureK = point.state.temperatureKelvin();
            minTemperature = qMin(minTemperature, point.temperatureK);
            maxTemperature = qMax(maxTemperature, point.temperatureK);
        }

        double surfaceGravity = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            const double radiusMeters = radiusKm * 1000.0;
            const double planetMassKg = massEarths * kEarthMassKg;
            surfaceGravity = kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);
        }
        const double gravity = (surfaceGravity > 0.0) ? surfaceGravity : 9.80665;
        for (auto &point : surfaceGrid_.points()) {
            // Инициализируем поверхностное давление (на уровне рельефа), чтобы дальше
            // переносить его ветром без пересчёта из всей атмосферы каждый тик.
            const double pressureAtm =
                AtmosphericPressureModel::pressureAtHeightAtm(atmospherePressureAtm,
                                                              point.heightKm * 1000.0,
                                                              point.temperatureK,
                                                              atmosphere,
                                                              gravity);
            // Храним поверхностное давление в точке для последующего переноса по ветру.
            point.pressureAtm = qMax(0.0, pressureAtm);
        }
        updateSurfaceWindField(atmosphere, atmospherePressureAtm, dayLengthDays, surfaceGravity);

        applySurfaceGridToViews();
        if (minTemperature <= maxTemperature) {
            applySurfaceTemperatureRangeToViews(minTemperature, maxTemperature);
            updateSurfaceTemperatureLegend(true, minTemperature, maxTemperature);
        } else {
            updateSurfaceTemperatureLegend(false, 0.0, 0.0);
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

        const auto material = currentMaterial();
        if (!material) {
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
        if (advectedPressures.size() == surfaceGrid_.points().size()) {
            for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
                // Фиксируем перенесённое поверхностное давление для UI и динамики.
                surfaceGrid_.points()[i].pressureAtm = advectedPressures.at(i);
            }
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
        for (int i = 0; i < surfaceGrid_.points().size(); ++i) {
            auto &point = surfaceGrid_.points()[i];
            point.state.setTemperatureKelvin(advectedTemperatures.at(i));
            // Температура после переноса хранится как поверхностная величина.
            point.temperatureK = point.state.temperatureKelvin();
            minTemperature = qMin(minTemperature, point.temperatureK);
            maxTemperature = qMax(maxTemperature, point.temperatureK);
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

        if (!hasSolarConstant_) {
            clearTemperatureSegments();
            return;
        }

        if (planetComboBox_->currentIndex() < 0) {
            clearTemperatureSegments();
            return;
        }

        const double semiMajorAxis = planetComboBox_->currentData(kRoleSemiMajorAxis).toDouble();
        if (semiMajorAxis <= 0.0) {
            clearTemperatureSegments();
            return;
        }

        // referenceDistanceAU — расстояние, на котором была рассчитана lastSolarConstant_.
        const double referenceDistanceAU = lastSolarConstantDistanceAU_;
        if (referenceDistanceAU <= 0.0) {
            clearTemperatureSegments();
            return;
        }

        const double dayLength = planetComboBox_->currentData(kRoleDayLength).toDouble();
        if (dayLength <= 0.0) {
            clearTemperatureSegments();
            return;
        }

        const double eccentricity = planetComboBox_->currentData(kRoleEccentricity).toDouble();
        const double obliquity = planetComboBox_->currentData(kRoleObliquity).toDouble();
        const double perihelionArgument =
            planetComboBox_->currentData(kRolePerihelionArgument).toDouble();

        const auto material = currentMaterial();
        if (!material) {
            clearTemperatureSegments();
            return;
        }

        const SubsurfaceModelSettings subsurfaceSettings = buildSubsurfaceSettings();
        const int latitudePointCount = latitudePoints();
        const RotationMode rotationMode =
            static_cast<RotationMode>(planetComboBox_->currentData(kRoleRotationMode).toInt());
        AtmosphereComposition atmosphere;
        const QVariant atmosphereValue = planetComboBox_->currentData(kRoleAtmosphere);
        if (atmosphereValue.isValid()) {
            atmosphere = atmosphereValue.value<AtmosphereComposition>();
        }
        const double massEarths = planetComboBox_->currentData(kRoleMassEarths).toDouble();
        const double radiusKm = planetComboBox_->currentData(kRoleRadiusKm).toDouble();
        const double greenhouseOpacity =
            planetComboBox_->currentData(kRoleGreenhouseOpacity).toDouble();
        const HeightSourceType heightSourceType =
            static_cast<HeightSourceType>(planetComboBox_->currentData(kRoleHeightSourceType)
                                              .toInt());
        const QString heightmapPath =
            planetComboBox_->currentData(kRoleHeightmapPath).toString();
        const double heightmapScaleKm =
            planetComboBox_->currentData(kRoleHeightmapScaleKm).toDouble();
        const quint32 heightSeed =
            planetComboBox_->currentData(kRoleHeightSeed).toUInt();
        const bool useContinentsHeight =
            planetComboBox_->currentData(kRoleUseContinentsHeight).toBool();
        double atmospherePressureAtm = 0.0;
        double surfaceGravity = 0.0;
        if (massEarths > 0.0 && radiusKm > 0.0) {
            atmospherePressureAtm = atmosphere.totalPressureAtm(massEarths, radiusKm);
            const double radiusMeters = radiusKm * 1000.0;
            const double planetMassKg = massEarths * kEarthMassKg;
            // g = G * M / R^2.
            surfaceGravity =
                kGravitationalConstant * planetMassKg / (radiusMeters * radiusMeters);
        }
        // Атмосферу считаем присутствующей, если есть масса газов или ненулевое давление:
        // это позволяет явно переключать формулы поверхности/атмосферы.
        const bool hasAtmosphere =
            atmosphere.totalMassGigatons() > 0.0 || atmospherePressureAtm > 0.0;
        const int segmentCount = 12;
        const TemperatureCacheKey cacheKey{lastSolarConstant_,
                                            material->id,
                                            atmosphereSignature(atmosphere),
                                            atmospherePressureAtm,
                                            surfaceGravity,
                                            greenhouseOpacity,
                                            dayLength,
                                            referenceDistanceAU,
                                            semiMajorAxis,
                                            eccentricity,
                                            obliquity,
                                            perihelionArgument,
                                            // Радиус влияет на столбовую плотность и водный баланс, поэтому он
                                            // должен входить в ключ кэша температурных расчётов.
                                            radiusKm,
                                            heightSourceType,
                                            heightmapPath,
                                            heightmapScaleKm,
                                            heightSeed,
                                            useContinentsHeight,
                                            subsurfaceSettings.layerCount,
                                            subsurfaceSettings.topLayerThicknessMeters,
                                            subsurfaceSettings.bottomDepthMeters,
                                            subsurfaceSettings.bottomBoundary,
                                            latitudePointCount,
                                            segmentCount,
                                            rotationMode};
        const AtmosphereComposition surfaceOnlyAtmosphere;
        const TemperatureCacheKey surfaceOnlyCacheKey{lastSolarConstant_,
                                                      material->id,
                                                      atmosphereSignature(surfaceOnlyAtmosphere),
                                                      0.0,
                                                      surfaceGravity,
                                                      0.0,
                                                      dayLength,
                                                      referenceDistanceAU,
                                                      semiMajorAxis,
                                                      eccentricity,
                                                      obliquity,
                                                      perihelionArgument,
                                                      // Радиус влияет на столбовую плотность и водный баланс, поэтому он
                                                      // должен входить в ключ кэша температурных расчётов.
                                                      radiusKm,
                                                      heightSourceType,
                                                      heightmapPath,
                                                      heightmapScaleKm,
                                                      heightSeed,
                                                      useContinentsHeight,
                                                      subsurfaceSettings.layerCount,
                                                      subsurfaceSettings.topLayerThicknessMeters,
                                                      subsurfaceSettings.bottomDepthMeters,
                                                      subsurfaceSettings.bottomBoundary,
                                                      latitudePointCount,
                                                      segmentCount,
                                                      rotationMode};
        auto startSurfaceOnlyCalculation =
            [this,
             material,
             dayLength,
             rotationMode,
             surfaceOnlyAtmosphere,
             surfaceGravity,
             radiusKm,
             subsurfaceSettings,
             heightSourceType,
             heightmapPath,
             heightmapScaleKm,
             heightSeed,
             useContinentsHeight,
             latitudePointCount,
             referenceDistanceAU,
             obliquity,
             perihelionArgument,
             surfaceOnlyCacheKey](int requestId,
                                  const std::shared_ptr<std::atomic_bool> &cancelFlag,
                                  const QVector<OrbitSegment> &segments) {
                SurfaceTemperatureCalculator surfaceOnlyCalculator(lastSolarConstant_,
                                                                  *material,
                                                                  dayLength,
                                                                  rotationMode,
                                                                  surfaceOnlyAtmosphere,
                                                                  0.0,
                                                                  0.0,
                                                                  surfaceGravity,
                                                                  radiusKm,
                                                                  false,
                                                                  8,
                                                                  heightSourceType,
                                                                  heightmapPath,
                                                                  heightmapScaleKm,
                                                                  heightSeed,
                                                                  useContinentsHeight,
                                                                  subsurfaceSettings);
                startTemperatureElapsedUi(requestId, QPointer<QProgressDialog>());
                auto *surfaceWatcher =
                    new QFutureWatcher<QVector<QVector<TemperatureRangePoint>>>(this);
                connect(surfaceWatcher,
                        &QFutureWatcher<QVector<QVector<TemperatureRangePoint>>>::finished,
                        this,
                        [this, surfaceWatcher, requestId, cancelFlag, surfaceOnlyCacheKey]() {
                            surfaceWatcher->deleteLater();
                            const bool shouldFinalizeElapsed =
                                requestId == temperatureRequestId_ &&
                                (!temperatureProgressDialog_ ||
                                 !temperatureProgressDialog_->isVisible());
                            if (shouldFinalizeElapsed && temperatureUiTimer_) {
                                temperatureUiTimer_->stop();
                            }
                            if (requestId != temperatureRequestId_ || cancelFlag->load()) {
                                return;
                            }
                            const auto result = surfaceWatcher->result();
                            if (result.isEmpty()) {
                                return;
                            }
                            lastTemperatureSegmentsSurfaceOnly_ = result;
                            temperatureCacheSurfaceOnly_.insert(
                                surfaceOnlyCacheKey,
                                TemperatureCacheEntry{lastOrbitSegments_, result});
                            rebuildSurfaceOnlyTemperatureSummary();
                            updateSurfaceGridTemperatures();
                            if (shouldFinalizeElapsed && temperatureElapsedLabel_) {
                                const qint64 elapsedSeconds =
                                    temperatureElapsed_.elapsed() / 1000;
                                const int minutes = static_cast<int>(elapsedSeconds / 60);
                                const int seconds = static_cast<int>(elapsedSeconds % 60);
                                const QString pauseSuffix = temperaturePauseFlag_.load()
                                    ? QStringLiteral(" (пауза)")
                                    : QString();
                                temperatureElapsedLabel_->setText(
                                    QStringLiteral("Прошло: %1:%2%3")
                                        .arg(minutes, 2, 10, QLatin1Char('0'))
                                        .arg(seconds, 2, 10, QLatin1Char('0'))
                                        .arg(pauseSuffix));
                            }
                        });

                const std::atomic_bool *pauseFlag = &temperaturePauseFlag_;
                auto mapSegmentSurfaceOnly =
                    [surfaceOnlyCalculator,
                     latitudePointCount,
                     cancelFlag,
                     pauseFlag,
                     referenceDistanceAU,
                     obliquity,
                     perihelionArgument](const OrbitSegment &segment) {
                        if (cancelFlag && cancelFlag->load()) {
                            return QVector<TemperatureRangePoint>{};
                        }
                        auto noopProgress = [](int, int) {};
                        return surfaceOnlyCalculator.temperatureRangesForOrbitSegment(
                            segment,
                            referenceDistanceAU,
                            obliquity,
                            perihelionArgument,
                            latitudePointCount,
                            noopProgress,
                            cancelFlag.get(),
                            pauseFlag);
                    };

                surfaceWatcher->setFuture(QtConcurrent::run(
                    [segments, mapSegmentSurfaceOnly, cancelFlag]() {
                        QVector<QVector<TemperatureRangePoint>> results;
                        results.reserve(segments.size());
                        for (const auto &segment : segments) {
                            if (cancelFlag && cancelFlag->load()) {
                                break;
                            }
                            results.push_back(mapSegmentSurfaceOnly(segment));
                        }
                        return results;
                    }));
            };
        const auto cached = temperatureCache_.constFind(cacheKey);
        if (cached != temperatureCache_.constEnd()) {
            // Кэш нужен для быстрого переключения сегментов/планет без повторного расчёта.
            lastOrbitSegments_ = cached->orbitSegments;
            lastTemperatureSegments_ = cached->temperatureSegments;
            lastTemperatureUsesAtmosphere_ = hasAtmosphere;
            rebuildTemperatureSummary();
            updateSegmentComboBox();
            updateTemperaturePlotForSelectedSegment();
            if (!hasAtmosphere) {
                lastTemperatureSegmentsSurfaceOnly_ = lastTemperatureSegments_;
                rebuildSurfaceOnlyTemperatureSummary();
                temperatureCacheSurfaceOnly_.insert(
                    surfaceOnlyCacheKey,
                    TemperatureCacheEntry{lastOrbitSegments_, lastTemperatureSegments_});
            } else {
                const auto surfaceCached =
                    temperatureCacheSurfaceOnly_.constFind(surfaceOnlyCacheKey);
                if (surfaceCached != temperatureCacheSurfaceOnly_.constEnd()) {
                    lastTemperatureSegmentsSurfaceOnly_ = surfaceCached->temperatureSegments;
                    rebuildSurfaceOnlyTemperatureSummary();
                } else {
                    lastTemperatureSegmentsSurfaceOnly_.clear();
                    temperatureSummarySurfaceOnly_.clear();
                    const int requestId = ++temperatureRequestId_;
                    temperatureCancelFlag_ = std::make_shared<std::atomic_bool>(false);
                    const auto cancelFlag = temperatureCancelFlag_;
                    resetTemperatureUiState();
                    startSurfaceOnlyCalculation(requestId, cancelFlag, lastOrbitSegments_);
                }
            }
            updateSurfaceGridTemperatures();
            return;
        }

        lastTemperatureSegments_.clear();
        lastTemperatureSegmentsSurfaceOnly_.clear();
        temperatureSummarySurfaceOnly_.clear();

        const int requestId = ++temperatureRequestId_;
        temperatureCancelFlag_ = std::make_shared<std::atomic_bool>(false);
        const auto cancelFlag = temperatureCancelFlag_;
        const int totalLatitudes = latitudePointCount;
        const int totalProgress = totalLatitudes * segmentCount;

        resetTemperatureUiState();

        auto *progressDialog = ensureTemperatureProgressDialog();
        progressDialog->setRange(0, totalProgress);
        progressDialog->setValue(0);
        progressDialog->show();

        QPointer<QProgressDialog> dialogGuard(progressDialog);
        QPointer<SolarCalculatorWidget> widgetGuard(this);
        startTemperatureElapsedUi(requestId, dialogGuard);
        // Вызов из фонового потока: защищаемся от удаления виджета во время вычисления.
        auto progressCallback = [widgetGuard, dialogGuard, requestId](int processed, int total) {
            if (!widgetGuard) {
                return;
            }
            QMetaObject::invokeMethod(
                widgetGuard.data(),
                [widgetGuard, dialogGuard, processed, total, requestId]() {
                    if (!widgetGuard || !dialogGuard ||
                        requestId != widgetGuard->temperatureRequestId_) {
                        return;
                    }
                    dialogGuard->setMaximum(total);
                    dialogGuard->setValue(processed);
                },
                Qt::QueuedConnection);
        };

        // При наличии атмосферы включаем расширенную модель с парниковым слоем и циркуляцией.
        SurfaceTemperatureCalculator calculator(lastSolarConstant_, *material, dayLength,
                                                rotationMode, atmosphere, greenhouseOpacity,
                                                atmospherePressureAtm,
                                                surfaceGravity,
                                                radiusKm,
                                                hasAtmosphere,
                                                8,
                                                heightSourceType,
                                                heightmapPath,
                                                heightmapScaleKm,
                                                heightSeed,
                                                useContinentsHeight,
                                                subsurfaceSettings);
        auto *watcher = new QFutureWatcher<QVector<QVector<TemperatureRangePoint>>>(this);
        connect(watcher, &QFutureWatcher<QVector<QVector<TemperatureRangePoint>>>::finished, this,
                [this,
                 watcher,
                 requestId,
                 cancelFlag,
                 cacheKey,
                 surfaceOnlyCacheKey,
                 hasAtmosphere]() {
            watcher->deleteLater();
            if (requestId == temperatureRequestId_) {
                resetTemperatureUiState();
                temperatureProgressDialog_->hide();
            }
            if (requestId != temperatureRequestId_ || cancelFlag->load()) {
                return;
            }
            const auto result = watcher->result();
            if (result.isEmpty()) {
                return;
            }
            lastTemperatureSegments_ = result;
            temperatureCache_.insert(cacheKey, TemperatureCacheEntry{lastOrbitSegments_, result});
            rebuildTemperatureSummary();
            updateSegmentComboBox();
            updateTemperaturePlotForSelectedSegment();
            if (!hasAtmosphere) {
                lastTemperatureSegmentsSurfaceOnly_ = result;
                temperatureCacheSurfaceOnly_.insert(
                    surfaceOnlyCacheKey,
                    TemperatureCacheEntry{lastOrbitSegments_, result});
                rebuildSurfaceOnlyTemperatureSummary();
            }
            updateSurfaceGridTemperatures();
        });

        const OrbitSegmentCalculator orbitCalculator(semiMajorAxis, eccentricity);
        lastOrbitSegments_ = orbitCalculator.segments(segmentCount);
        updateSegmentComboBox();
        const auto segments = lastOrbitSegments_;
        const std::atomic_bool *pauseFlag = &temperaturePauseFlag_;

        auto progressCounter = std::make_shared<std::atomic_int>(0);
        auto segmentProgress = [progressCounter, progressCallback, totalProgress, cancelFlag](int, int) {
            if (cancelFlag && cancelFlag->load()) {
                return;
            }
            const int processed = progressCounter->fetch_add(1) + 1;
            progressCallback(processed, totalProgress);
        };

        auto mapSegment = [calculator,
                           latitudePointCount,
                           segmentProgress,
                           cancelFlag,
                           pauseFlag,
                           referenceDistanceAU,
                           obliquity,
                           perihelionArgument](const OrbitSegment &segment) {
            if (cancelFlag && cancelFlag->load()) {
                return QVector<TemperatureRangePoint>{};
            }
            return calculator.temperatureRangesForOrbitSegment(segment,
                                                               referenceDistanceAU,
                                                               obliquity,
                                                               perihelionArgument,
                                                               latitudePointCount,
                                                               segmentProgress,
                                                               cancelFlag.get(),
                                                               pauseFlag);
        };

        if (hasAtmosphere) {
            startSurfaceOnlyCalculation(requestId, cancelFlag, segments);
        }

        // В Qt 5.15 возникают проблемы с выводом типов в mappedReduced для сложных лямбд.
        // Запускаем фоновый расчет вручную, сохраняя порядок сегментов.
        watcher->setFuture(QtConcurrent::run([segments, mapSegment, cancelFlag]() {
            QVector<QVector<TemperatureRangePoint>> results;
            results.reserve(segments.size());
            for (const auto &segment : segments) {
                if (cancelFlag && cancelFlag->load()) {
                    break;
                }
                results.push_back(mapSegment(segment));
            }
            return results;
        }));
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
        clearTemperatureCache();
        updateTemperaturePlot();
    }
};
}  // namespace

ArgumentsParseResult parseParametersFromArguments(const QCoreApplication &app,
                                                  QTextStream &output,
                                                  BinarySystemParameters &parameters,
                                                  int &precision) {
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

    parser.addOptions({radiusOption, temperatureOption, distanceOption, radius2Option, temperature2Option,
                       distance2Option, precisionOption});
    parser.process(app);

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
    const ArgumentsParseResult argsResult = parseParametersFromArguments(app, output, parameters, precision);

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
