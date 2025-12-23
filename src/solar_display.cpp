#include "orbit_segment_calculator.h"
#include "planet_presets.h"
#include "solar_calculator.h"
#include "solar_display.h"
#include "surface_temperature_calculator.h"
#include "surface_temperature_plot.h"

#include <QtCore/QCommandLineOption>
#include <QtCore/QCommandLineParser>
#include <QtCore/QLocale>
#include <QtCore/QPointer>
#include <QtCore/QSignalBlocker>
#include <QtCore/QThread>
#include <QtCore/QThreadPool>
// #include <QtCore/QOverload>
#include <QtGlobal>
#include <QtCore/QtMath>
#include <QtCore/QSet>
#include <QtCore/QFutureWatcher>
#include <QtGui/QDoubleValidator>
#include <QtConcurrent/QtConcurrent>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QProgressDialog>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStyle>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

#include <atomic>
#include <cmath>
#include <functional>
#include <limits>

namespace {
constexpr int kRoleSemiMajorAxis = Qt::UserRole;
constexpr int kRoleIsCustom = Qt::UserRole + 1;
constexpr int kRolePlanetName = Qt::UserRole + 2;
constexpr int kRoleMaterialId = Qt::UserRole + 3;
constexpr int kRoleDayLength = Qt::UserRole + 4;
constexpr int kRoleEccentricity = Qt::UserRole + 5;
constexpr int kRoleObliquity = Qt::UserRole + 6;
constexpr int kRolePerihelionArgument = Qt::UserRole + 7;

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
            setPlanetPresets(solarSystemPresets());
        });

        addPresetButton(QStringLiteral("Сладкое Небо"), [this, applyPrimary, applySecondary]() {
            applyPrimary(StellarParameters{0.3761, 2576.0, 1.0});
            applySecondary(StellarParameters{0.3741, 2349.0, 1.0});
            resetSolarConstant();
            setPlanetPresets(sweetSkyPresets());
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
        connect(calculateButton, &QPushButton::clicked, this, [this]() { onCalculateRequested(); });

        resultLabel_ = new QLabel(
            QStringLiteral("Введите параметры и нажмите \"Рассчитать\"."), this);
        resultLabel_->setWordWrap(true);

        planetComboBox_ = new QComboBox(this);
        planetSemiMajorAxisLabel_ = new QLabel(QStringLiteral("—"), this);
        planetDayLengthLabel_ = new QLabel(QStringLiteral("—"), this);
        planetEccentricityLabel_ = new QLabel(QStringLiteral("—"), this);
        planetObliquityLabel_ = new QLabel(QStringLiteral("—"), this);
        planetPerihelionArgumentLabel_ = new QLabel(QStringLiteral("—"), this);
        materialComboBox_ = new QComboBox(this);
        populateMaterials();
        addPlanetButton_ = new QPushButton(QStringLiteral("Добавить"), this);
        deletePlanetButton_ = new QPushButton(this);
        deletePlanetButton_->setIcon(style()->standardIcon(QStyle::SP_TrashIcon));
        deletePlanetButton_->setToolTip(QStringLiteral("Удалить планету"));
        deletePlanetButton_->setVisible(false);

        auto *planetSelectorLayout = new QHBoxLayout();
        planetSelectorLayout->addWidget(planetComboBox_);
        planetSelectorLayout->addWidget(deletePlanetButton_);

        auto *planetFormLayout = new QFormLayout();
        planetFormLayout->addRow(QStringLiteral("Планета:"), planetSelectorLayout);
        planetFormLayout->addRow(QStringLiteral("Большая полуось (а.е.):"), planetSemiMajorAxisLabel_);
        planetFormLayout->addRow(QStringLiteral("Длина суток (земн. дни):"), planetDayLengthLabel_);
        planetFormLayout->addRow(QStringLiteral("Эксцентриситет орбиты:"), planetEccentricityLabel_);
        planetFormLayout->addRow(QStringLiteral("Наклон оси (°):"), planetObliquityLabel_);
        planetFormLayout->addRow(QStringLiteral("Аргумент перицентра (°):"),
                                 planetPerihelionArgumentLabel_);
        planetFormLayout->addRow(QStringLiteral("Материал поверхности:"), materialComboBox_);
        planetFormLayout->addRow(QStringLiteral("Солнечная постоянная (Вт/м²):"), resultLabel_);
        planetFormLayout->addRow(QString(), addPlanetButton_);
        auto *planetGroupBox = new QGroupBox(QStringLiteral("Планеты"), this);
        planetGroupBox->setLayout(planetFormLayout);

        auto *starsPanelLayout = new QVBoxLayout();
        starsPanelLayout->addWidget(primaryGroupBox_);
        starsPanelLayout->addWidget(secondStarCheckBox_);
        starsPanelLayout->addWidget(secondaryGroupBox_);
        auto *starsPanel = new QGroupBox(QStringLiteral("Панель звезд"), this);
        starsPanel->setLayout(starsPanelLayout);

        temperaturePlot_ = new SurfaceTemperaturePlot(this);
        auto *plotGroupBox = new QGroupBox(QStringLiteral("Температурный профиль"), this);
        auto *plotLayout = new QVBoxLayout(plotGroupBox);
        auto *segmentLayout = new QHBoxLayout();
        segmentComboBox_ = new QComboBox(plotGroupBox);
        segmentComboBox_->setEnabled(false);
        segmentLayout->addWidget(new QLabel(QStringLiteral("Сегмент орбиты:"), plotGroupBox));
        segmentLayout->addWidget(segmentComboBox_, 1);
        plotLayout->addLayout(segmentLayout);
        plotLayout->addWidget(temperaturePlot_);
        plotGroupBox->setLayout(plotLayout);

        auto *leftLayout = new QVBoxLayout();
        leftLayout->addLayout(presetsLayout);
        leftLayout->addWidget(starsPanel);
        leftLayout->addWidget(planetGroupBox);
        leftLayout->addWidget(calculateButton);
        leftLayout->addStretch();

        temperatureProgressDialog_ = new QProgressDialog(this);
        temperatureProgressDialog_->setWindowTitle(QStringLiteral("Расчет температур"));
        temperatureProgressDialog_->setLabelText(QStringLiteral("Вычисление температурного профиля..."));
        temperatureProgressDialog_->setCancelButtonText(QStringLiteral("Отмена"));
        temperatureProgressDialog_->setWindowModality(Qt::WindowModal);
        temperatureProgressDialog_->setAutoClose(true);
        temperatureProgressDialog_->setAutoReset(true);
        temperatureProgressDialog_->hide();

        connect(temperatureProgressDialog_, &QProgressDialog::canceled, this, [this]() {
            cancelTemperatureCalculation();
        });

        auto *layout = new QHBoxLayout(this);
        layout->addLayout(leftLayout, 0);
        layout->addWidget(plotGroupBox, 1);

        setLayout(layout);
        resize(480, 360);
        // Оставляем один поток под UI, чтобы параллельные вычисления не блокировали интерфейс.
        QThreadPool::globalInstance()->setMaxThreadCount(
            qMax(1, QThread::idealThreadCount() - 1));

        connect(planetComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int) {
            cancelTemperatureCalculation();
            updatePlanetSemiMajorAxisLabel();
            updatePlanetDayLengthLabel();
            updatePlanetOrbitLabels();
            syncMaterialWithPlanet();
            updatePlanetActions();
            if (hasPrimaryInputs() && (!secondStarCheckBox_->isChecked() || hasSecondaryInputs())) {
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
            updatePlanetActions();
            updateTemperaturePlot();
        });

        connect(materialComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int) {
            syncPlanetMaterialWithSelection();
            updateTemperaturePlot();
        });

        connect(segmentComboBox_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
                [this](int) { updateTemperaturePlotForSelectedSegment(); });

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
    QLabel *planetEccentricityLabel_ = nullptr;
    QLabel *planetObliquityLabel_ = nullptr;
    QLabel *planetPerihelionArgumentLabel_ = nullptr;
    QComboBox *materialComboBox_ = nullptr;
    QPushButton *addPlanetButton_ = nullptr;
    QPushButton *deletePlanetButton_ = nullptr;

    QLabel *resultLabel_ = nullptr;
    SurfaceTemperaturePlot *temperaturePlot_ = nullptr;
    QComboBox *segmentComboBox_ = nullptr;
    QProgressDialog *temperatureProgressDialog_ = nullptr;
    std::shared_ptr<std::atomic_bool> temperatureCancelFlag_;
    int temperatureRequestId_ = 0;
    int precision_ = kDefaultPrecision;
    QSet<QString> presetPlanetNames_;
    double lastSolarConstant_ = 0.0;
    bool hasSolarConstant_ = false;
    QVector<OrbitSegment> lastOrbitSegments_;
    QVector<QVector<TemperatureRangePoint>> lastTemperatureSegments_;

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
        updatePlanetOrbitLabels();
        syncMaterialWithPlanet();
        updatePlanetActions();
    }

    void clearPlanetPresets() {
        const QSignalBlocker blocker(planetComboBox_);
        planetComboBox_->clear();
        presetPlanetNames_.clear();
        planetSemiMajorAxisLabel_->setText(QStringLiteral("—"));
        planetDayLengthLabel_->setText(QStringLiteral("—"));
        planetEccentricityLabel_->setText(QStringLiteral("—"));
        planetObliquityLabel_->setText(QStringLiteral("—"));
        planetPerihelionArgumentLabel_->setText(QStringLiteral("—"));
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
        planetComboBox_->setItemData(index, isCustom, kRoleIsCustom);
        planetComboBox_->setItemData(index, planet.name, kRolePlanetName);
        planetComboBox_->setItemData(index, planet.surfaceMaterialId, kRoleMaterialId);
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

        auto *formLayout = new QFormLayout(&dialog);
        formLayout->addRow(QStringLiteral("Имя:"), nameInput);
        formLayout->addRow(QStringLiteral("Большая полуось (а.е.):"), axisInput);
        formLayout->addRow(QStringLiteral("Длина суток (земн. дни):"), dayLengthInput);
        formLayout->addRow(QStringLiteral("Эксцентриситет:"), eccentricityInput);
        formLayout->addRow(QStringLiteral("Наклон оси (°):"), obliquityInput);
        formLayout->addRow(QStringLiteral("Аргумент перицентра (°):"), perihelionArgumentInput);

        auto *materialInput = new QComboBox(&dialog);
        for (const auto &material : surfaceMaterials()) {
            materialInput->addItem(material.name, material.id);
        }
        formLayout->addRow(QStringLiteral("Материал поверхности:"), materialInput);

        auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
        formLayout->addWidget(buttons);

        connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
        connect(buttons, &QDialogButtonBox::accepted, &dialog,
                [&dialog, nameInput, axisInput, dayLengthInput, eccentricityInput,
                 obliquityInput, perihelionArgumentInput, materialInput, this]() {
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

            const int existingIndex = findPlanetIndexByName(name);
            const QString materialId = materialInput->currentData().toString();
            PlanetPreset preset{name, axis, dayLength, eccentricity, obliquity,
                                perihelionArgument, materialId};
            if (existingIndex >= 0) {
                if (!isCustomPlanetIndex(existingIndex)) {
                    showInputError(QStringLiteral("Нельзя заменить планету из пресета."));
                    return;
                }
                planetComboBox_->setItemText(existingIndex, formatPlanetName(preset));
                planetComboBox_->setItemData(existingIndex, axis, kRoleSemiMajorAxis);
                planetComboBox_->setItemData(existingIndex, dayLength, kRoleDayLength);
                planetComboBox_->setItemData(existingIndex, eccentricity, kRoleEccentricity);
                planetComboBox_->setItemData(existingIndex, obliquity, kRoleObliquity);
                planetComboBox_->setItemData(existingIndex, perihelionArgument,
                                             kRolePerihelionArgument);
                planetComboBox_->setItemData(existingIndex, true, kRoleIsCustom);
                planetComboBox_->setItemData(existingIndex, name, kRolePlanetName);
                planetComboBox_->setItemData(existingIndex, materialId, kRoleMaterialId);
                planetComboBox_->setCurrentIndex(existingIndex);
            } else {
                addPlanetItem(preset, true);
                planetComboBox_->setCurrentIndex(planetComboBox_->count() - 1);
            }

            updatePlanetSemiMajorAxisLabel();
            updatePlanetDayLengthLabel();
            updatePlanetOrbitLabels();
            syncMaterialWithPlanet();
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

    void syncPlanetMaterialWithSelection() {
        const int index = planetComboBox_->currentIndex();
        if (index < 0) {
            return;
        }
        planetComboBox_->setItemData(index, materialComboBox_->currentData(), kRoleMaterialId);
    }

    void updateSegmentComboBox() {
        const QSignalBlocker blocker(segmentComboBox_);
        const int previousIndex = segmentComboBox_->currentIndex();
        segmentComboBox_->clear();

        if (lastOrbitSegments_.isEmpty()) {
            segmentComboBox_->setEnabled(false);
            return;
        }

        for (const auto &segment : lastOrbitSegments_) {
            segmentComboBox_->addItem(formatSegmentLabel(segment));
        }

        segmentComboBox_->setEnabled(true);
        if (previousIndex >= 0 && previousIndex < segmentComboBox_->count()) {
            segmentComboBox_->setCurrentIndex(previousIndex);
        } else {
            segmentComboBox_->setCurrentIndex(0);
        }
    }

    void updateTemperaturePlotForSelectedSegment() {
        if (segmentComboBox_->currentIndex() < 0 ||
            segmentComboBox_->currentIndex() >= lastTemperatureSegments_.size()) {
            temperaturePlot_->clearSeries();
            return;
        }

        const int index = segmentComboBox_->currentIndex();
        const QString label = segmentComboBox_->itemText(index);
        temperaturePlot_->setTemperatureSeries(lastTemperatureSegments_.at(index), label);
    }

    void clearTemperatureSegments() {
        lastOrbitSegments_.clear();
        lastTemperatureSegments_.clear();
        segmentComboBox_->clear();
        segmentComboBox_->setEnabled(false);
        temperaturePlot_->clearSeries();
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

        lastTemperatureSegments_.clear();

        const int requestId = ++temperatureRequestId_;
        temperatureCancelFlag_ = std::make_shared<std::atomic_bool>(false);
        const auto cancelFlag = temperatureCancelFlag_;
        const int stepDegrees = 1;
        const int segmentCount = 12;
        const int totalLatitudes = 180 / stepDegrees + 1;
        const int totalProgress = totalLatitudes * segmentCount;

        temperatureProgressDialog_->setRange(0, totalProgress);
        temperatureProgressDialog_->setValue(0);
        temperatureProgressDialog_->show();

        QPointer<QProgressDialog> dialogGuard(temperatureProgressDialog_);
        QPointer<SolarCalculatorWidget> widgetGuard(this);
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

        SurfaceTemperatureCalculator calculator(lastSolarConstant_, *material, dayLength);
        auto *watcher = new QFutureWatcher<QVector<QVector<TemperatureRangePoint>>>(this);
        connect(watcher, &QFutureWatcher<QVector<QVector<TemperatureRangePoint>>>::finished, this,
                [this, watcher, requestId, cancelFlag]() {
            watcher->deleteLater();
            if (requestId == temperatureRequestId_) {
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
            updateSegmentComboBox();
            updateTemperaturePlotForSelectedSegment();
        });

        const OrbitSegmentCalculator orbitCalculator(semiMajorAxis, eccentricity);
        lastOrbitSegments_ = orbitCalculator.segments(segmentCount);
        updateSegmentComboBox();
        const auto segments = lastOrbitSegments_;

        auto progressCounter = std::make_shared<std::atomic_int>(0);
        auto segmentProgress = [progressCounter, progressCallback, totalProgress, cancelFlag](int, int) {
            if (cancelFlag && cancelFlag->load()) {
                return;
            }
            const int processed = progressCounter->fetch_add(1) + 1;
            progressCallback(processed, totalProgress);
        };

        auto mapSegment = [calculator,
                           stepDegrees,
                           segmentProgress,
                           cancelFlag,
                           semiMajorAxis,
                           obliquity,
                           perihelionArgument](const OrbitSegment &segment) {
            if (cancelFlag && cancelFlag->load()) {
                return QVector<TemperatureRangePoint>{};
            }
            return calculator.temperatureRangesForOrbitSegment(segment,
                                                               semiMajorAxis,
                                                               obliquity,
                                                               perihelionArgument,
                                                               stepDegrees,
                                                               segmentProgress,
                                                               cancelFlag.get());
        };

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
        if (temperatureProgressDialog_) {
            temperatureProgressDialog_->hide();
        }
    }

    void resetSolarConstant() {
        hasSolarConstant_ = false;
        lastSolarConstant_ = 0.0;
        resultLabel_->setText(
            QStringLiteral("Введите параметры и нажмите \"Рассчитать\"."));
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
