#include "surface_point_status_dialog.h"

#include "atmospheric_column.h"

#include <QtWidgets/QFormLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtGlobal>

#include <algorithm>

namespace {
QString formatNumber(double value, int precision = 2) {
    return QString::number(value, 'f', precision);
}
} // namespace

SurfacePointStatusDialog::SurfacePointStatusDialog(QWidget *parent)
    : QDialog(parent) {
    setWindowTitle(QStringLiteral("Параметры точки поверхности"));

    auto *layout = new QFormLayout(this);

    latitudeValueLabel_ = new QLabel(QStringLiteral("—"), this);
    longitudeValueLabel_ = new QLabel(QStringLiteral("—"), this);
    temperatureValueLabel_ = new QLabel(QStringLiteral("—"), this);
    airTemperatureValueLabel_ = new QLabel(QStringLiteral("—"), this);
    heightValueLabel_ = new QLabel(QStringLiteral("—"), this);
    pressureValueLabel_ = new QLabel(QStringLiteral("—"), this);
    windValueLabel_ = new QLabel(QStringLiteral("—"), this);
    insolationValueLabel_ = new QLabel(QStringLiteral("—"), this);
    shortwaveSurfaceValueLabel_ = new QLabel(QStringLiteral("—"), this);
    longwaveUpValueLabel_ = new QLabel(QStringLiteral("—"), this);
    longwaveDownValueLabel_ = new QLabel(QStringLiteral("—"), this);
    surfaceAirFluxValueLabel_ = new QLabel(QStringLiteral("—"), this);
    subsurfaceFluxInValueLabel_ = new QLabel(QStringLiteral("—"), this);
    subsurfaceFluxOutValueLabel_ = new QLabel(QStringLiteral("—"), this);
    bottomBoundaryValueLabel_ = new QLabel(QStringLiteral("—"), this);
    bottomTemperatureValueLabel_ = new QLabel(QStringLiteral("—"), this);
    materialValueLabel_ = new QLabel(QStringLiteral("—"), this);
    tileAreaValueLabel_ = new QLabel(QStringLiteral("—"), this);
    tileEdgeLengthValueLabel_ = new QLabel(QStringLiteral("—"), this);

    profileTable_ = new QTableWidget(this);
    profileTable_->setColumnCount(4);
    profileTable_->setHorizontalHeaderLabels({
        QStringLiteral("Слой"),
        QStringLiteral("Высота (км)"),
        QStringLiteral("Температура (K)"),
        QStringLiteral("Давление (атм)")
    });
    profileTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    profileTable_->setSelectionMode(QAbstractItemView::NoSelection);
    profileTable_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    profileTable_->verticalHeader()->setVisible(false);

    subsurfaceTable_ = new QTableWidget(this);
    subsurfaceTable_->setColumnCount(3);
    subsurfaceTable_->setHorizontalHeaderLabels({
        QStringLiteral("Слой"),
        QStringLiteral("Глубина/толщина (м)"),
        QStringLiteral("Температура (K)")
    });
    subsurfaceTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    subsurfaceTable_->setSelectionMode(QAbstractItemView::NoSelection);
    subsurfaceTable_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    subsurfaceTable_->verticalHeader()->setVisible(false);

    auto *profilesWidget = new QWidget(this);
    auto *profilesLayout = new QHBoxLayout(profilesWidget);
    profilesLayout->setContentsMargins(0, 0, 0, 0);

    auto *atmosphereLayout = new QVBoxLayout();
    atmosphereLayout->setContentsMargins(0, 0, 0, 0);
    atmosphereLayout->addWidget(new QLabel(QStringLiteral("Профиль атмосферы по слоям:"), this));
    atmosphereLayout->addWidget(profileTable_);

    auto *subsurfaceLayout = new QVBoxLayout();
    subsurfaceLayout->setContentsMargins(0, 0, 0, 0);
    subsurfaceLayout->addWidget(new QLabel(QStringLiteral("Подповерхностный профиль:"), this));
    subsurfaceLayout->addWidget(subsurfaceTable_);

    profilesLayout->addLayout(atmosphereLayout);
    profilesLayout->addLayout(subsurfaceLayout);

    layout->addRow(QStringLiteral("Широта (°):"), latitudeValueLabel_);
    layout->addRow(QStringLiteral("Долгота (°):"), longitudeValueLabel_);
    layout->addRow(QStringLiteral("Температура поверхности (K):"), temperatureValueLabel_);
    layout->addRow(QStringLiteral("Температура воздуха (K):"), airTemperatureValueLabel_);
    layout->addRow(QStringLiteral("Высота (км):"), heightValueLabel_);
    layout->addRow(QStringLiteral("Поверхностное давление (атм):"), pressureValueLabel_);
    layout->addRow(QStringLiteral("Скорость ветра (м/с):"), windValueLabel_);
    layout->addRow(QStringLiteral("Солнечная инсоляция (Вт/м²):"), insolationValueLabel_);
    layout->addRow(QStringLiteral("SW у поверхности (Вт/м²):"), shortwaveSurfaceValueLabel_);
    layout->addRow(QStringLiteral("LW вверх (Вт/м²):"), longwaveUpValueLabel_);
    layout->addRow(QStringLiteral("LW вниз (Вт/м²):"), longwaveDownValueLabel_);
    layout->addRow(QStringLiteral("Поток поверхность→воздух (Вт/м²):"), surfaceAirFluxValueLabel_);
    layout->addRow(QStringLiteral("Поток в грунт (Вт/м²):"), subsurfaceFluxInValueLabel_);
    layout->addRow(QStringLiteral("Поток из грунта (Вт/м²):"), subsurfaceFluxOutValueLabel_);
    const QString bottomBoundaryTooltip = QStringLiteral(
        "Фиксированная температура задаёт «глубинный резервуар», который может "
        "подогревать профиль.");
    auto *bottomBoundaryLabel = new QLabel(QStringLiteral("Bottom boundary:"), this);
    bottomBoundaryLabel->setToolTip(bottomBoundaryTooltip);
    bottomBoundaryValueLabel_->setToolTip(bottomBoundaryTooltip);
    layout->addRow(bottomBoundaryLabel, bottomBoundaryValueLabel_);
    auto *bottomTemperatureLabel = new QLabel(QStringLiteral("Bottom temperature (K):"), this);
    bottomTemperatureLabel->setToolTip(bottomBoundaryTooltip);
    bottomTemperatureValueLabel_->setToolTip(bottomBoundaryTooltip);
    layout->addRow(bottomTemperatureLabel, bottomTemperatureValueLabel_);
    layout->addRow(QStringLiteral("Материал:"), materialValueLabel_);
    layout->addRow(QStringLiteral("Площадь тайла (км²):"), tileAreaValueLabel_);
    layout->addRow(QStringLiteral("Средняя длина ребра (км):"), tileEdgeLengthValueLabel_);
    layout->addRow(profilesWidget);
}

void SurfacePointStatusDialog::setPoint(const SurfacePoint &point,
                                        double tileAreaKm2,
                                        double tileEdgeLengthKm) {
    latitudeValueLabel_->setText(formatNumber(point.latitudeDeg));
    longitudeValueLabel_->setText(formatNumber(point.longitudeDeg));
    temperatureValueLabel_->setText(formatNumber(point.temperatureK));
    airTemperatureValueLabel_->setText(formatNumber(point.airTemperatureK));
    heightValueLabel_->setText(formatNumber(point.heightKm));
    pressureValueLabel_->setText(formatNumber(point.pressureAtm, 3));
    windValueLabel_->setText(formatNumber(point.windSpeedMps));
    insolationValueLabel_->setText(formatNumber(point.solarFluxWPerM2));
    shortwaveSurfaceValueLabel_->setText(formatNumber(point.shortwaveSurfaceWPerM2));
    longwaveUpValueLabel_->setText(formatNumber(point.longwaveUpWPerM2));
    longwaveDownValueLabel_->setText(formatNumber(point.longwaveDownWPerM2));
    surfaceAirFluxValueLabel_->setText(formatNumber(point.surfaceAirFluxWPerM2));
    // subsurfaceFluxWPerM2 положителен при потоке в грунт.
    const double fluxIntoGround = qMax(0.0, point.subsurfaceFluxWPerM2);
    const double fluxFromGround = qMax(0.0, -point.subsurfaceFluxWPerM2);
    subsurfaceFluxInValueLabel_->setText(formatNumber(fluxIntoGround));
    subsurfaceFluxOutValueLabel_->setText(formatNumber(fluxFromGround));
    materialValueLabel_->setText(point.materialId.isEmpty() ? QStringLiteral("—") : point.materialId);
    tileAreaValueLabel_->setText(formatNumber(tileAreaKm2));
    tileEdgeLengthValueLabel_->setText(formatNumber(tileEdgeLengthKm));

    const auto &solver = point.state.solver();
    const auto bottomBoundary = solver.bottomBoundary();
    const bool isFixedBoundary =
        bottomBoundary == SubsurfaceBottomBoundaryCondition::FixedTemperature;
    bottomBoundaryValueLabel_->setText(
        isFixedBoundary ? QStringLiteral("FixedTemperature") : QStringLiteral("Insulating"));
    bottomTemperatureValueLabel_->setText(
        isFixedBoundary ? formatNumber(solver.bottomTemperatureKelvin()) : QStringLiteral("—"));

    if (!subsurfaceTable_) {
        return;
    }

    // Используем solver() у SurfacePointState, потому что он хранит локальную
    // тепловую инерцию и «грунтовые» слои для конкретной точки поверхности.
    const auto &temperatures = solver.temperatures();
    const auto &layerDepths = solver.layerDepthsMeters();
    const auto &layerThicknesses = solver.layerThicknessesMeters();
    const int rows = std::min({temperatures.size(), layerDepths.size(), layerThicknesses.size()});
    if (rows == 0) {
        clearSubsurfaceProfile();
        return;
    }

    subsurfaceTable_->setRowCount(rows);
    for (int i = 0; i < rows; ++i) {
        auto *layerItem = new QTableWidgetItem(QString::number(i));
        layerItem->setFlags(layerItem->flags() & ~Qt::ItemIsEditable);
        subsurfaceTable_->setItem(i, 0, layerItem);

        const QString depthThickness = QStringLiteral("%1 / %2")
                                           .arg(formatNumber(layerDepths.at(i), 3))
                                           .arg(formatNumber(layerThicknesses.at(i), 3));
        auto *depthItem = new QTableWidgetItem(depthThickness);
        depthItem->setFlags(depthItem->flags() & ~Qt::ItemIsEditable);
        subsurfaceTable_->setItem(i, 1, depthItem);

        auto *temperatureItem = new QTableWidgetItem(formatNumber(temperatures.at(i), 2));
        temperatureItem->setFlags(temperatureItem->flags() & ~Qt::ItemIsEditable);
        subsurfaceTable_->setItem(i, 2, temperatureItem);
    }
}

void SurfacePointStatusDialog::setAtmosphereProfile(const AtmosphericColumn *column) {
    if (!profileTable_) {
        return;
    }
    if (!column || column->layers().isEmpty()) {
        clearAtmosphereProfile();
        return;
    }

    const auto &layers = column->layers();
    const QString invalidMarker = QStringLiteral("—");
    bool hasValidLayer = false;
    for (const auto &layer : layers) {
        if (layer.temperatureKelvin() > 0.0 && layer.pressureAtm() > 0.0) {
            hasValidLayer = true;
            break;
        }
    }
    if (!hasValidLayer) {
        clearAtmosphereProfile();
        return;
    }

    profileTable_->setRowCount(layers.size());
    for (int i = 0; i < layers.size(); ++i) {
        const AtmosphericLayerState &layer = layers.at(i);
        const double heightKm = layer.heightMeters() / 1000.0;
        const bool isValidLayer = layer.temperatureKelvin() > 0.0 && layer.pressureAtm() > 0.0;

        auto *layerItem = new QTableWidgetItem(QString::number(i));
        layerItem->setFlags(layerItem->flags() & ~Qt::ItemIsEditable);
        profileTable_->setItem(i, 0, layerItem);

        auto *heightItem = new QTableWidgetItem(formatNumber(heightKm, 3));
        heightItem->setFlags(heightItem->flags() & ~Qt::ItemIsEditable);
        profileTable_->setItem(i, 1, heightItem);

        auto *temperatureItem = new QTableWidgetItem(
            isValidLayer ? formatNumber(layer.temperatureKelvin(), 2) : invalidMarker);
        temperatureItem->setFlags(temperatureItem->flags() & ~Qt::ItemIsEditable);
        profileTable_->setItem(i, 2, temperatureItem);

        auto *pressureItem = new QTableWidgetItem(
            isValidLayer ? formatNumber(layer.pressureAtm(), 4) : invalidMarker);
        pressureItem->setFlags(pressureItem->flags() & ~Qt::ItemIsEditable);
        profileTable_->setItem(i, 3, pressureItem);
    }
}

void SurfacePointStatusDialog::clearPoint() {
    latitudeValueLabel_->setText(QStringLiteral("—"));
    longitudeValueLabel_->setText(QStringLiteral("—"));
    temperatureValueLabel_->setText(QStringLiteral("—"));
    airTemperatureValueLabel_->setText(QStringLiteral("—"));
    heightValueLabel_->setText(QStringLiteral("—"));
    pressureValueLabel_->setText(QStringLiteral("—"));
    windValueLabel_->setText(QStringLiteral("—"));
    insolationValueLabel_->setText(QStringLiteral("—"));
    shortwaveSurfaceValueLabel_->setText(QStringLiteral("—"));
    longwaveUpValueLabel_->setText(QStringLiteral("—"));
    longwaveDownValueLabel_->setText(QStringLiteral("—"));
    surfaceAirFluxValueLabel_->setText(QStringLiteral("—"));
    subsurfaceFluxInValueLabel_->setText(QStringLiteral("—"));
    subsurfaceFluxOutValueLabel_->setText(QStringLiteral("—"));
    bottomBoundaryValueLabel_->setText(QStringLiteral("—"));
    bottomTemperatureValueLabel_->setText(QStringLiteral("—"));
    materialValueLabel_->setText(QStringLiteral("—"));
    tileAreaValueLabel_->setText(QStringLiteral("—"));
    tileEdgeLengthValueLabel_->setText(QStringLiteral("—"));
    clearAtmosphereProfile();
    clearSubsurfaceProfile();
}

void SurfacePointStatusDialog::clearAtmosphereProfile() {
    if (!profileTable_) {
        return;
    }
    profileTable_->clearContents();
    profileTable_->setRowCount(0);
}

void SurfacePointStatusDialog::clearSubsurfaceProfile() {
    if (!subsurfaceTable_) {
        return;
    }
    subsurfaceTable_->clearContents();
    subsurfaceTable_->setRowCount(0);
}
