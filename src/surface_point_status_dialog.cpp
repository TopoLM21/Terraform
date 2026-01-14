#include "surface_point_status_dialog.h"

#include <QtWidgets/QFormLayout>
#include <QtWidgets/QLabel>

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
    materialValueLabel_ = new QLabel(QStringLiteral("—"), this);
    tileAreaValueLabel_ = new QLabel(QStringLiteral("—"), this);
    tileEdgeLengthValueLabel_ = new QLabel(QStringLiteral("—"), this);

    layout->addRow(QStringLiteral("Широта (°):"), latitudeValueLabel_);
    layout->addRow(QStringLiteral("Долгота (°):"), longitudeValueLabel_);
    layout->addRow(QStringLiteral("Температура поверхности (K):"), temperatureValueLabel_);
    layout->addRow(QStringLiteral("Температура воздуха (K):"), airTemperatureValueLabel_);
    layout->addRow(QStringLiteral("Высота (км):"), heightValueLabel_);
    layout->addRow(QStringLiteral("Поверхностное давление (атм):"), pressureValueLabel_);
    layout->addRow(QStringLiteral("Скорость ветра (м/с):"), windValueLabel_);
    layout->addRow(QStringLiteral("Солнечная инсоляция (Вт/м²):"), insolationValueLabel_);
    layout->addRow(QStringLiteral("Материал:"), materialValueLabel_);
    layout->addRow(QStringLiteral("Площадь тайла (км²):"), tileAreaValueLabel_);
    layout->addRow(QStringLiteral("Средняя длина ребра (км):"), tileEdgeLengthValueLabel_);
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
    materialValueLabel_->setText(point.materialId.isEmpty() ? QStringLiteral("—") : point.materialId);
    tileAreaValueLabel_->setText(formatNumber(tileAreaKm2));
    tileEdgeLengthValueLabel_->setText(formatNumber(tileEdgeLengthKm));
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
    materialValueLabel_->setText(QStringLiteral("—"));
    tileAreaValueLabel_->setText(QStringLiteral("—"));
    tileEdgeLengthValueLabel_->setText(QStringLiteral("—"));
}
