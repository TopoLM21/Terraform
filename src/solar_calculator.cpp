#include "solar_calculator.h"

#include <QtWidgets/QApplication>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <QtGui/QDoubleValidator>
#include <cmath>
#include <limits>

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kSolarRadiusMeters = 6.957e8;              // meters
constexpr double kAstronomicalUnitMeters = 1.496e11;         // meters
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;  // W·m−2·K−4
}

double SolarCalculator::solarConstant(const StellarParameters &parameters) {
    const double radiusMeters = parameters.radiusInSolarRadii * kSolarRadiusMeters;
    const double distanceMeters = parameters.distanceInAU * kAstronomicalUnitMeters;

    // Полная светимость звезды по закону Стефана — Больцмана: L = 4πR²σT⁴.
    const double luminosity = 4.0 * kPi * radiusMeters * radiusMeters *
                              kStefanBoltzmannConstant *
                              std::pow(parameters.temperatureKelvin, 4.0);

    // Плотность потока уменьшается пропорционально площади сферы с радиусом расстояния
    // до планеты: F = L / (4πd²).
    return luminosity / (4.0 * kPi * distanceMeters * distanceMeters);
}

class SolarCalculatorWidget : public QWidget {
public:
    explicit SolarCalculatorWidget(QWidget *parent = nullptr) : QWidget(parent) {
        setWindowTitle(QStringLiteral("Солнечная постоянная"));

        auto *validator = new QDoubleValidator(0.0001, std::numeric_limits<double>::max(), 6, this);
        validator->setNotation(QDoubleValidator::StandardNotation);

        radiusInput_ = new QLineEdit(this);
        radiusInput_->setPlaceholderText(QStringLiteral("Например, 1.0"));
        radiusInput_->setValidator(validator);

        temperatureInput_ = new QLineEdit(this);
        temperatureInput_->setPlaceholderText(QStringLiteral("Например, 5772"));
        temperatureInput_->setValidator(validator);

        distanceInput_ = new QLineEdit(this);
        distanceInput_->setPlaceholderText(QStringLiteral("Например, 1.0"));
        distanceInput_->setValidator(validator);

        auto *formLayout = new QFormLayout();
        formLayout->addRow(QStringLiteral("Радиус звезды, R☉:"), radiusInput_);
        formLayout->addRow(QStringLiteral("Температура, K:"), temperatureInput_);
        formLayout->addRow(QStringLiteral("Расстояние, а.е.:"), distanceInput_);

        auto *calculateButton = new QPushButton(QStringLiteral("Рассчитать"), this);
        connect(calculateButton, &QPushButton::clicked, this, &SolarCalculatorWidget::onCalculateRequested);

        resultLabel_ = new QLabel(QStringLiteral("Введите параметры и нажмите \"Рассчитать\"."), this);
        resultLabel_->setWordWrap(true);

        auto *layout = new QVBoxLayout(this);
        layout->addLayout(formLayout);
        layout->addWidget(calculateButton);
        layout->addWidget(resultLabel_);

        setLayout(layout);
        resize(380, 220);
    }

private:
    void onCalculateRequested() {
        bool ok = false;
        const double radius = radiusInput_->text().toDouble(&ok);
        if (!ok || radius <= 0.0) {
            showInputError(QStringLiteral("Укажите положительный радиус звезды."));
            return;
        }

        const double temperature = temperatureInput_->text().toDouble(&ok);
        if (!ok || temperature <= 0.0) {
            showInputError(QStringLiteral("Укажите положительную температуру."));
            return;
        }

        const double distance = distanceInput_->text().toDouble(&ok);
        if (!ok || distance <= 0.0) {
            showInputError(QStringLiteral("Укажите положительное расстояние."));
            return;
        }

        const StellarParameters parameters{radius, temperature, distance};
        const double flux = SolarCalculator::solarConstant(parameters);

        resultLabel_->setText(
            QStringLiteral("Солнечная постоянная у планеты: %1 Вт/м²").arg(flux, 0, 'g', 12));
    }

    void showInputError(const QString &message) {
        QMessageBox::warning(this, QStringLiteral("Некорректный ввод"), message);
    }

    QLineEdit *radiusInput_ = nullptr;
    QLineEdit *temperatureInput_ = nullptr;
    QLineEdit *distanceInput_ = nullptr;
    QLabel *resultLabel_ = nullptr;
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    SolarCalculatorWidget widget;
    widget.show();
    return app.exec();
}

