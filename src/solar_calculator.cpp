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

ArgumentsParseResult parseParametersFromArguments(const QCoreApplication &app,
                                                 QTextStream &output,
                                                 StellarParameters &parameters,
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
                                     QStringLiteral("Расстояние до планеты в а.е."),
                                     QStringLiteral("value"));
    QCommandLineOption precisionOption({QStringLiteral("p"), QStringLiteral("precision")},
                                       QStringLiteral("Количество значащих цифр в выводе."),
                                       QStringLiteral("digits"),
                                       QString::number(kDefaultPrecision));

    parser.addOptions({radiusOption, temperatureOption, distanceOption, precisionOption});
    parser.process(app);

    const auto parsePositive = [&](const QCommandLineOption &option, double &target,
                                   const QString &name) -> bool {
        const QString valueString = parser.value(option);
        bool ok = false;
        const double parsed = valueString.toDouble(&ok);
        if (!ok || !std::isfinite(parsed) || parsed <= 0.0) {
            output << "Значение для " << name
                   << " должно быть положительным числом." << Qt::endl;
            return false;
        }
        target = parsed;
        return true;
    };

    const bool parametersProvided = parser.isSet(radiusOption) ||
                                    parser.isSet(temperatureOption) ||
                                    parser.isSet(distanceOption);
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

    if (!parametersProvided) {
        return ArgumentsParseResult::None;
    }

    if (!parser.isSet(radiusOption) || !parser.isSet(temperatureOption) ||
        !parser.isSet(distanceOption)) {
        output << "Для неинтерактивного режима укажите все параметры: --radius,"
               << " --temperature и --distance." << Qt::endl;
        return ArgumentsParseResult::Failure;
    }

        radiusInput_ = new QLineEdit(this);
        radiusInput_->setPlaceholderText(QStringLiteral("Например, 1.0"));
        radiusInput_->setValidator(validator);

        temperatureInput_ = new QLineEdit(this);
        temperatureInput_->setPlaceholderText(QStringLiteral("Например, 5772"));
        temperatureInput_->setValidator(validator);

        distanceInput_ = new QLineEdit(this);
        distanceInput_->setPlaceholderText(QStringLiteral("Например, 1.0"));
        distanceInput_->setValidator(validator);

void promptAndComputeSolarConstant(QTextStream &input, QTextStream &output,
                                   const int precision) {
    output.setRealNumberPrecision(precision);
    output.setRealNumberNotation(QTextStream::SmartNotation);

    const auto readValue = [&](const QString &prompt, double &value) -> bool {
        output << prompt << Qt::endl;
        output.flush();
        input >> value;

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

    StellarParameters parameters{};
    int precision = kDefaultPrecision;
    const ArgumentsParseResult argsResult =
        parseParametersFromArguments(app, output, parameters, precision);

    switch (argsResult) {
    case ArgumentsParseResult::Success: {
        output.setRealNumberPrecision(precision);
        output.setRealNumberNotation(QTextStream::SmartNotation);
        const double flux = SolarCalculator::solarConstant(parameters);

        resultLabel_->setText(
            QStringLiteral("Солнечная постоянная у планеты: %1 Вт/м²").arg(flux, 0, 'g', 12));
    }

    void showInputError(const QString &message) {
        QMessageBox::warning(this, QStringLiteral("Некорректный ввод"), message);
    }

    promptAndComputeSolarConstant(input, output, precision);
    return 0;
}
