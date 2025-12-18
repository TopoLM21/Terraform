#include "solar_calculator.h"

#include <QtCore/QCommandLineOption>
#include <QtCore/QCommandLineParser>
#include <QtGui/QDoubleValidator>
#include <QtWidgets/QApplication>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

#include <cmath>
#include <limits>

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kSolarRadiusMeters = 6.957e8;              // meters
constexpr double kAstronomicalUnitMeters = 1.496e11;         // meters
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;  // W·m−2·K−4

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

        distanceInput_ = new QLineEdit(this);
        distanceInput_->setPlaceholderText(QStringLiteral("Например, 1.0"));
        distanceInput_->setValidator(validator);

        auto *formLayout = new QFormLayout();
        formLayout->addRow(QStringLiteral("Радиус звезды (в R☉):"), radiusInput_);
        formLayout->addRow(QStringLiteral("Температура поверхности (K):"), temperatureInput_);
        formLayout->addRow(QStringLiteral("Расстояние до планеты (а.е.):"), distanceInput_);

        auto *calculateButton = new QPushButton(QStringLiteral("Рассчитать"), this);
        connect(calculateButton, &QPushButton::clicked, this, [this]() { onCalculateRequested(); });

        resultLabel_ = new QLabel(
            QStringLiteral("Введите параметры и нажмите \"Рассчитать\"."), this);
        resultLabel_->setWordWrap(true);

        auto *layout = new QVBoxLayout(this);
        layout->addLayout(formLayout);
        layout->addWidget(calculateButton);
        layout->addWidget(resultLabel_);

        setLayout(layout);
        resize(420, 240);
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
            showInputError(QStringLiteral("Укажите положительное расстояние до планеты."));
            return;
        }

        StellarParameters parameters{radius, temperature, distance};
        const double flux = SolarCalculator::solarConstant(parameters);

        resultLabel_->setText(
            QStringLiteral("Солнечная постоянная у планеты: %1 Вт/м²")
                .arg(flux, 0, 'g', precision_));
    }

    void showInputError(const QString &message) {
        QMessageBox::warning(this, QStringLiteral("Некорректный ввод"), message);
    }

    QLineEdit *radiusInput_ = nullptr;
    QLineEdit *temperatureInput_ = nullptr;
    QLineEdit *distanceInput_ = nullptr;
    QLabel *resultLabel_ = nullptr;
    int precision_ = kDefaultPrecision;
};
}  // namespace

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

    const bool parametersProvided = parser.isSet(radiusOption) || parser.isSet(temperatureOption) ||
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

    if (!parser.isSet(radiusOption) || !parser.isSet(temperatureOption) || !parser.isSet(distanceOption)) {
        output << "Для неинтерактивного режима укажите все параметры: --radius, --temperature и --distance."
               << Qt::endl;
        return ArgumentsParseResult::Failure;
    }

    if (!parsePositive(radiusOption, parameters.radiusInSolarRadii, QStringLiteral("radius")) ||
        !parsePositive(temperatureOption, parameters.temperatureKelvin, QStringLiteral("temperature")) ||
        !parsePositive(distanceOption, parameters.distanceInAU, QStringLiteral("distance"))) {
        return ArgumentsParseResult::Failure;
    }

    return ArgumentsParseResult::Success;
}

void promptAndComputeSolarConstant(QTextStream &input, QTextStream &output, int precision) {
    output.setRealNumberPrecision(precision);
    output.setRealNumberNotation(QTextStream::SmartNotation);

    StellarParameters parameters{};

    output << "Введите радиус звезды (в солнечных радиусах):" << Qt::endl;
    output.flush();
    input >> parameters.radiusInSolarRadii;
    if (input.status() != QTextStream::Ok || parameters.radiusInSolarRadii <= 0.0) {
        output << "Радиус должен быть положительным числом." << Qt::endl;
        return;
    }

    output << "Введите температуру поверхности (в К):" << Qt::endl;
    output.flush();
    input >> parameters.temperatureKelvin;
    if (input.status() != QTextStream::Ok || parameters.temperatureKelvin <= 0.0) {
        output << "Температура должна быть положительным числом." << Qt::endl;
        return;
    }

    output << "Введите расстояние до планеты (в а.е.):" << Qt::endl;
    output.flush();
    input >> parameters.distanceInAU;
    if (input.status() != QTextStream::Ok || parameters.distanceInAU <= 0.0) {
        output << "Расстояние должно быть положительным числом." << Qt::endl;
        return;
    }

    const double flux = SolarCalculator::solarConstant(parameters);
    output << "Солнечная постоянная у планеты: " << flux << " Вт/м²" << Qt::endl;
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QTextStream output(stdout);
    QTextStream input(stdin);

    StellarParameters parameters{};
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
