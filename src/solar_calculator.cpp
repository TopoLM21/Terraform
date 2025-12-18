#include "solar_calculator.h"

#include <QtCore/QCoreApplication>
#include <QtCore/QCommandLineOption>
#include <QtCore/QCommandLineParser>
#include <QtCore/QString>
#include <cmath>

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kSolarRadiusMeters = 6.957e8;              // meters
constexpr double kAstronomicalUnitMeters = 1.496e11;         // meters
constexpr double kStefanBoltzmannConstant = 5.670374419e-8;  // W·m−2·K−4
}

double SolarCalculator::solarConstant(const StellarParameters &parameters) {
    const double radiusMeters = parameters.radiusInSolarRadii * kSolarRadiusMeters;
    const double distanceMeters = parameters.distanceInAU * kAstronomicalUnitMeters;

    const double luminosity = 4.0 * kPi * radiusMeters * radiusMeters *
                              kStefanBoltzmannConstant *
                              std::pow(parameters.temperatureKelvin, 4.0);

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

    double radius = 0.0;
    double temperature = 0.0;
    double distance = 0.0;

    if (!parsePositive(radiusOption, radius, QStringLiteral("радиуса звезды")) ||
        !parsePositive(temperatureOption, temperature, QStringLiteral("температуры")) ||
        !parsePositive(distanceOption, distance, QStringLiteral("расстояния"))) {
        return ArgumentsParseResult::Failure;
    }

    parameters = StellarParameters{radius, temperature, distance};
    return ArgumentsParseResult::Success;
}

void promptAndComputeSolarConstant(QTextStream &input, QTextStream &output,
                                   const int precision) {
    output.setRealNumberPrecision(precision);
    output.setRealNumberNotation(QTextStream::SmartNotation);

    const auto readValue = [&](const QString &prompt, double &value) -> bool {
        output << prompt << Qt::endl;
        output.flush();
        input >> value;

        if (input.status() != QTextStream::Ok || !std::isfinite(value) || value <= 0.0) {
            output << "Некорректный ввод: требуется положительное число." << Qt::endl;
            output.flush();
            return false;
        }
        return true;
    };

    double radius = 0.0;
    if (!readValue("Введите радиус звезды в солнечных радиусах:", radius)) {
        return;
    }

    double temperature = 0.0;
    if (!readValue("Введите температуру поверхности в кельвинах:", temperature)) {
        return;
    }

    double distance = 0.0;
    if (!readValue("Введите расстояние до планеты в астрономических единицах:", distance)) {
        return;
    }

    const StellarParameters parameters{radius, temperature, distance};
    const double flux = SolarCalculator::solarConstant(parameters);

    output << "Солнечная постоянная у планеты: " << flux << " Вт/м^2" << Qt::endl;
    output.flush();
}

int main(int argc, char *argv[]) {
    QCoreApplication app(argc, argv);
    QTextStream input(stdin);
    QTextStream output(stdout);

    StellarParameters parameters{};
    int precision = kDefaultPrecision;
    const ArgumentsParseResult argsResult =
        parseParametersFromArguments(app, output, parameters, precision);

    switch (argsResult) {
    case ArgumentsParseResult::Success: {
        output.setRealNumberPrecision(precision);
        output.setRealNumberNotation(QTextStream::SmartNotation);
        const double flux = SolarCalculator::solarConstant(parameters);
        output << "Солнечная постоянная у планеты: " << flux << " Вт/м^2"
               << Qt::endl;
        output.flush();
        return 0;
    }
    case ArgumentsParseResult::Failure:
        return 1;
    case ArgumentsParseResult::None:
        break;
    }

    promptAndComputeSolarConstant(input, output, precision);
    return 0;
}
