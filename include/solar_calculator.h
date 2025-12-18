#pragma once

#include <QtCore/QTextStream>
#include <QtCore/QCoreApplication>

constexpr int kDefaultPrecision = 6;

struct StellarParameters {
    double radiusInSolarRadii;
    double temperatureKelvin;
    double distanceInAU;
};

class SolarCalculator {
public:
    static double solarConstant(const StellarParameters &parameters);
};

enum class ArgumentsParseResult {
    None,
    Success,
    Failure,
};

ArgumentsParseResult parseParametersFromArguments(const QCoreApplication &app,
                                                  QTextStream &output,
                                                  StellarParameters &parameters,
                                                  int &precision);
void promptAndComputeSolarConstant(QTextStream &input, QTextStream &output,
                                   int precision);
