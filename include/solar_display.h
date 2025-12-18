#pragma once

#include "solar_calculator.h"

#include <QtCore/QCoreApplication>
#include <QtCore/QTextStream>

constexpr int kDefaultPrecision = 6;

enum class ArgumentsParseResult {
    None,
    Success,
    Failure,
};

ArgumentsParseResult parseParametersFromArguments(const QCoreApplication &app,
                                                  QTextStream &output,
                                                  BinarySystemParameters &parameters,
                                                  int &precision);
void promptAndComputeSolarConstant(QTextStream &input, QTextStream &output,
                                   int precision);
