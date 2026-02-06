#include "OceanAlbedoModel.h"

#include <QtCore/QtGlobal>

#include <cmath>

namespace {
constexpr double kWaterRefractiveIndex = 1.333;
constexpr double kMaxOceanAlbedo = 0.95;
}

double OceanAlbedoModel::liquidAlbedo(double cosZenith, double baseAlbedo) {
    const double clampedCos = qBound(0.0, cosZenith, 1.0);
    // Используем аппроксимацию Френеля (Schlick):
    // R(θ) = R0 + (1 - R0) * (1 - cosθ)^5, где R0 = ((n-1)/(n+1))^2.
    // Выбранный закон хорошо передает резкий рост отражения при касательном
    // освещении и остаётся устойчивым в диапазоне 0..1 без тяжёлых вычислений.
    const double r0 = std::pow((kWaterRefractiveIndex - 1.0) / (kWaterRefractiveIndex + 1.0), 2.0);
    const double fresnel = r0 + (1.0 - r0) * std::pow(1.0 - clampedCos, 5.0);
    // Базовое альбедо учитывает волну/пену и диффузную составляющую океана,
    // поэтому берём максимум между табличным значением и Френелем.
    // Сверху ограничиваем диапазон, чтобы избежать нереалистичной единицы
    // при сильном наклонном освещении.
    return qBound(qMax(0.0, baseAlbedo), fresnel, kMaxOceanAlbedo);
}

double OceanAlbedoModel::albedoForPhase(PhaseModel::Phase phase, double cosZenith) {
    const auto &table = PhaseModel::albedoTable();
    switch (phase) {
    case PhaseModel::Phase::Ice:
        return table.ice;
    case PhaseModel::Phase::Liquid:
        return liquidAlbedo(cosZenith, table.liquid);
    case PhaseModel::Phase::Vapor:
        return table.vapor;
    }
    return table.liquid;
}
