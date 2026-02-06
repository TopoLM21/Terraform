#pragma once

#include "PhaseModel.h"

class OceanAlbedoModel {
public:
    // Углозависимое альбедо воды с учетом приближения Френеля (Schlick).
    static double liquidAlbedo(double cosZenith, double baseAlbedo);

    // Альбедо океана: для льда используем табличное значение,
    // для жидкости — углозависимое.
    static double albedoForPhase(PhaseModel::Phase phase, double cosZenith);
};
