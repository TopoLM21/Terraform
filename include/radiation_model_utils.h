#pragma once

#include "radiation_model.h"

// Преобразует оптическую толщину в коэффициент прохождения длинноволнового излучения
// с учетом выбранной радиационной модели.
double longwaveTransmissionForOpticalDepth(double opticalDepth, RadiationModelType type);

// Переводит пользовательскую "парниковую непрозрачность" (долю поглощения LW)
// в эквивалентную оптическую толщину для выбранного режима модели.
double opticalDepthFromGreenhouseOpacity(double greenhouseOpacity, RadiationModelType type);
