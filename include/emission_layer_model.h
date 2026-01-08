#pragma once

// Модель эффективной температуры длинноволнового эмиссионного слоя.
// В серой атмосфере температура связана с оптической толщиной как
// T(τ)^4 = (3/4) * T_eff^4 * (τ + 2/3). Слой с τ≈1 считается «видимым космосу».

double opticalDepthFromGreenhouseOpacity(double greenhouseOpacity);

double emissionLayerTemperature(double surfaceTemperatureKelvin, double opticalDepth);
