#include "radiation_model_utils.h"

#include <QtGlobal>

#include <cmath>

double longwaveTransmissionForOpticalDepth(double opticalDepth, RadiationModelType type) {
    if (opticalDepth <= 0.0) {
        return 1.0;
    }
    if (type == RadiationModelType::Fast) {
        // Закон Бугера-Ламберта для одного слоя.
        return std::exp(-opticalDepth);
    }
    // Двухпотоковая аппроксимация для серой атмосферы.
    return 1.0 / (1.0 + 0.75 * opticalDepth);
}

double opticalDepthFromGreenhouseOpacity(double greenhouseOpacity, RadiationModelType type) {
    const double opacity = qBound(0.0, greenhouseOpacity, 0.999);
    if (opacity <= 0.0) {
        return 0.0;
    }
    const double transmission = qMax(1e-6, 1.0 - opacity);
    if (type == RadiationModelType::Fast) {
        // Для закона Бугера-Ламберта: T = exp(-tau).
        return -std::log(transmission);
    }
    // Для Eddington-приближения: T = 1 / (1 + 0.75 * tau).
    return (1.0 / transmission - 1.0) / 0.75;
}
