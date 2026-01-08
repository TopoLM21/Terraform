#include "surface_realistic_color.h"

#include "height_color_scale.h"

#include <QtCore/QString>
#include <QtMath>

namespace {
QColor blendColors(const QColor &from, const QColor &to, double ratio) {
    const double t = qBound(0.0, ratio, 1.0);
    return QColor(qRound(from.red() + (to.red() - from.red()) * t),
                  qRound(from.green() + (to.green() - from.green()) * t),
                  qRound(from.blue() + (to.blue() - from.blue()) * t));
}

bool isOceanMaterial(const QString &materialId) {
    return materialId == QLatin1String("ocean");
}

QColor baseLandColorForMaterial(const QString &materialId) {
    if (materialId == QLatin1String("forest")) {
        return QColor(52, 112, 66);
    }
    if (materialId == QLatin1String("desert")) {
        return QColor(210, 182, 120);
    }
    if (materialId == QLatin1String("ice")) {
        return QColor(220, 235, 245);
    }
    if (materialId == QLatin1String("rocky")) {
        return QColor(150, 120, 90);
    }
    if (materialId == QLatin1String("metal")) {
        return QColor(140, 140, 150);
    }
    if (materialId.startsWith(QLatin1String("regolith"))) {
        return QColor(135, 125, 115);
    }
    return QColor();
}

QColor oceanColorForDepth(double heightKm, double minHeightKm) {
    const double maxDepthKm = qMax(1.0, -minHeightKm);
    const double depthRatio = qBound(0.0, -heightKm / maxDepthKm, 1.0);
    const QColor shallow(40, 140, 200);
    const QColor deep(8, 28, 70);
    return blendColors(shallow, deep, depthRatio);
}
} // namespace

QColor realisticSurfaceColor(const SurfacePoint &point,
                             double minHeightKm,
                             double maxHeightKm) {
    const double heightKm = point.heightKm;
    const bool treatAsOcean = isOceanMaterial(point.materialId) ||
        (heightKm <= 0.0 && maxHeightKm > 0.0);

    QColor baseColor;
    if (treatAsOcean) {
        baseColor = oceanColorForDepth(heightKm, minHeightKm);
    } else {
        baseColor = baseLandColorForMaterial(point.materialId);
        if (!baseColor.isValid()) {
            const double heightRatio =
                (maxHeightKm > 0.0) ? qBound(0.0, heightKm / maxHeightKm, 1.0) : 0.5;
            baseColor = heightColorForRatio(heightRatio);
        }

        if (heightKm > 0.0 && maxHeightKm > 0.0) {
            // Подсветка вершин, чтобы высокогорные области читались реалистичнее.
            const double elevation = qBound(0.0, heightKm / maxHeightKm, 1.0);
            baseColor = blendColors(baseColor, QColor(235, 235, 235), elevation * 0.5);
        }
    }

    const double freezePointK = 273.15;
    const double iceBlendRangeK = 15.0;
    if (point.temperatureK < freezePointK) {
        // Плавный переход к снегу/льду по мере снижения температуры.
        const double iceRatio =
            qBound(0.0, (freezePointK - point.temperatureK) / iceBlendRangeK, 1.0);
        baseColor = blendColors(baseColor, QColor(220, 235, 245), iceRatio);
    }

    return baseColor;
}
