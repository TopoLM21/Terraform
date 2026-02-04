#include "surface_realistic_color.h"

#include "planet_presets.h"

#include <QHash>
#include <QtGlobal>
#include <QtCore/QString>
#include <cmath>

namespace {
const QHash<QString, QColor> &materialBaseColors() {
    static QHash<QString, QColor> cached;
    if (cached.isEmpty()) {
        const auto materials = surfaceMaterials();
        cached.reserve(materials.size());
        for (const auto &material : materials) {
            if (material.baseColor.isValid()) {
                cached.insert(material.id, material.baseColor);
            }
        }
    }
    return cached;
}

QColor baseColorForMaterialId(const QString &materialId) {
    const auto &colors = materialBaseColors();
    const auto it = colors.constFind(materialId);
    return it != colors.cend() ? it.value() : QColor();
}

QColor applyHeightTint(const QColor &baseColor, double normalizedHeight) {
    const QColor lowAltitudeTint(210, 190, 170);
    const QColor highAltitudeTint(230, 235, 255);
    // Используем кривую pow с гаммой > 1, чтобы усилить различия в средних высотах.
    // Диапазон смешивания 15–35% выбран как компромисс: оттенок заметнее,
    // но материалы остаются различимыми даже на контрастных текстурах.
    const double gamma = 1.4;
    const double curvedHeight = std::pow(qBound(0.0, normalizedHeight, 1.0), gamma);
    const double centered = (curvedHeight - 0.5) * 2.0;
    const double mixFactor = 0.15 + 0.20 * std::abs(centered);
    const QColor tint = centered >= 0.0 ? highAltitudeTint : lowAltitudeTint;
    const double invMix = 1.0 - mixFactor;
    return QColor::fromRgbF(invMix * baseColor.redF() + mixFactor * tint.redF(),
                            invMix * baseColor.greenF() + mixFactor * tint.greenF(),
                            invMix * baseColor.blueF() + mixFactor * tint.blueF(),
                            baseColor.alphaF());
}

QColor applyVegetationTint(const QColor &baseColor, double vegetationBlend) {
    const QColor vegetationTint(52, 120, 68);
    // Смешиваем не более чем на 55%, чтобы материал оставался различимым.
    const double mixFactor = qBound(0.0, vegetationBlend * 0.55, 0.55);
    const double invMix = 1.0 - mixFactor;
    return QColor::fromRgbF(invMix * baseColor.redF() + mixFactor * vegetationTint.redF(),
                            invMix * baseColor.greenF() + mixFactor * vegetationTint.greenF(),
                            invMix * baseColor.blueF() + mixFactor * vegetationTint.blueF(),
                            baseColor.alphaF());
}

QColor applySnowTint(const QColor &baseColor, double snowKgPerM2) {
    if (snowKgPerM2 <= 0.0) {
        return baseColor;
    }
    const QColor snowTint(240, 245, 255);
    // Приближение: 30 кг/м² (~3 см водного эквивалента) достаточно для сплошного покрова.
    const double fullCoverKgPerM2 = 30.0;
    const double mixFactor = qBound(0.0, snowKgPerM2 / fullCoverKgPerM2, 1.0);
    const double invMix = 1.0 - mixFactor;
    return QColor::fromRgbF(invMix * baseColor.redF() + mixFactor * snowTint.redF(),
                            invMix * baseColor.greenF() + mixFactor * snowTint.greenF(),
                            invMix * baseColor.blueF() + mixFactor * snowTint.blueF(),
                            baseColor.alphaF());
}
} // namespace

QColor realisticSurfaceColor(const SurfacePoint &point,
                             double minHeightKm,
                             double maxHeightKm) {
    const double heightRangeKm = maxHeightKm - minHeightKm;
    double normalizedHeight = 0.0;
    if (heightRangeKm > 0.0) {
        normalizedHeight = qBound(0.0, (point.heightKm - minHeightKm) / heightRangeKm, 1.0);
    }

    const QString effectiveMaterialId =
        (point.materialId == QLatin1String("ocean") && point.waterPhase == PhaseModel::Phase::Ice)
            ? QStringLiteral("ice")
            : point.materialId;
    const QColor materialColor = baseColorForMaterialId(effectiveMaterialId);
    if (materialColor.isValid()) {
        const QColor heightTinted = applyHeightTint(materialColor, normalizedHeight);
        if (point.materialId != QLatin1String("ocean")) {
            const double vegetationBlend =
                qBound(0.0, point.vegetationFraction * point.vegetationBlendMask, 1.0);
            const QColor vegetated = applyVegetationTint(heightTinted, vegetationBlend);
            return applySnowTint(vegetated, point.snowKgPerM2);
        }
        return heightTinted;
    }

    // Если материал не распознан, используем нейтральный тон,
    // чтобы реалистичный режим зависел только от типа материала, а не от высоты/температуры.
    const QColor fallbackMaterial = baseColorForMaterialId(QStringLiteral("rocky"));
    const QColor fallbackBase = fallbackMaterial.isValid() ? fallbackMaterial : QColor(128, 128, 128);
    const QColor heightTinted = applyHeightTint(fallbackBase, normalizedHeight);
    if (point.materialId != QLatin1String("ocean")) {
        const double vegetationBlend =
            qBound(0.0, point.vegetationFraction * point.vegetationBlendMask, 1.0);
        const QColor vegetated = applyVegetationTint(heightTinted, vegetationBlend);
        return applySnowTint(vegetated, point.snowKgPerM2);
    }
    return heightTinted;
}
