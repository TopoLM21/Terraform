#include "surface_realistic_color.h"

#include "planet_presets.h"

#include <QHash>
#include <QtCore/QString>

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
} // namespace

QColor realisticSurfaceColor(const SurfacePoint &point,
                             double minHeightKm,
                             double maxHeightKm) {
    Q_UNUSED(minHeightKm);
    Q_UNUSED(maxHeightKm);

    const QColor materialColor = baseColorForMaterialId(point.materialId);
    if (materialColor.isValid()) {
        return materialColor;
    }

    // Если материал не распознан, используем нейтральный тон,
    // чтобы реалистичный режим зависел только от типа материала, а не от высоты/температуры.
    const QColor fallbackMaterial = baseColorForMaterialId(QStringLiteral("rocky"));
    return fallbackMaterial.isValid() ? fallbackMaterial : QColor(128, 128, 128);
}
