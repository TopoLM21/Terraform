#include "StarOcclusionGeometry.h"

#include <QLineF>

StarOcclusionGeometry::Result StarOcclusionGeometry::evaluate(const QPointF &screenCenter,
                                                              double sphereRadiusPx,
                                                              const QVector3D &screenDirection,
                                                              double starRadiusPx) const {
    Result result;
    result.starRadiusPx = starRadiusPx;
    result.starCenter = QPointF(screenCenter.x() + screenDirection.x() * sphereRadiusPx,
                                screenCenter.y() - screenDirection.y() * sphereRadiusPx);

    if (sphereRadiusPx <= 0.0 || starRadiusPx <= 0.0) {
        result.hidden = true;
        return result;
    }

    // Звезда позади сферы (z <= 0 в координатах экрана) — полностью скрываем.
    if (screenDirection.z() <= 0.0f) {
        result.hidden = true;
        return result;
    }

    const double centerDistance = QLineF(screenCenter, result.starCenter).length();

    // Сравниваем угловое расстояние с радиусом диска в экранных координатах:
    // расстояние < радиуса означает, что звезда находится за планетой.
    if (centerDistance <= sphereRadiusPx - starRadiusPx) {
        result.hidden = true;
        return result;
    }

    if (centerDistance >= sphereRadiusPx + starRadiusPx) {
        return result;
    }

    // Частичное перекрытие: рисуем только ту часть диска, что остаётся вне планеты.
    QPainterPath starPath;
    starPath.addEllipse(result.starCenter, starRadiusPx, starRadiusPx);
    QPainterPath planetPath;
    planetPath.addEllipse(screenCenter, sphereRadiusPx, sphereRadiusPx);
    result.clipPath = starPath.subtracted(planetPath);
    result.hasClip = true;
    return result;
}
