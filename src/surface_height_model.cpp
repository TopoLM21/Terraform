#include "surface_height_model.h"

#include "planet_presets.h"

#include <QVector3D>
#include <QtMath>

#include <limits>

namespace {
constexpr double kPi = 3.14159265358979323846;

double degreesToRadians(double degrees) {
    return degrees * kPi / 180.0;
}

QVector3D unitSpherePoint(double latitudeDeg, double longitudeDeg) {
    const double latRad = degreesToRadians(latitudeDeg);
    const double lonRad = degreesToRadians(longitudeDeg);
    const double cosLat = qCos(latRad);
    return QVector3D(static_cast<float>(cosLat * qCos(lonRad)),
                     static_cast<float>(qSin(latRad)),
                     static_cast<float>(cosLat * qSin(lonRad)));
}

quint32 hash3(int x, int y, int z) {
    quint32 h = 2166136261u;
    h = (h ^ static_cast<quint32>(x)) * 16777619u;
    h = (h ^ static_cast<quint32>(y)) * 16777619u;
    h = (h ^ static_cast<quint32>(z)) * 16777619u;
    h ^= h >> 13;
    h *= 1274126177u;
    return h;
}

double randomValue(int x, int y, int z) {
    const quint32 h = hash3(x, y, z);
    return (static_cast<double>(h) / static_cast<double>(std::numeric_limits<quint32>::max())) * 2.0 - 1.0;
}

double smoothstep(double t) {
    return t * t * (3.0 - 2.0 * t);
}

double lerp(double a, double b, double t) {
    return a + (b - a) * t;
}

double valueNoise(const QVector3D &p) {
    const int ix = static_cast<int>(qFloor(p.x()));
    const int iy = static_cast<int>(qFloor(p.y()));
    const int iz = static_cast<int>(qFloor(p.z()));
    const double fx = smoothstep(p.x() - ix);
    const double fy = smoothstep(p.y() - iy);
    const double fz = smoothstep(p.z() - iz);

    const double v000 = randomValue(ix, iy, iz);
    const double v100 = randomValue(ix + 1, iy, iz);
    const double v010 = randomValue(ix, iy + 1, iz);
    const double v110 = randomValue(ix + 1, iy + 1, iz);
    const double v001 = randomValue(ix, iy, iz + 1);
    const double v101 = randomValue(ix + 1, iy, iz + 1);
    const double v011 = randomValue(ix, iy + 1, iz + 1);
    const double v111 = randomValue(ix + 1, iy + 1, iz + 1);

    const double x00 = lerp(v000, v100, fx);
    const double x10 = lerp(v010, v110, fx);
    const double x01 = lerp(v001, v101, fx);
    const double x11 = lerp(v011, v111, fx);
    const double y0 = lerp(x00, x10, fy);
    const double y1 = lerp(x01, x11, fy);
    return lerp(y0, y1, fz);
}

double fbmNoise(const QVector3D &p, double baseFrequency, int octaves) {
    double sum = 0.0;
    double amplitude = 1.0;
    double frequency = baseFrequency;
    double amplitudeSum = 0.0;
    for (int i = 0; i < octaves; ++i) {
        sum += valueNoise(p * static_cast<float>(frequency)) * amplitude;
        amplitudeSum += amplitude;
        amplitude *= 0.5;
        frequency *= 2.0;
    }
    return (amplitudeSum > 0.0) ? sum / amplitudeSum : 0.0;
}

double ridgedFbmNoise(const QVector3D &p, double baseFrequency, int octaves) {
    double sum = 0.0;
    double amplitude = 1.0;
    double frequency = baseFrequency;
    double amplitudeSum = 0.0;
    for (int i = 0; i < octaves; ++i) {
        const double n = valueNoise(p * static_cast<float>(frequency));
        const double ridge = 1.0 - qAbs(n);
        sum += ridge * ridge * amplitude;
        amplitudeSum += amplitude;
        amplitude *= 0.5;
        frequency *= 2.1;
    }
    return (amplitudeSum > 0.0) ? sum / amplitudeSum : 0.0;
}
} // namespace

SurfaceHeightModel::SurfaceHeightModel() = default;

SurfaceHeightModel::SurfaceHeightModel(HeightSourceType sourceType,
                                       const QString &heightmapPath,
                                       double heightmapScaleKm)
    : sourceType_(sourceType) {
    if (sourceType_ == HeightSourceType::HeightmapEquirectangular) {
        if (!heightmap_.loadFromFile(heightmapPath, heightmapScaleKm)) {
            sourceType_ = HeightSourceType::Procedural;
        }
    }
}

double SurfaceHeightModel::heightKmAt(double latitudeDeg, double longitudeDeg) const {
    if (sourceType_ == HeightSourceType::HeightmapEquirectangular && heightmap_.isValid()) {
        return heightmap_.heightKmAt(latitudeDeg, longitudeDeg);
    }

    const QVector3D p = unitSpherePoint(latitudeDeg, longitudeDeg);

    const double continentA = ridgedFbmNoise(p, 0.7, 4);
    const double continentB = ridgedFbmNoise(p, 1.3, 3);
    const double continentC = 0.5 + 0.5 * fbmNoise(p, 0.9, 5);
    const double continents = 0.55 * continentA + 0.30 * continentB + 0.15 * continentC;

    const double detail = fbmNoise(p, 6.5, 5);

    // Нормируем суммарный сигнал: крупный рельеф доминирует, детализация лишь подчёркивает форму.
    const double normalized = qBound(-1.0, (continents * 2.0 - 1.0) + detail * 0.15, 1.0);
    // Масштабируем в километры, чтобы получить реалистичный диапазон высот (-9..+9 км).
    return normalized * 9.0;
}
