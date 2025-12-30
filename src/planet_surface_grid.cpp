#include "planet_surface_grid.h"
#include "surface_height_model.h"

#include <algorithm>

#include <QHash>
#include <QVector3D>
#include <QtMath>

namespace {
constexpr double kPi = 3.14159265358979323846;
const double kGoldenAngle = kPi * (3.0 - qSqrt(5.0));

double radiansToDegrees(double radians) {
    return radians * 180.0 / kPi;
}

double degreesToRadians(double degrees) {
    return degrees * kPi / 180.0;
}

struct GridVertex {
    QVector3D position;
    QVector<int> adjacentCells;
};

QVector3D normalized(const QVector3D &v) {
    const float length = v.length();
    if (qFuzzyIsNull(length)) {
        return QVector3D();
    }
    return v / length;
}

QVector3D latLonToCartesian(double latitudeDeg, double longitudeDeg) {
    const double latRad = degreesToRadians(latitudeDeg);
    const double lonRad = degreesToRadians(longitudeDeg);
    const double cosLat = qCos(latRad);
    return QVector3D(static_cast<float>(cosLat * qSin(lonRad)),
                     static_cast<float>(qSin(latRad)),
                     static_cast<float>(cosLat * qCos(lonRad)));
}

SurfaceVertex cartesianToLatLon(const QVector3D &v) {
    const QVector3D n = normalized(v);
    const double lat = qAsin(qBound(-1.0, static_cast<double>(n.y()), 1.0));
    const double lon = qAtan2(static_cast<double>(n.x()), static_cast<double>(n.z()));
    return SurfaceVertex{radiansToDegrees(lat), radiansToDegrees(lon)};
}

struct Face {
    int a = 0;
    int b = 0;
    int c = 0;
};

int midpointIndex(int i, int j, QVector<QVector3D> *vertices, QHash<quint64, int> *cache) {
    const int minIndex = qMin(i, j);
    const int maxIndex = qMax(i, j);
    const quint64 key = (static_cast<quint64>(minIndex) << 32) | static_cast<quint64>(maxIndex);
    if (cache->contains(key)) {
        return cache->value(key);
    }
    const QVector3D midpoint = normalized(vertices->at(minIndex) + vertices->at(maxIndex));
    const int index = vertices->size();
    vertices->push_back(midpoint);
    cache->insert(key, index);
    return index;
}

double signedAngleOnPlane(const QVector3D &base,
                          const QVector3D &tangent,
                          const QVector3D &v) {
    const QVector3D vPlane = v - base * QVector3D::dotProduct(base, v);
    const QVector3D tPlane = tangent - base * QVector3D::dotProduct(base, tangent);
    const QVector3D cross = QVector3D::crossProduct(tPlane, vPlane);
    const double sinValue = QVector3D::dotProduct(base, cross);
    const double cosValue = QVector3D::dotProduct(tPlane, vPlane);
    return qAtan2(sinValue, cosValue);
}
} // namespace

void PlanetSurfaceGrid::setRadiusKm(double radiusKm) {
    radiusKm_ = radiusKm;
    if (points_.isEmpty()) {
        return;
    }
    const double surfaceAreaKm2 = 4.0 * kPi * radiusKm_ * radiusKm_;
    pointAreaKm2_ = surfaceAreaKm2 / static_cast<double>(points_.size());
    for (auto &point : points_) {
        point.radiusKm = radiusKm_;
    }
}

double PlanetSurfaceGrid::radiusKm() const {
    return radiusKm_;
}

void PlanetSurfaceGrid::generateFibonacciPoints(int pointCount) {
    points_.clear();
    cells_.clear();
    if (pointCount <= 0 || radiusKm_ <= 0.0) {
        pointAreaKm2_ = 0.0;
        return;
    }

    points_.reserve(pointCount);
    for (int i = 0; i < pointCount; ++i) {
        // Фибоначчиева сфера: равномерное распределение по площади.
        // t задаёт равномерный шаг по sin(lat), а golden angle избегает слипания точек.
        const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(pointCount);
        const double latitude = qAsin(1.0 - 2.0 * t);
        const double longitude = qAtan2(qSin(kGoldenAngle * i), qCos(kGoldenAngle * i));

        SurfacePoint point;
        point.latitudeDeg = radiansToDegrees(latitude);
        point.longitudeDeg = radiansToDegrees(longitude);
        point.radiusKm = radiusKm_;
        points_.push_back(point);
    }

    const double surfaceAreaKm2 = 4.0 * kPi * radiusKm_ * radiusKm_;
    pointAreaKm2_ = surfaceAreaKm2 / static_cast<double>(points_.size());
    applyHeightModel();
}

void PlanetSurfaceGrid::generateIcosahedronGrid(int subdivisionLevel) {
    points_.clear();
    cells_.clear();
    if (radiusKm_ <= 0.0 || subdivisionLevel < 0) {
        pointAreaKm2_ = 0.0;
        return;
    }

    rebuildIcosahedronCells(subdivisionLevel);
    if (points_.isEmpty()) {
        pointAreaKm2_ = 0.0;
        return;
    }

    const double surfaceAreaKm2 = 4.0 * kPi * radiusKm_ * radiusKm_;
    pointAreaKm2_ = surfaceAreaKm2 / static_cast<double>(points_.size());
    applyHeightModel();
}

int PlanetSurfaceGrid::pointCount() const {
    return points_.size();
}

double PlanetSurfaceGrid::pointAreaKm2() const {
    return pointAreaKm2_;
}

const QVector<SurfacePoint> &PlanetSurfaceGrid::points() const {
    return points_;
}

QVector<SurfacePoint> &PlanetSurfaceGrid::points() {
    return points_;
}

const SurfacePoint *PlanetSurfaceGrid::pointAt(int index) const {
    if (index < 0 || index >= points_.size()) {
        return nullptr;
    }
    return &points_[index];
}

const QVector<SurfaceCell> &PlanetSurfaceGrid::cells() const {
    return cells_;
}

const SurfaceCell *PlanetSurfaceGrid::cellAt(int index) const {
    if (index < 0 || index >= cells_.size()) {
        return nullptr;
    }
    return &cells_[index];
}

void PlanetSurfaceGrid::setHeightSource(HeightSourceType sourceType,
                                        const QString &heightmapPath,
                                        double heightmapScaleKm,
                                        quint32 heightSeed,
                                        bool useContinentsHeight) {
    heightSourceType_ = sourceType;
    heightmapPath_ = heightmapPath;
    heightmapScaleKm_ = heightmapScaleKm;
    heightSeed_ = heightSeed;
    useContinentsHeight_ = useContinentsHeight;
}

void PlanetSurfaceGrid::rebuildIcosahedronCells(int subdivisionLevel) {
    QVector<QVector3D> vertices;
    QVector<Face> faces;

    const double phi = (1.0 + qSqrt(5.0)) * 0.5;
    vertices.reserve(12);
    vertices.push_back(normalized(QVector3D(-1.0f, static_cast<float>(phi), 0.0f)));
    vertices.push_back(normalized(QVector3D(1.0f, static_cast<float>(phi), 0.0f)));
    vertices.push_back(normalized(QVector3D(-1.0f, static_cast<float>(-phi), 0.0f)));
    vertices.push_back(normalized(QVector3D(1.0f, static_cast<float>(-phi), 0.0f)));
    vertices.push_back(normalized(QVector3D(0.0f, -1.0f, static_cast<float>(phi))));
    vertices.push_back(normalized(QVector3D(0.0f, 1.0f, static_cast<float>(phi))));
    vertices.push_back(normalized(QVector3D(0.0f, -1.0f, static_cast<float>(-phi))));
    vertices.push_back(normalized(QVector3D(0.0f, 1.0f, static_cast<float>(-phi))));
    vertices.push_back(normalized(QVector3D(static_cast<float>(phi), 0.0f, -1.0f)));
    vertices.push_back(normalized(QVector3D(static_cast<float>(phi), 0.0f, 1.0f)));
    vertices.push_back(normalized(QVector3D(static_cast<float>(-phi), 0.0f, -1.0f)));
    vertices.push_back(normalized(QVector3D(static_cast<float>(-phi), 0.0f, 1.0f)));

    faces = {
        {0, 11, 5},
        {0, 5, 1},
        {0, 1, 7},
        {0, 7, 10},
        {0, 10, 11},
        {1, 5, 9},
        {5, 11, 4},
        {11, 10, 2},
        {10, 7, 6},
        {7, 1, 8},
        {3, 9, 4},
        {3, 4, 2},
        {3, 2, 6},
        {3, 6, 8},
        {3, 8, 9},
        {4, 9, 5},
        {2, 4, 11},
        {6, 2, 10},
        {8, 6, 7},
        {9, 8, 1},
    };

    for (int level = 0; level < subdivisionLevel; ++level) {
        QVector<Face> refined;
        refined.reserve(faces.size() * 4);
        QHash<quint64, int> midpointCache;
        midpointCache.reserve(faces.size() * 3);

        for (const Face &face : faces) {
            const int ab = midpointIndex(face.a, face.b, &vertices, &midpointCache);
            const int bc = midpointIndex(face.b, face.c, &vertices, &midpointCache);
            const int ca = midpointIndex(face.c, face.a, &vertices, &midpointCache);
            refined.push_back({face.a, ab, ca});
            refined.push_back({face.b, bc, ab});
            refined.push_back({face.c, ca, bc});
            refined.push_back({ab, bc, ca});
        }
        faces.swap(refined);
    }

    points_.reserve(faces.size());
    cells_.reserve(faces.size());

    QVector<GridVertex> gridVertices;
    gridVertices.resize(vertices.size());
    for (int i = 0; i < vertices.size(); ++i) {
        gridVertices[i].position = vertices[i];
    }

    for (int faceIndex = 0; faceIndex < faces.size(); ++faceIndex) {
        const Face &face = faces[faceIndex];
        const QVector3D centroid =
            normalized(vertices[face.a] + vertices[face.b] + vertices[face.c]);
        const SurfaceVertex latLon = cartesianToLatLon(centroid);

        SurfacePoint point;
        point.latitudeDeg = latLon.latitudeDeg;
        point.longitudeDeg = latLon.longitudeDeg;
        point.radiusKm = radiusKm_;
        points_.push_back(point);

        gridVertices[face.a].adjacentCells.push_back(faceIndex);
        gridVertices[face.b].adjacentCells.push_back(faceIndex);
        gridVertices[face.c].adjacentCells.push_back(faceIndex);
    }

    for (int faceIndex = 0; faceIndex < faces.size(); ++faceIndex) {
        const Face &face = faces[faceIndex];
        const QVector3D centroid =
            normalized(vertices[face.a] + vertices[face.b] + vertices[face.c]);
        const QVector3D radial = centroid;
        QVector3D tangent = normalized(QVector3D::crossProduct(radial, QVector3D(0.0f, 1.0f, 0.0f)));
        if (tangent.lengthSquared() < 1e-6f) {
            // Если нормаль почти совпадает с осью Y, берём запасной базис,
            // чтобы определить направление обхода вершин на касательной плоскости.
            tangent = normalized(QVector3D::crossProduct(radial, QVector3D(1.0f, 0.0f, 0.0f)));
        }
        QVector<int> cornerIndices = {face.a, face.b, face.c};

        QVector<SurfaceVertex> polygon;
        polygon.reserve(3);
        QVector<QPair<double, SurfaceVertex>> ordered;
        ordered.reserve(3);

        // Для треугольной геодезической сетки каждая ячейка - сферический треугольник
        // вокруг центра грани. Вершины полигона - исходные вершины icosahedron,
        // упорядоченные вокруг нормали грани.
        for (int vertexIndex : cornerIndices) {
            const SurfaceVertex latLon = cartesianToLatLon(vertices[vertexIndex]);
            const double angle = signedAngleOnPlane(radial, tangent, vertices[vertexIndex]);
            ordered.push_back({angle, latLon});
        }
        std::sort(ordered.begin(), ordered.end(), [](const auto &a, const auto &b) {
            return a.first < b.first;
        });
        for (const auto &entry : ordered) {
            polygon.push_back(entry.second);
        }

        SurfaceCell cell;
        cell.pointIndex = faceIndex;
        cell.polygon = polygon;
        cells_.push_back(cell);
    }
}

void PlanetSurfaceGrid::applyHeightModel() {
    SurfaceHeightModel heightModel(heightSourceType_, heightmapPath_, heightmapScaleKm_,
                                   heightSeed_, useContinentsHeight_);
    for (auto &point : points_) {
        point.heightKm = heightModel.heightKmAt(point.latitudeDeg, point.longitudeDeg);
    }
}
