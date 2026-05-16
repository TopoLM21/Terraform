#include "ocean_current_solver.h"
#include "fluids/PhaseModel.h"

#include <QtCore/QtMath>
#include <cmath>
#include <algorithm>

namespace {

constexpr int kNeighborCount = 6;
// Доля ветра, передающаяся океанскому течению (Экмановский перенос, упрощённо).
constexpr double kWindDragFraction = 0.03;
// Угловой сдвиг течения относительно ветра (Экман, ~45° на поверхности).
constexpr double kEkmanDeflectionRad = M_PI / 4.0;
// Горизонтальная турбулентная диффузивность тепла в океане (м²/с).
// Порядок 10³–10⁴ м²/с для крупномасштабного вихревого переноса.
constexpr double kHorizontalDiffusivityM2PerS = 2000.0;
// Экспоненциальный масштаб затухания течения с глубиной (м).
constexpr double kCurrentDepthScaleMeters = 0.5;

double angularDistanceRad(double lat1Rad, double lon1Rad,
                          double lat2Rad, double lon2Rad) {
    const double dlat = lat2Rad - lat1Rad;
    const double dlon = lon2Rad - lon1Rad;
    const double a = std::sin(dlat / 2.0) * std::sin(dlat / 2.0) +
                     std::cos(lat1Rad) * std::cos(lat2Rad) *
                         std::sin(dlon / 2.0) * std::sin(dlon / 2.0);
    return 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
}

} // namespace

void OceanCurrentSolver::ensureNeighbors(const PlanetSurfaceGrid &grid) const {
    if (cachedGrid_ == &grid && qFuzzyCompare(cachedRadiusKm_, grid.radiusKm())) {
        return;
    }
    cachedGrid_ = &grid;
    cachedRadiusKm_ = grid.radiusKm();

    const int n = grid.points().size();
    neighborIndices_.resize(n);

    // Вычислим 3D координаты на единичной сфере.
    QVector<double> cx(n), cy(n), cz(n);
    for (int i = 0; i < n; ++i) {
        const auto &p = grid.points().at(i);
        const double latRad = p.latitudeRadians;
        const double lonRad = p.longitudeRadians;
        cx[i] = std::cos(latRad) * std::cos(lonRad);
        cy[i] = std::cos(latRad) * std::sin(lonRad);
        cz[i] = std::sin(latRad);
    }

    for (int i = 0; i < n; ++i) {
        QVector<QPair<double, int>> distances;
        distances.reserve(n);
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            const double dx = cx[j] - cx[i];
            const double dy = cy[j] - cy[i];
            const double dz = cz[j] - cz[i];
            const double dist = dx * dx + dy * dy + dz * dz;
            distances.push_back({dist, j});
        }
        std::partial_sort(distances.begin(),
                          distances.begin() + qMin(kNeighborCount, distances.size()),
                          distances.end(),
                          [](const auto &a, const auto &b) { return a.first < b.first; });
        neighborIndices_[i].resize(qMin(kNeighborCount, distances.size()));
        for (int k = 0; k < neighborIndices_[i].size(); ++k) {
            neighborIndices_[i][k] = distances[k].second;
        }
    }
}

void OceanCurrentSolver::step(PlanetSurfaceGrid &grid,
                               double dayLengthSeconds,
                               bool isRetrograde,
                               double dtSeconds) {
    if (grid.points().isEmpty() || dtSeconds <= 0.0) {
        return;
    }
    ensureNeighbors(grid);

    const int n = grid.points().size();
    const double radiusMeters = grid.radiusKm() * 1000.0;
    const QLatin1String oceanId("ocean");

    // Параметр Кориолиса: f = 2Ω sin(lat).
    const double omega = (dayLengthSeconds > 0.0)
        ? 2.0 * M_PI / dayLengthSeconds
        : 0.0;
    const double omegaSigned = isRetrograde ? -omega : omega;

    // ── Фаза 1: вычисляем течения из ветра ──
    for (int i = 0; i < n; ++i) {
        auto &point = grid.points()[i];
        // Течения только для жидкого океана; замёрзший — нет течений.
        if (point.materialId != oceanId ||
            point.waterPhase == PhaseModel::Phase::Ice) {
            point.oceanCurrentEastMps = 0.0;
            point.oceanCurrentNorthMps = 0.0;
            point.oceanCurrentSpeedMps = 0.0;
            continue;
        }
        // Экмановский поток: течение = kWindDragFraction × |ветер|, повёрнут на 45°.
        const double windU = point.windEastMps;
        const double windV = point.windNorthMps;
        // Направление сдвига зависит от полушария (знак параметра Кориолиса).
        const double f = 2.0 * omegaSigned * point.sinLatitude;
        const double deflection = (f >= 0.0) ? -kEkmanDeflectionRad : kEkmanDeflectionRad;
        const double cosD = std::cos(deflection);
        const double sinD = std::sin(deflection);
        point.oceanCurrentEastMps = kWindDragFraction * (windU * cosD - windV * sinD);
        point.oceanCurrentNorthMps = kWindDragFraction * (windU * sinD + windV * cosD);
        point.oceanCurrentSpeedMps = std::hypot(point.oceanCurrentEastMps,
                                                 point.oceanCurrentNorthMps);
    }

    // ── Фаза 2: горизонтальная диффузия тепла между соседними океанскими тайлами ──
    // Для каждого слоя подповерхности: T_new = T + sum_j[ D * (T_j - T) / L² ] * dt
    // Коэффициент D уменьшается с глубиной.
    // Работаем через буфер, чтобы не зависеть от порядка обхода.

    // Собираем индексы жидких океанских точек (замёрзшие не участвуют).
    QVector<int> oceanIndices;
    oceanIndices.reserve(n);
    for (int i = 0; i < n; ++i) {
        const auto &pt = grid.points().at(i);
        if (pt.materialId == oceanId &&
            pt.waterPhase != PhaseModel::Phase::Ice) {
            oceanIndices.push_back(i);
        }
    }
    if (oceanIndices.isEmpty()) {
        return;
    }

    // Определяем число слоёв из первой океанской точки.
    const int layerCount = grid.points().at(oceanIndices.first())
                               .state.solver().temperatures().size();
    if (layerCount <= 0) {
        return;
    }

    const QVector<double> &layerDepths =
        grid.points().at(oceanIndices.first()).state.solver().layerDepthsMeters();

    // Буфер для новых температур: [oceanIndex][layer].
    QVector<QVector<double>> newTemps(oceanIndices.size());
    for (int oi = 0; oi < oceanIndices.size(); ++oi) {
        newTemps[oi] = grid.points().at(oceanIndices[oi]).state.solver().temperatures();
    }

    // Карта: глобальный индекс → позиция в oceanIndices (-1 если не океан).
    QVector<int> oceanLookup(n, -1);
    for (int oi = 0; oi < oceanIndices.size(); ++oi) {
        oceanLookup[oceanIndices[oi]] = oi;
    }

    for (int oi = 0; oi < oceanIndices.size(); ++oi) {
        const int idx = oceanIndices[oi];
        const auto &point = grid.points().at(idx);
        const auto &neighbors = (idx < neighborIndices_.size())
            ? neighborIndices_[idx]
            : point.neighborIndices;
        const QVector<double> &myTemps = point.state.solver().temperatures();

        // Считаем среднее расстояние до соседей.
        int oceanNeighborCount = 0;
        for (int ni = 0; ni < neighbors.size(); ++ni) {
            const int nIdx = neighbors[ni];
            if (nIdx < 0 || nIdx >= n) continue;
            if (oceanLookup[nIdx] < 0) continue;
            oceanNeighborCount++;
        }
        if (oceanNeighborCount == 0) continue;

        for (int layer = 0; layer < layerCount; ++layer) {
            // Глубина центра слоя.
            const double depth = (layer < layerDepths.size()) ? layerDepths[layer] : 0.0;
            // Диффузивность убывает с глубиной.
            const double depthFactor = std::exp(-depth / kCurrentDepthScaleMeters);
            const double diffusivity = kHorizontalDiffusivityM2PerS * depthFactor;

            double sumFlux = 0.0;
            for (int ni = 0; ni < neighbors.size(); ++ni) {
                const int nIdx = neighbors[ni];
                if (nIdx < 0 || nIdx >= n) continue;
                if (oceanLookup[nIdx] < 0) continue;

                const auto &neighborPoint = grid.points().at(nIdx);
                const QVector<double> &neighborTemps =
                    neighborPoint.state.solver().temperatures();
                if (layer >= neighborTemps.size()) continue;

                // Расстояние между центрами тайлов.
                const double angDist = angularDistanceRad(
                    point.latitudeRadians, point.longitudeRadians,
                    neighborPoint.latitudeRadians, neighborPoint.longitudeRadians);
                const double distanceMeters = angDist * radiusMeters;
                if (distanceMeters < 1.0) continue;

                // Диффузионный поток: D × ΔT / L²
                const double dT = neighborTemps[layer] - myTemps[layer];
                sumFlux += diffusivity * dT / (distanceMeters * distanceMeters);
            }

            // ΔT = sumFlux × dt, ограничиваем чтобы не было нестабильности.
            double deltaT = sumFlux * dtSeconds;
            deltaT = qBound(-5.0, deltaT, 5.0);
            newTemps[oi][layer] += deltaT;
            // Не позволяем уйти ниже точки замерзания для глубинных слоёв.
            if (newTemps[oi][layer] < 200.0) {
                newTemps[oi][layer] = 200.0;
            }
        }
    }

    // Применяем новые температуры.
    for (int oi = 0; oi < oceanIndices.size(); ++oi) {
        const int idx = oceanIndices[oi];
        auto &point = grid.points()[idx];
        point.state.mutableSolver().setTemperatures(newTemps[oi]);
        // Обновляем поверхностную температуру из верхнего слоя.
        point.temperatureK = point.state.solver().surfaceTemperatureKelvin();
    }
}
