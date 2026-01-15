#include "StarSystemTopView.h"

#include "star_color.h"

#include <QtCore/QtMath>
#include <QtGui/QMouseEvent>
#include <QtGui/QPainter>
#include <QtGui/QPainterPath>

#include <cmath>

namespace {
constexpr int kOrbitSteps = 240;
constexpr int kMarginPixels = 18;
constexpr qreal kSunRadius = 6.0;
constexpr qreal kPlanetRadius = 3.5;
constexpr qreal kPlanetRadiusSelected = 5.0;
constexpr qreal kOrbitWidth = 1.2;
constexpr qreal kOrbitWidthSelected = 2.2;
constexpr qreal kMinStarRadius = 4.0;
constexpr qreal kMaxStarRadius = 16.0;
constexpr int kKeplerIterations = 8;

QPointF orbitPointAu(double semiMajorAxis,
                    double eccentricity,
                    double perihelionArgumentRadians,
                    double trueAnomalyRadians) {
    // Полярное уравнение эллипса относительно фокуса (где находится звезда):
    // r = a * (1 - e^2) / (1 + e * cos(ν)), ν — истинная аномалия.
    const double radius =
        semiMajorAxis * (1.0 - eccentricity * eccentricity) /
        (1.0 + eccentricity * qCos(trueAnomalyRadians));
    const double x = radius * qCos(trueAnomalyRadians);
    const double y = radius * qSin(trueAnomalyRadians);

    // Аргумент перицентра поворачивает орбиту в плоскости экрана.
    const double cosArg = qCos(perihelionArgumentRadians);
    const double sinArg = qSin(perihelionArgumentRadians);
    return QPointF(x * cosArg - y * sinArg, x * sinArg + y * cosArg);
}

QPointF toScreen(const QPointF &positionAu, const QPointF &center, double scale) {
    // Переводим астрономические единицы в пиксели и инвертируем ось Y,
    // чтобы орбита не выглядела отражённой из-за экранных координат.
    return QPointF(center.x() + positionAu.x() * scale,
                   center.y() - positionAu.y() * scale);
}

double normalizeAngleRadians(double angle) {
    double wrapped = std::fmod(angle, 2.0 * M_PI);
    if (wrapped < 0.0) {
        wrapped += 2.0 * M_PI;
    }
    return wrapped;
}

double solveEccentricAnomalyRadians(double meanAnomaly, double eccentricity) {
    // Решаем уравнение Кеплера: M = E - e * sin(E).
    // Используем итерации Ньютона, что достаточно для малых/умеренных эксцентриситетов.
    const double normalizedM = normalizeAngleRadians(meanAnomaly);
    double eccentricAnomaly = (eccentricity < 0.8) ? normalizedM : M_PI;
    for (int i = 0; i < kKeplerIterations; ++i) {
        const double f = eccentricAnomaly - eccentricity * qSin(eccentricAnomaly) - normalizedM;
        const double fPrime = 1.0 - eccentricity * qCos(eccentricAnomaly);
        if (qFuzzyIsNull(fPrime)) {
            break;
        }
        eccentricAnomaly -= f / fPrime;
    }
    return eccentricAnomaly;
}

double trueAnomalyFromMeanAnomalyRadians(double meanAnomaly, double eccentricity) {
    const double eccentricAnomaly = solveEccentricAnomalyRadians(meanAnomaly, eccentricity);
    // Преобразование эксцентрической аномалии в истинную:
    // tan(ν/2) = sqrt((1+e)/(1-e)) * tan(E/2).
    const double factor = qSqrt((1.0 + eccentricity) / (1.0 - eccentricity));
    return 2.0 * std::atan2(factor * qSin(eccentricAnomaly / 2.0),
                            qCos(eccentricAnomaly / 2.0));
}

QPointF orbitPointAuInclined(double semiMajorAxis,
                             double eccentricity,
                             double perihelionArgumentRadians,
                             double inclinationRadians,
                             double trueAnomalyRadians) {
    // Сначала строим положение в плоскости орбиты, затем поворачиваем на аргумент
    // перицентра (ω) в плоскости и наклоняем орбиту относительно экрана.
    const double radius =
        semiMajorAxis * (1.0 - eccentricity * eccentricity) /
        (1.0 + eccentricity * qCos(trueAnomalyRadians));
    const double xOrbital = radius * qCos(trueAnomalyRadians);
    const double yOrbital = radius * qSin(trueAnomalyRadians);

    const double cosArg = qCos(perihelionArgumentRadians);
    const double sinArg = qSin(perihelionArgumentRadians);
    const double xRotated = xOrbital * cosArg - yOrbital * sinArg;
    const double yRotated = xOrbital * sinArg + yOrbital * cosArg;

    // Наклонение i задаёт поворот вокруг оси X. В проекции на экран мы видим
    // сжатие координаты Y на cos(i), а координата Z уходит из плоскости.
    const double yInclined = yRotated * qCos(inclinationRadians);
    return QPointF(xRotated, yInclined);
}

QColor starColorForParameters(double temperatureKelvin, double radiusSolar) {
    QColor starColor = starColorFromTemperature(temperatureKelvin);
    // Яркость диска слегка зависит от радиуса звезды (визуальный акцент, не фотометрия).
    const double brightnessFactor = qBound(0.75, 0.9 + 0.1 * qSqrt(radiusSolar), 1.25);
    starColor.setRedF(qMin(1.0, starColor.redF() * brightnessFactor));
    starColor.setGreenF(qMin(1.0, starColor.greenF() * brightnessFactor));
    starColor.setBlueF(qMin(1.0, starColor.blueF() * brightnessFactor));
    return starColor;
}
}

StarSystemTopView::StarSystemTopView(QWidget *parent)
    : QWidget(parent) {
    setMinimumSize(240, 240);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

void StarSystemTopView::setPlanets(const QVector<PlanetOrbit> &planets) {
    planets_ = planets;
    if (selectedIndex_ >= planets_.size()) {
        selectedIndex_ = -1;
    }
    update();
}

void StarSystemTopView::setStarParameters(double temperatureKelvin, double radiusSolar) {
    if (temperatureKelvin <= 0.0 || radiusSolar <= 0.0) {
        return;
    }
    if (qFuzzyCompare(starTemperatureKelvin_, temperatureKelvin) &&
        qFuzzyCompare(starRadiusSolar_, radiusSolar)) {
        return;
    }
    starTemperatureKelvin_ = temperatureKelvin;
    starRadiusSolar_ = radiusSolar;
    update();
}

void StarSystemTopView::setSecondaryStarParameters(double temperatureKelvin, double radiusSolar) {
    if (temperatureKelvin <= 0.0 || radiusSolar <= 0.0) {
        if (hasSecondaryStar_) {
            hasSecondaryStar_ = false;
            secondaryStarTemperatureKelvin_ = 0.0;
            secondaryStarRadiusSolar_ = 0.0;
            update();
        }
        return;
    }
    if (hasSecondaryStar_ &&
        qFuzzyCompare(secondaryStarTemperatureKelvin_, temperatureKelvin) &&
        qFuzzyCompare(secondaryStarRadiusSolar_, radiusSolar)) {
        return;
    }
    hasSecondaryStar_ = true;
    secondaryStarTemperatureKelvin_ = temperatureKelvin;
    secondaryStarRadiusSolar_ = radiusSolar;
    update();
}

void StarSystemTopView::setBinarySystemParameters(const BinaryOrbitParameters &parameters) {
    const bool enabled = parameters.semiMajorAxisAu > 0.0 && parameters.periodDays > 0.0;
    if (!enabled && !hasBinarySystem_) {
        return;
    }
    if (enabled &&
        qFuzzyCompare(binaryParameters_.semiMajorAxisAu, parameters.semiMajorAxisAu) &&
        qFuzzyCompare(binaryParameters_.periodDays, parameters.periodDays) &&
        qFuzzyCompare(binaryParameters_.eccentricity, parameters.eccentricity) &&
        qFuzzyCompare(binaryParameters_.inclinationDegrees, parameters.inclinationDegrees) &&
        qFuzzyCompare(binaryParameters_.argumentPericenterDegreesA,
                      parameters.argumentPericenterDegreesA) &&
        qFuzzyCompare(binaryParameters_.argumentPericenterDegreesB,
                      parameters.argumentPericenterDegreesB) &&
        hasBinarySystem_ == enabled) {
        return;
    }
    binaryParameters_ = parameters;
    hasBinarySystem_ = enabled;
    update();
}

void StarSystemTopView::setBinarySystemElapsedDays(double elapsedDays) {
    if (qFuzzyCompare(binaryElapsedDays_, elapsedDays)) {
        return;
    }
    binaryElapsedDays_ = elapsedDays;
    if (hasBinarySystem_) {
        update();
    }
}

void StarSystemTopView::setSelectedIndex(int index) {
    if (selectedIndex_ == index) {
        return;
    }
    selectedIndex_ = index;
    update();
}

int StarSystemTopView::selectedIndex() const {
    return selectedIndex_;
}

void StarSystemTopView::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event);
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);

    const QRectF bounds = rect();
    const QPointF center = bounds.center();
    const double maxDistanceAu = maxOrbitDistanceAu();
    const double radiusPixels = qMax(1.0, qMin(bounds.width(), bounds.height()) / 2.0 - kMarginPixels);
    const double scale = radiusPixels / maxDistanceAu;

    painter.fillRect(bounds, Qt::black);
    painter.setPen(Qt::NoPen);
    if (hasBinarySystem_ && hasSecondaryStar_) {
        // Для бинарной системы вычисляем истинную аномалию по сидерическому периоду,
        // а затем строим положение каждой звезды относительно барицентра.
        const double meanAnomaly =
            2.0 * M_PI * (binaryElapsedDays_ / binaryParameters_.periodDays);
        const double trueAnomaly =
            trueAnomalyFromMeanAnomalyRadians(meanAnomaly, binaryParameters_.eccentricity);
        const double inclinationRadians =
            qDegreesToRadians(binaryParameters_.inclinationDegrees);

        const QPointF primaryPositionAu = orbitPointAuInclined(
            binaryParameters_.semiMajorAxisAu,
            binaryParameters_.eccentricity,
            qDegreesToRadians(binaryParameters_.argumentPericenterDegreesA),
            inclinationRadians,
            trueAnomaly);
        const QPointF secondaryPositionAu = orbitPointAuInclined(
            binaryParameters_.semiMajorAxisAu,
            binaryParameters_.eccentricity,
            qDegreesToRadians(binaryParameters_.argumentPericenterDegreesB),
            inclinationRadians,
            trueAnomaly);

        const QPointF primaryScreen = toScreen(primaryPositionAu, center, scale);
        const QPointF secondaryScreen = toScreen(secondaryPositionAu, center, scale);

        painter.setBrush(starColorForParameters(starTemperatureKelvin_, starRadiusSolar_));
        const qreal primaryRadiusPixels =
            qBound(kMinStarRadius, kSunRadius * qSqrt(starRadiusSolar_), kMaxStarRadius);
        painter.drawEllipse(primaryScreen, primaryRadiusPixels, primaryRadiusPixels);

        painter.setBrush(
            starColorForParameters(secondaryStarTemperatureKelvin_, secondaryStarRadiusSolar_));
        const qreal secondaryRadiusPixels =
            qBound(kMinStarRadius, kSunRadius * qSqrt(secondaryStarRadiusSolar_), kMaxStarRadius);
        painter.drawEllipse(secondaryScreen, secondaryRadiusPixels, secondaryRadiusPixels);
    } else {
        painter.setBrush(starColorForParameters(starTemperatureKelvin_, starRadiusSolar_));
        const qreal starRadiusPixels =
            qBound(kMinStarRadius, kSunRadius * qSqrt(starRadiusSolar_), kMaxStarRadius);
        painter.drawEllipse(center, starRadiusPixels, starRadiusPixels);
    }

    const QBrush baseBrush = painter.brush();
    QVector<QPointF> planetPositions = planetScreenPositions(bounds.size().toSize());

    for (int i = 0; i < planets_.size(); ++i) {
        const PlanetOrbit &planet = planets_[i];
        if (planet.semiMajorAxisAu <= 0.0) {
            continue;
        }

        const bool selected = (i == selectedIndex_);
        const QColor orbitColor = selected ? QColor(102, 178, 255) : QColor(255, 255, 255);
        QPen orbitPen(orbitColor, selected ? kOrbitWidthSelected : kOrbitWidth);
        orbitPen.setCosmetic(true);
        painter.setPen(orbitPen);
        painter.setBrush(Qt::NoBrush);

        QPainterPath orbitPath;
        for (int step = 0; step <= kOrbitSteps; ++step) {
            const double angle = (2.0 * M_PI * step) / kOrbitSteps;
            const QPointF orbitPoint = orbitPointAu(planet.semiMajorAxisAu,
                                                    planet.eccentricity,
                                                    qDegreesToRadians(planet.perihelionArgumentDegrees),
                                                    angle);
            const QPointF screenPoint = toScreen(orbitPoint, center, scale);
            if (step == 0) {
                orbitPath.moveTo(screenPoint);
            } else {
                orbitPath.lineTo(screenPoint);
            }
        }
        painter.drawPath(orbitPath);
        painter.setBrush(baseBrush);

        if (i < planetPositions.size()) {
            painter.setPen(Qt::NoPen);
            painter.setBrush(selected ? QColor(90, 190, 255) : QColor(200, 200, 200));
            const qreal pointRadius = selected ? kPlanetRadiusSelected : kPlanetRadius;
            painter.drawEllipse(planetPositions[i], pointRadius, pointRadius);
        }
    }
}

void StarSystemTopView::mousePressEvent(QMouseEvent *event) {
    if (event->button() != Qt::LeftButton) {
        QWidget::mousePressEvent(event);
        return;
    }

    const QVector<QPointF> positions = planetScreenPositions(size());
    const qreal clickRadius = kPlanetRadiusSelected + 4.0;
    const qreal clickRadiusSquared = clickRadius * clickRadius;
    for (int i = 0; i < positions.size(); ++i) {
        const QPointF delta = positions[i] - event->pos();
        if (QPointF::dotProduct(delta, delta) <= clickRadiusSquared) {
            selectedIndex_ = i;
            update();
            emit planetSelected(i);
            return;
        }
    }

    QWidget::mousePressEvent(event);
}

QVector<QPointF> StarSystemTopView::planetScreenPositions(const QSize &size) const {
    QVector<QPointF> positions;
    positions.reserve(planets_.size());

    const QPointF center(size.width() / 2.0, size.height() / 2.0);
    const double maxDistanceAu = maxOrbitDistanceAu();
    const double radiusPixels = qMax(1.0, qMin(size.width(), size.height()) / 2.0 - kMarginPixels);
    const double scale = radiusPixels / maxDistanceAu;

    for (const PlanetOrbit &planet : planets_) {
        if (planet.semiMajorAxisAu <= 0.0) {
            positions.push_back(center);
            continue;
        }
        // Положение на орбите определяется истинной аномалией.
        const QPointF orbitPoint = orbitPointAu(planet.semiMajorAxisAu,
                                                planet.eccentricity,
                                                qDegreesToRadians(planet.perihelionArgumentDegrees),
                                                planet.trueAnomalyRadians);
        positions.push_back(toScreen(orbitPoint, center, scale));
    }

    return positions;
}

double StarSystemTopView::maxOrbitDistanceAu() const {
    double maxDistance = 1.0;
    if (hasBinarySystem_ && binaryParameters_.semiMajorAxisAu > 0.0) {
        // Апогей орбит каждой звезды вокруг барицентра.
        const double binaryFarthest =
            binaryParameters_.semiMajorAxisAu *
            (1.0 + qMax(0.0, binaryParameters_.eccentricity));
        maxDistance = qMax(maxDistance, binaryFarthest);
    }
    for (const PlanetOrbit &planet : planets_) {
        if (planet.semiMajorAxisAu <= 0.0) {
            continue;
        }
        // Апогей (a * (1 + e)) задаёт максимальное расстояние от звезды,
        // поэтому используем его для масштабирования.
        const double farthest = planet.semiMajorAxisAu * (1.0 + qMax(0.0, planet.eccentricity));
        maxDistance = qMax(maxDistance, farthest);
    }
    return maxDistance;
}
