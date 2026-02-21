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
    const double radius =
        semiMajorAxis * (1.0 - eccentricity * eccentricity) /
        (1.0 + eccentricity * qCos(trueAnomalyRadians));
    const double x = radius * qCos(trueAnomalyRadians);
    const double y = radius * qSin(trueAnomalyRadians);

    const double cosArg = qCos(perihelionArgumentRadians);
    const double sinArg = qSin(perihelionArgumentRadians);
    return QPointF(x * cosArg - y * sinArg, x * sinArg + y * cosArg);
}

QPointF toScreen(const QPointF &positionAu, const QPointF &center, double scale) {
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
    const double factor = qSqrt((1.0 + eccentricity) / (1.0 - eccentricity));
    return 2.0 * std::atan2(factor * qSin(eccentricAnomaly / 2.0),
                            qCos(eccentricAnomaly / 2.0));
}

QPointF orbitPointAuInclined(double semiMajorAxis,
                             double eccentricity,
                             double perihelionArgumentRadians,
                             double inclinationRadians,
                             double trueAnomalyRadians) {
    const double radius =
        semiMajorAxis * (1.0 - eccentricity * eccentricity) /
        (1.0 + eccentricity * qCos(trueAnomalyRadians));
    const double xOrbital = radius * qCos(trueAnomalyRadians);
    const double yOrbital = radius * qSin(trueAnomalyRadians);

    const double cosArg = qCos(perihelionArgumentRadians);
    const double sinArg = qSin(perihelionArgumentRadians);
    const double xRotated = xOrbital * cosArg - yOrbital * sinArg;
    const double yRotated = xOrbital * sinArg + yOrbital * cosArg;
    const double yInclined = yRotated * qCos(inclinationRadians);
    return QPointF(xRotated, yInclined);
}

QColor starColorForParameters(double temperatureKelvin, double radiusSolar) {
    QColor starColor = starColorFromTemperature(temperatureKelvin);
    const double brightnessFactor = qBound(0.75, 0.9 + 0.1 * qSqrt(radiusSolar), 1.25);
    starColor.setRedF(qMin(1.0, starColor.redF() * brightnessFactor));
    starColor.setGreenF(qMin(1.0, starColor.greenF() * brightnessFactor));
    starColor.setBlueF(qMin(1.0, starColor.blueF() * brightnessFactor));
    return starColor;
}

const StarSystemTopView::BinarySubsystem *findBinarySubsystem(
    const QVector<StarSystemTopView::BinarySubsystem> &subsystems,
    int hostPlanetIndex) {
    for (const StarSystemTopView::BinarySubsystem &subsystem : subsystems) {
        if (subsystem.hostPlanetIndex == hostPlanetIndex && subsystem.orbit.semiMajorAxisAu > 0.0 &&
            subsystem.orbit.periodDays > 0.0) {
            return &subsystem;
        }
    }
    return nullptr;
}

QPair<QPointF, QPointF> binaryStarOffsets(const StarSystemTopView::BinarySubsystem &subsystem) {
    const double meanAnomaly = 2.0 * M_PI * (subsystem.elapsedDays / subsystem.orbit.periodDays);
    const double trueAnomaly =
        trueAnomalyFromMeanAnomalyRadians(meanAnomaly, subsystem.orbit.eccentricity);
    const double inclinationRadians = qDegreesToRadians(subsystem.orbit.inclinationDegrees);

    const QPointF primaryOffset = orbitPointAuInclined(
        subsystem.orbit.semiMajorAxisAu,
        subsystem.orbit.eccentricity,
        qDegreesToRadians(subsystem.orbit.argumentPericenterDegreesA),
        inclinationRadians,
        trueAnomaly);
    const QPointF secondaryOffset = orbitPointAuInclined(
        subsystem.orbit.semiMajorAxisAu,
        subsystem.orbit.eccentricity,
        qDegreesToRadians(subsystem.orbit.argumentPericenterDegreesB),
        inclinationRadians,
        trueAnomaly);
    return qMakePair(primaryOffset, secondaryOffset);
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

void StarSystemTopView::setBinarySubsystems(const QVector<BinarySubsystem> &subsystems) {
    binarySubsystems_ = subsystems;
    update();
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

    // Центральная звезда системы.
    painter.setBrush(starColorForParameters(starTemperatureKelvin_, starRadiusSolar_));
    const qreal starRadiusPixels =
        qBound(kMinStarRadius, kSunRadius * qSqrt(starRadiusSolar_), kMaxStarRadius);
    painter.drawEllipse(center, starRadiusPixels, starRadiusPixels);

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

        if (i >= planetPositions.size()) {
            continue;
        }

        const QPointF planetScreenPos = planetPositions[i];
        const QPointF planetAuPos((planetScreenPos.x() - center.x()) / scale,
                                  -(planetScreenPos.y() - center.y()) / scale);

        const BinarySubsystem *subsystem = findBinarySubsystem(binarySubsystems_, i);
        if (subsystem) {
            const QPair<QPointF, QPointF> pairOffsets = binaryStarOffsets(*subsystem);
            const QPointF primaryScreen = toScreen(planetAuPos + pairOffsets.first, center, scale);
            const QPointF secondaryScreen = toScreen(planetAuPos + pairOffsets.second, center, scale);

            painter.setPen(Qt::NoPen);
            painter.setBrush(starColorForParameters(subsystem->primaryTemperatureKelvin,
                                                    qMax(0.05, subsystem->primaryRadiusSolar)));
            const qreal primaryRadiusPixels =
                qBound(kMinStarRadius - 1.0,
                       kSunRadius * qSqrt(qMax(0.05, subsystem->primaryRadiusSolar)),
                       kMaxStarRadius - 1.0);
            painter.drawEllipse(primaryScreen, primaryRadiusPixels, primaryRadiusPixels);

            painter.setBrush(starColorForParameters(subsystem->secondaryTemperatureKelvin,
                                                    qMax(0.05, subsystem->secondaryRadiusSolar)));
            const qreal secondaryRadiusPixels =
                qBound(kMinStarRadius - 1.0,
                       kSunRadius * qSqrt(qMax(0.05, subsystem->secondaryRadiusSolar)),
                       kMaxStarRadius - 1.0);
            painter.drawEllipse(secondaryScreen, secondaryRadiusPixels, secondaryRadiusPixels);
        } else {
            painter.setPen(Qt::NoPen);
            painter.setBrush(selected ? QColor(90, 190, 255) : QColor(200, 200, 200));
            const qreal pointRadius = selected ? kPlanetRadiusSelected : kPlanetRadius;
            painter.drawEllipse(planetScreenPos, pointRadius, pointRadius);
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
    for (const PlanetOrbit &planet : planets_) {
        if (planet.semiMajorAxisAu <= 0.0) {
            continue;
        }
        const double farthest = planet.semiMajorAxisAu * (1.0 + qMax(0.0, planet.eccentricity));
        maxDistance = qMax(maxDistance, farthest);
    }
    for (const BinarySubsystem &subsystem : binarySubsystems_) {
        if (subsystem.hostPlanetIndex < 0 || subsystem.hostPlanetIndex >= planets_.size()) {
            continue;
        }
        const PlanetOrbit &host = planets_[subsystem.hostPlanetIndex];
        const double hostFarthest =
            host.semiMajorAxisAu * (1.0 + qMax(0.0, host.eccentricity));
        const double binaryFarthest =
            subsystem.orbit.semiMajorAxisAu * (1.0 + qMax(0.0, subsystem.orbit.eccentricity));
        maxDistance = qMax(maxDistance, hostFarthest + binaryFarthest);
    }
    return maxDistance;
}
