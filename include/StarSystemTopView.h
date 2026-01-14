#pragma once

#include <QtCore/QPointF>
#include <QtCore/QVector>
#include <QtCore/QString>
#include <QtWidgets/QWidget>

class StarSystemTopView : public QWidget {
    Q_OBJECT

public:
    struct PlanetOrbit {
        QString name;
        double semiMajorAxisAu = 0.0;
        double eccentricity = 0.0;
        double perihelionArgumentDegrees = 0.0;
        double trueAnomalyRadians = 0.0;
    };

    explicit StarSystemTopView(QWidget *parent = nullptr);

    void setPlanets(const QVector<PlanetOrbit> &planets);
    void setStarParameters(double temperatureKelvin, double radiusSolar);
    void setSelectedIndex(int index);
    int selectedIndex() const;

signals:
    void planetSelected(int index);

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;

private:
    QVector<QPointF> planetScreenPositions(const QSize &size) const;
    double maxOrbitDistanceAu() const;

    QVector<PlanetOrbit> planets_;
    int selectedIndex_ = -1;
    double starTemperatureKelvin_ = 5772.0;
    double starRadiusSolar_ = 1.0;
};
