#pragma once

#include <QtCore/QMetaType>
#include <QtCore/QString>
#include <QtCore/QVector>

struct GasSpec {
    QString id;
    QString displayName;
    // Молярная масса в г/моль.
    double molarMass;
    bool isGreenhouse;
};

struct GasFraction {
    QString id;
    // Масса газа в гигатоннах (Gt).
    double massGigatons = 0.0;
    // Доля от общей массы (0..1).
    double share = 0.0;
};

class AtmosphereComposition {
public:
    // Массы газов хранятся в гигатоннах (Gt).
    void setMassGigatons(const QString &gasId, double massGigatons);
    double massGigatons(const QString &gasId) const;

    QVector<GasFraction> fractions() const;

    double totalMassGigatons() const;
    double totalMassTons() const;
    double totalMassKg() const;
    double totalPressureAtm(double planetMassEarths, double radiusKm) const;

private:
    QVector<GasFraction> m_gases;
};

Q_DECLARE_METATYPE(AtmosphereComposition)

QVector<GasSpec> availableGases();

// Рассчитать давление атмосферы в атм.
// massTons: масса атмосферы в метрических тоннах.
// planetMassEarths: масса планеты в массах Земли (M⊕).
// radiusKm: радиус планеты в километрах.
// Формулы:
//  g = G * M / R^2
//  P = (m_atm * g) / (4 * π * R^2)
//  1 атм = 101325 Па
// Единицы:
//  m_atm [кг], M [кг], R [м], g [м/с^2], P [Па].
//
// Для случая, когда масса уже в кг, используйте calculatePressureAtmFromKg.
double calculatePressureAtm(double massTons, double planetMassEarths, double radiusKm);

double calculatePressureAtmFromKg(double massKg, double planetMassEarths, double radiusKm);

// Рассчитать массу атмосферы в кг по давлению в атм (обратная формула).
// Формула выводится из P = (m_atm * g) / (4 * π * R^2).
double calculateAtmosphereMassKgFromPressureAtm(double pressureAtm,
                                                double planetMassEarths,
                                                double radiusKm);
