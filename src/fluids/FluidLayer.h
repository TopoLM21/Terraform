#pragma once

class FluidLayer {
public:
    FluidLayer() = default;
    FluidLayer(double densityKgPerM3, double thicknessMeters);

    double densityKgPerM3() const;
    void setDensityKgPerM3(double densityKgPerM3);

    double thicknessMeters() const;
    void setThicknessMeters(double thicknessMeters);

    // Масса слоя на единицу площади: m/A = ρ * h.
    double columnMassKgPerM2() const;

    // Гидростатическое давление слоя: P = ρ * g * h.
    double pressurePascal(double gravityMps2) const;

private:
    double densityKgPerM3_ = 0.0;
    double thicknessMeters_ = 0.0;
};
