#pragma once

#include <QImage>
#include <QString>

class SurfaceHeightmap {
public:
    bool loadFromFile(const QString &path, double heightScaleKm);
    bool isValid() const;
    double heightKmAt(double latitudeDeg, double longitudeDeg) const;

private:
    quint16 sampleValue(int x, int y) const;

    QImage image_;
    double heightScaleKm_ = 0.0;
    bool isValid_ = false;
};
