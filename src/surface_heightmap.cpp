#include "surface_heightmap.h"

#include <QtMath>

#include <cmath>

namespace {
constexpr double kMaxU16 = 65535.0;
} // namespace

bool SurfaceHeightmap::loadFromFile(const QString &path, double heightScaleKm) {
    isValid_ = false;
    heightScaleKm_ = heightScaleKm;
    image_ = QImage();

    if (path.trimmed().isEmpty()) {
        return false;
    }

    QImage image(path);
    if (image.isNull()) {
        return false;
    }

    if (image.width() != image.height() * 2) {
        return false;
    }

    if (image.depth() < 16) {
        return false;
    }

    if (image.format() != QImage::Format_Grayscale16) {
        image = image.convertToFormat(QImage::Format_Grayscale16);
        if (image.isNull()) {
            return false;
        }
    }

    image_ = image;
    isValid_ = true;
    return true;
}

bool SurfaceHeightmap::isValid() const {
    return isValid_ && !image_.isNull();
}

quint16 SurfaceHeightmap::sampleValue(int x, int y) const {
    const int safeY = qBound(0, y, image_.height() - 1);
    const int safeX = (x % image_.width() + image_.width()) % image_.width();
    const auto *line =
        reinterpret_cast<const quint16 *>(image_.constScanLine(safeY));
    return line[safeX];
}

double SurfaceHeightmap::heightKmAt(double latitudeDeg, double longitudeDeg) const {
    if (!isValid()) {
        return 0.0;
    }

    const double lat = qBound(-90.0, latitudeDeg, 90.0);
    const double wrappedLon = std::fmod(longitudeDeg + 180.0, 360.0);
    const double lon = wrappedLon < 0.0 ? wrappedLon + 360.0 : wrappedLon;

    // Equirectangular развертка: широта идет сверху вниз, долготе соответствует полный оборот.
    const double u = lon / 360.0;
    const double v = (90.0 - lat) / 180.0;

    const double x = u * (image_.width() - 1);
    const double y = v * (image_.height() - 1);
    const int x0 = static_cast<int>(qFloor(x));
    const int y0 = static_cast<int>(qFloor(y));
    const int x1 = x0 + 1;
    const int y1 = y0 + 1;
    const double tx = x - x0;
    const double ty = y - y0;

    const double h00 = static_cast<double>(sampleValue(x0, y0));
    const double h10 = static_cast<double>(sampleValue(x1, y0));
    const double h01 = static_cast<double>(sampleValue(x0, y1));
    const double h11 = static_cast<double>(sampleValue(x1, y1));

    const double hx0 = h00 + (h10 - h00) * tx;
    const double hx1 = h01 + (h11 - h01) * tx;
    const double value = hx0 + (hx1 - hx0) * ty;

    // Нормируем диапазон 0..65535 в симметричный диапазон (-1..+1),
    // чтобы масштаб heightScaleKm задавал амплитуду рельефа в километрах.
    const double normalized = (value / kMaxU16) * 2.0 - 1.0;
    return normalized * heightScaleKm_;
}
