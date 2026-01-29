#include "heightmap_storage.h"

#include <QDataStream>
#include <QDir>
#include <QFile>
#include <QStandardPaths>

#include <cstring>

namespace {
constexpr int kMagicSize = 4;
constexpr char kMagic[kMagicSize] = {'H', 'M', 'A', 'P'};
constexpr quint16 kFormatVersion = 1;

qint64 headerByteSize() {
    return kMagicSize + static_cast<qint64>(sizeof(kFormatVersion)) +
           static_cast<qint64>(sizeof(quint32)) * 2 +
           static_cast<qint64>(sizeof(double));
}
} // namespace

QString HeightmapStorage::storageDir() const {
    return QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) +
           "/heightmaps";
}

QString HeightmapStorage::sanitizeName(const QString &planetName) {
    QString name = planetName.trimmed();
    name.replace(' ', '_');
    return name;
}

QString HeightmapStorage::makeKey(const QString &planetName, int subdivisionLevel) const {
    return QString("%1_subdiv%2")
        .arg(sanitizeName(planetName))
        .arg(subdivisionLevel);
}

QString HeightmapStorage::makeFilePath(const QString &planetName, int subdivisionLevel) const {
    const QString fileName = makeKey(planetName, subdivisionLevel) + ".h16";
    return QDir(storageDir()).filePath(fileName);
}

bool HeightmapStorage::saveHeightmap(const QImage &image,
                                     double scaleKm,
                                     const QString &planetName,
                                     int subdivisionLevel) const {
    if (image.isNull()) {
        return false;
    }

    if (image.width() != image.height() * 2) {
        return false;
    }

    QImage grayscale = image;
    if (grayscale.format() != QImage::Format_Grayscale16) {
        grayscale = grayscale.convertToFormat(QImage::Format_Grayscale16);
        if (grayscale.isNull()) {
            return false;
        }
    }

    QDir dir;
    const QString dirPath = storageDir();
    if (!dir.mkpath(dirPath)) {
        return false;
    }

    const QString filePath = makeFilePath(planetName, subdivisionLevel);
    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly)) {
        return false;
    }

    QDataStream stream(&file);
    stream.setByteOrder(QDataStream::LittleEndian);
    stream.setVersion(QDataStream::Qt_5_15);

    // Формат файла:
    //   - 4 байта: magic "HMAP"
    //   - quint16: версия
    //   - quint32: ширина
    //   - quint32: высота
    //   - double: масштаб рельефа (км)
    //   - далее массив quint16 (equirectangular 2:1)
    // 16-бит выбраны, чтобы сохранить точность высот без большого размера,
    // и чтобы напрямую работать с QImage::Format_Grayscale16.
    stream.writeRawData(kMagic, kMagicSize);
    stream << kFormatVersion;
    stream << static_cast<quint32>(grayscale.width());
    stream << static_cast<quint32>(grayscale.height());
    stream << scaleKm;

    for (int y = 0; y < grayscale.height(); ++y) {
        const auto *line =
            reinterpret_cast<const quint16 *>(grayscale.constScanLine(y));
        for (int x = 0; x < grayscale.width(); ++x) {
            stream << line[x];
        }
    }

    return stream.status() == QDataStream::Ok;
}

bool HeightmapStorage::loadHeightmap(QString *pathOut,
                                     double *scaleOut,
                                     const QString &planetName,
                                     int subdivisionLevel) const {
    const QString filePath = makeFilePath(planetName, subdivisionLevel);
    QFile file(filePath);
    if (!file.exists()) {
        return false;
    }

    if (!file.open(QIODevice::ReadOnly)) {
        return false;
    }

    QDataStream stream(&file);
    stream.setByteOrder(QDataStream::LittleEndian);
    stream.setVersion(QDataStream::Qt_5_15);

    char magic[kMagicSize] = {};
    if (stream.readRawData(magic, kMagicSize) != kMagicSize) {
        return false;
    }

    if (std::memcmp(magic, kMagic, kMagicSize) != 0) {
        return false;
    }

    quint16 version = 0;
    quint32 width = 0;
    quint32 height = 0;
    double scaleKm = 0.0;
    stream >> version >> width >> height >> scaleKm;

    if (stream.status() != QDataStream::Ok) {
        return false;
    }

    if (version != kFormatVersion) {
        return false;
    }

    if (width == 0 || height == 0 || width != height * 2) {
        return false;
    }

    const qint64 expectedSize = headerByteSize() +
        static_cast<qint64>(width) * static_cast<qint64>(height) *
            static_cast<qint64>(sizeof(quint16));
    if (file.size() < expectedSize) {
        return false;
    }

    if (pathOut) {
        *pathOut = filePath;
    }
    if (scaleOut) {
        *scaleOut = scaleKm;
    }

    return true;
}
