#pragma once

#include <QImage>
#include <QString>

class HeightmapStorage {
public:
    QString storageDir() const;
    QString makeKey(const QString &planetName, int subdivisionLevel) const;
    QString makeFilePath(const QString &planetName, int subdivisionLevel) const;

    bool saveHeightmap(const QImage &image,
                       double scaleKm,
                       const QString &planetName,
                       int subdivisionLevel) const;
    bool loadHeightmap(QString *pathOut,
                       double *scaleOut,
                       const QString &planetName,
                       int subdivisionLevel) const;
    bool loadHeightmapImage(QImage *imageOut,
                            double *scaleOut,
                            const QString &filePath) const;

private:
    static QString sanitizeName(const QString &planetName);
};
