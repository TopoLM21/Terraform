#pragma once

#include <QWidget>

class SurfaceRealisticScaleWidget : public QWidget {
    Q_OBJECT

public:
    explicit SurfaceRealisticScaleWidget(QWidget *parent = nullptr);

protected:
    void paintEvent(QPaintEvent *event) override;
};
