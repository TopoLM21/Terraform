#pragma once

#include "rotation_mode.h"

#include <QtWidgets/QWidget>

class ModeIllustrationWidget : public QWidget {
public:
    explicit ModeIllustrationWidget(QWidget *parent = nullptr);

    void setRotationMode(RotationMode mode);
    RotationMode rotationMode() const;

    QSize sizeHint() const override;

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    RotationMode rotationMode_ = RotationMode::Normal;
};
