#include "atmosphere_widget.h"

#include "atmosphere_chart_widget.h"

#include <QtCore/QLocale>
#include <QtCore/QSignalBlocker>
#include <QtGui/QDoubleValidator>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QStyledItemDelegate>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>

#include <limits>

namespace {
constexpr int kColumnGas = 0;
constexpr int kColumnMass = 1;
constexpr int kColumnPressure = 2;
constexpr int kColumnShare = 3;
constexpr double kKgPerGigaton = 1.0e12;
constexpr double kGramsPerKg = 1000.0;

class NumericColumnDelegate : public QStyledItemDelegate {
public:
    explicit NumericColumnDelegate(QObject *parent = nullptr)
        : QStyledItemDelegate(parent) {}

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const override {
        QWidget *editor = QStyledItemDelegate::createEditor(parent, option, index);
        auto *lineEdit = qobject_cast<QLineEdit *>(editor);
        if (lineEdit) {
            auto *validator = new QDoubleValidator(0.0, std::numeric_limits<double>::max(), 8, lineEdit);
            validator->setNotation(QDoubleValidator::StandardNotation);
            lineEdit->setValidator(validator);
        }
        return editor;
    }
};

QString formatNumber(double value, int decimals = 3) {
    return QLocale().toString(value, 'f', decimals);
}
} // namespace

AtmosphereWidget::AtmosphereWidget(QWidget *parent, bool showTable)
    : QGroupBox(QStringLiteral("Атмосфера"), parent) {
    gases_ = availableGases();

    table_ = new QTableWidget(this);
    table_->setColumnCount(4);
    table_->setHorizontalHeaderLabels({
        QStringLiteral("Газ"),
        QStringLiteral("Масса (Гт)"),
        QStringLiteral("Давление (атм)"),
        QStringLiteral("Доля (%)")
    });
    table_->setItemDelegateForColumn(kColumnMass, new NumericColumnDelegate(table_));
    table_->setItemDelegateForColumn(kColumnPressure, new NumericColumnDelegate(table_));
    table_->horizontalHeader()->setSectionResizeMode(kColumnGas, QHeaderView::Stretch);
    table_->horizontalHeader()->setSectionResizeMode(kColumnMass, QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setSectionResizeMode(kColumnPressure, QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setSectionResizeMode(kColumnShare, QHeaderView::ResizeToContents);
    table_->verticalHeader()->setVisible(false);
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->setEnabled(false);
    if (!showTable) {
        table_->setVisible(false);
    }

    populateTable();

    chartWidget_ = new AtmosphereChartWidget(gases_, this);

    totalMassLabel_ = new QLabel(QStringLiteral("—"), this);
    pressureLabel_ = new QLabel(QStringLiteral("—"), this);
    meanMolarMassLabel_ = new QLabel(QStringLiteral("—"), this);

    auto *summaryLayout = new QFormLayout();
    summaryLayout->addRow(QStringLiteral("Масса атмосферы:"), totalMassLabel_);
    summaryLayout->addRow(QStringLiteral("Давление (атм):"), pressureLabel_);
    summaryLayout->addRow(QStringLiteral("Средняя молекулярная масса (г/моль):"),
                          meanMolarMassLabel_);

    auto *layout = new QVBoxLayout(this);
    layout->addWidget(table_);
    layout->addWidget(chartWidget_, 1);
    layout->addLayout(summaryLayout);
    setLayout(layout);

    connect(table_, &QTableWidget::itemChanged, this, [this](QTableWidgetItem *item) {
        if (!item) {
            return;
        }
        if (item->column() == kColumnMass) {
            normalizeMassItem(item->row());
            updateAllPressures();
            updateAllShares();
            updateSummary();
            if (chartWidget_) {
                chartWidget_->setComposition(composition(true));
            }
            return;
        }
        if (item->column() == kColumnPressure) {
            normalizePressureItem(item->row());
            updateAllShares();
            updateSummary();
            if (chartWidget_) {
                chartWidget_->setComposition(composition(true));
            }
        }
    });

    updateAllShares();
    updateAllPressures();
    updateSummary();
    if (chartWidget_) {
        chartWidget_->setComposition(composition(true));
    }
}

void AtmosphereWidget::setPlanetParameters(double massEarths, double radiusKm) {
    if (massEarths <= 0.0 || radiusKm <= 0.0) {
        clearPlanetParameters();
        return;
    }
    planetMassEarths_ = massEarths;
    planetRadiusKm_ = radiusKm;
    hasPlanetParameters_ = true;
    table_->setEnabled(true);
    updateAllPressures();
    updateSummary();
    if (chartWidget_) {
        chartWidget_->setPlanetParameters(massEarths, radiusKm);
    }
}

void AtmosphereWidget::setComposition(const AtmosphereComposition &newComposition) {
    const QSignalBlocker blocker(table_);
    for (int row = 0; row < gases_.size(); ++row) {
        const auto &gas = gases_.at(row);
        const double massGigatons = newComposition.massGigatons(gas.id);
        if (auto *item = table_->item(row, kColumnMass)) {
            item->setText(formatNumber(massGigatons, 3));
        }
    }
    updateAllShares();
    updateAllPressures();
    updateSummary();
    if (chartWidget_) {
        chartWidget_->setComposition(newComposition);
    }
}

void AtmosphereWidget::clearPlanetParameters() {
    planetMassEarths_ = 0.0;
    planetRadiusKm_ = 0.0;
    hasPlanetParameters_ = false;
    table_->setEnabled(false);
    updateAllPressures();
    updateSummary();
    if (chartWidget_) {
        chartWidget_->clearPlanetParameters();
    }
}

AtmosphereComposition AtmosphereWidget::composition(bool includeZeroes) const {
    AtmosphereComposition composition;
    for (int row = 0; row < table_->rowCount(); ++row) {
        const double massGigatons = rowMassGigatons(row);
        if (!includeZeroes && massGigatons <= 0.0) {
            continue;
        }
        composition.setMassGigatons(gases_.at(row).id, massGigatons);
    }
    return composition;
}

void AtmosphereWidget::populateTable() {
    const QSignalBlocker blocker(table_);
    table_->setRowCount(gases_.size());

    for (int row = 0; row < gases_.size(); ++row) {
        const auto &gas = gases_.at(row);
        auto *gasItem = new QTableWidgetItem(gas.displayName);
        gasItem->setFlags(gasItem->flags() & ~Qt::ItemIsEditable);
        table_->setItem(row, kColumnGas, gasItem);

        auto *massItem = new QTableWidgetItem(QStringLiteral("0"));
        table_->setItem(row, kColumnMass, massItem);

        auto *pressureItem = new QTableWidgetItem(QStringLiteral("—"));
        table_->setItem(row, kColumnPressure, pressureItem);

        auto *shareItem = new QTableWidgetItem(QStringLiteral("0"));
        shareItem->setFlags(shareItem->flags() & ~Qt::ItemIsEditable);
        table_->setItem(row, kColumnShare, shareItem);
    }
}

void AtmosphereWidget::updateAllShares() {
    double totalMassGigatons = 0.0;
    for (int row = 0; row < table_->rowCount(); ++row) {
        totalMassGigatons += rowMassGigatons(row);
    }

    const QSignalBlocker blocker(table_);
    for (int row = 0; row < table_->rowCount(); ++row) {
        const double massGigatons = rowMassGigatons(row);
        const double share = totalMassGigatons > 0.0 ? (massGigatons / totalMassGigatons) * 100.0 : 0.0;
        if (auto *item = table_->item(row, kColumnShare)) {
            item->setText(formatNumber(share, 2));
        }
    }
}

void AtmosphereWidget::updateAllPressures() {
    const QSignalBlocker blocker(table_);
    for (int row = 0; row < table_->rowCount(); ++row) {
        auto *pressureItem = table_->item(row, kColumnPressure);
        if (!pressureItem) {
            continue;
        }
        if (!hasPlanetParameters_) {
            pressureItem->setText(QStringLiteral("—"));
            continue;
        }

        const double massKg = rowMassGigatons(row) * kKgPerGigaton;
        const double pressureAtm = calculatePressureAtmFromKg(massKg, planetMassEarths_, planetRadiusKm_);
        pressureItem->setText(formatNumber(pressureAtm, 4));
    }
}

void AtmosphereWidget::updateSummary() {
    double totalMassGigatons = 0.0;
    double totalMassKg = 0.0;
    double totalMassGrams = 0.0;
    double totalMoles = 0.0;

    for (int row = 0; row < table_->rowCount(); ++row) {
        const double massGigatons = rowMassGigatons(row);
        // Массы хранятся в гигатоннах (Gt), для физических расчетов переводим в кг и г.
        const double massKg = massGigatons * kKgPerGigaton;
        const double massGrams = massKg * kGramsPerKg;

        totalMassGigatons += massGigatons;
        totalMassKg += massKg;
        totalMassGrams += massGrams;

        const double molarMass = gases_.at(row).molarMass;
        if (molarMass > 0.0) {
            totalMoles += massGrams / molarMass;
        }
    }

    const QString massText = QStringLiteral("%1 Гт (%2 т)")
                                 .arg(formatNumber(totalMassGigatons, 3))
                                 .arg(QLocale().toString(totalMassGigatons * 1.0e9, 'g', 6));
    totalMassLabel_->setText(massText);

    if (hasPlanetParameters_) {
        // Давление на поверхности рассчитывается по модели атмосферы:
        // P = (m_atm * g) / (4 * π * R^2), где g = G * M / R^2.
        const double pressureAtm = calculatePressureAtmFromKg(totalMassKg, planetMassEarths_, planetRadiusKm_);
        pressureLabel_->setText(formatNumber(pressureAtm, 4));
    } else {
        pressureLabel_->setText(QStringLiteral("—"));
    }

    if (totalMoles > 0.0) {
        // Средняя молекулярная масса рассчитана через общее число молей
        // (масса -> моли) и суммарную массу смеси.
        const double meanMolarMass = totalMassGrams / totalMoles;
        meanMolarMassLabel_->setText(formatNumber(meanMolarMass, 3));
    } else {
        meanMolarMassLabel_->setText(QStringLiteral("—"));
    }
}

void AtmosphereWidget::normalizeMassItem(int row) {
    auto *item = table_->item(row, kColumnMass);
    if (!item) {
        return;
    }

    const double value = parseMassText(item->text());
    const QSignalBlocker blocker(table_);
    item->setText(formatNumber(value, 3));
}

void AtmosphereWidget::normalizePressureItem(int row) {
    auto *item = table_->item(row, kColumnPressure);
    if (!item) {
        return;
    }

    const double pressureAtm = parsePressureText(item->text());
    if (!hasPlanetParameters_) {
        const QSignalBlocker blocker(table_);
        item->setText(QStringLiteral("—"));
        return;
    }

    // Обратный пересчет выполняется из давления в атм в массу в кг:
    // P [атм] -> P [Па], m_atm [кг] = (P * 4πR^2) / g, где g = G*M/R^2.
    const double massKg = calculateAtmosphereMassKgFromPressureAtm(
        pressureAtm, planetMassEarths_, planetRadiusKm_);
    const double massGigatons = massKg / kKgPerGigaton;

    const QSignalBlocker blocker(table_);
    item->setText(formatNumber(pressureAtm, 4));
    if (auto *massItem = table_->item(row, kColumnMass)) {
        massItem->setText(formatNumber(massGigatons, 3));
    }
    updateAllPressures();
}

double AtmosphereWidget::parseMassText(const QString &text) const {
    const QString trimmed = text.trimmed();
    if (trimmed.isEmpty()) {
        return 0.0;
    }
    bool ok = false;
    const double value = QLocale().toDouble(trimmed, &ok);
    if (!ok || value < 0.0) {
        return 0.0;
    }
    return value;
}

double AtmosphereWidget::parsePressureText(const QString &text) const {
    const QString trimmed = text.trimmed();
    if (trimmed.isEmpty()) {
        return 0.0;
    }
    bool ok = false;
    const double value = QLocale().toDouble(trimmed, &ok);
    if (!ok || value < 0.0) {
        return 0.0;
    }
    return value;
}

double AtmosphereWidget::rowMassGigatons(int row) const {
    const auto *item = table_->item(row, kColumnMass);
    if (!item) {
        return 0.0;
    }
    return parseMassText(item->text());
}
