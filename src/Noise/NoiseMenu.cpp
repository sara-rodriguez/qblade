#include "NoiseMenu.h"

#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>

#include "NoiseModule.h"
#include "NoiseSimulation.h"

#include "../Store.h" //Sara
#include "NoiseCreatorDialog.h"//Sara


NoiseMenu::NoiseMenu(QMainWindow *parent, NoiseModule *module)
    : QMenu (parent)
{
    m_module = module;

    setTitle (tr("Noise Simulation"));
    connect (this, SIGNAL(aboutToShow()), SLOT(onAboutToShow()));

    m_modelValidityHint = new QAction("Model Validity Hint", this);
    connect(m_modelValidityHint, SIGNAL(triggered()), this, SLOT(onModelValidityHint()));
    addAction(m_modelValidityHint);

    m_exportNoise = new QAction("Export current Noise Simulation", this);
    connect(m_exportNoise, SIGNAL(triggered()), this, SLOT(onExportNoise()));
    addAction(m_exportNoise);

    //Sara
    m_exportqs3DNoiseLog = new QAction("Export current Quasi 3D Noise Log", this);
    connect(m_exportqs3DNoiseLog, SIGNAL(triggered()), this, SLOT(onExportqs3DNoiseLog()));
    addAction(m_exportqs3DNoiseLog);

    m_exportqs3DNoise_blade = new QAction("Export current Quasi 3D Noise Blade", this);
    connect(m_exportqs3DNoise_blade, SIGNAL(triggered()), this, SLOT(onExportqs3DNoise_blade()));
    addAction(m_exportqs3DNoise_blade);

    m_exportqs3DNoise_rotor_loops = new QAction("Export current Quasi 3D Noise Rotor", this);
    connect(m_exportqs3DNoise_rotor_loops, SIGNAL(triggered()), this, SLOT(onExportqs3DNoise_rotor()));
    addAction(m_exportqs3DNoise_rotor_loops);
    //Sara
}

void NoiseMenu::onAboutToShow() {
    const bool simulationAvailable = (  m_module->getShownSimulation() != NULL);

    //Sara
    NoiseCalculation *pNoiseCalculation = (NoiseCalculation *) g_mainFrame->m_pBEM;
    int index = pNoiseCalculation->UserSelectionMode();

    if (index==0){
        m_exportNoise->setEnabled(simulationAvailable);
        m_exportqs3DNoiseLog->setEnabled(false);
        m_exportqs3DNoise_blade->setEnabled(false);
        m_exportqs3DNoise_rotor_loops->setEnabled(false);
    }

    else if (index==1){
        m_exportNoise->setEnabled(simulationAvailable);
        m_exportqs3DNoiseLog->setEnabled(simulationAvailable);
        m_exportqs3DNoise_blade->setEnabled(simulationAvailable);
        m_exportqs3DNoise_rotor_loops->setEnabled(false);
    }

    else if (index==2){
        m_exportNoise->setEnabled(simulationAvailable);
        m_exportqs3DNoiseLog->setEnabled(simulationAvailable);
        m_exportqs3DNoise_blade->setEnabled(simulationAvailable);
        m_exportqs3DNoise_rotor_loops->setEnabled(simulationAvailable);
    }
    //Sara
}

void NoiseMenu::onExportNoise() {
    QString fileName = m_module->getShownSimulation()->getName() + ".txt";
    fileName.replace(' ', '_');
    fileName = QFileDialog::getSaveFileName(g_mainFrame, "Export Noise Simulation",
                                            g_mainFrame->m_ExportLastDirName + QDir::separator() + fileName,
                                            "Text File (*.txt)");
    if (!fileName.endsWith(".txt")) {
        fileName.append(".txt");
    }

    QFile file (fileName);
    g_mainFrame->m_ExportLastDirName = QFileInfo(file).absolutePath();
    if (file.open(QIODevice::WriteOnly)) {
        QTextStream fileStream (&file);
        m_module->getShownSimulation()->exportCalculation(fileStream);
    }
    file.close();
}

//Sara
void NoiseMenu::onExportqs3DNoiseLog() {
    QString fileName = m_module->getShownSimulation()->getName() + "-qs3D-log.csv";
    fileName.replace(' ', '_');
    fileName = QFileDialog::getSaveFileName(g_mainFrame, "Export Quasi 3D Noise Log",
                                            g_mainFrame->m_ExportLastDirName + QDir::separator() + fileName,
                                            "Comma Separated Values (*.csv)");
    if (!fileName.endsWith(".csv")) {
        fileName.append(".csv");
    }

    QFile file (fileName);
    g_mainFrame->m_ExportLastDirName = QFileInfo(file).absolutePath();
    if (file.open(QIODevice::WriteOnly)) {
        QTextStream fileStream (&file);
        m_module->getShownSimulation()->exportqs3DLog(fileStream);
    }
    file.close();
}

void NoiseMenu::onExportqs3DNoise_blade() {
    QString fileName = m_module->getShownSimulation()->getName() + "-qs3D-blade.csv";
    fileName.replace(' ', '_');
    fileName = QFileDialog::getSaveFileName(g_mainFrame, "Export quasi 3D Noise Blade Simulation",
                                            g_mainFrame->m_ExportLastDirName + QDir::separator() + fileName,
                                            "Comma Separated Values (*.csv)");
    if (!fileName.endsWith(".csv")) {
        fileName.append(".csv");
    }

    QFile file (fileName);
    g_mainFrame->m_ExportLastDirName = QFileInfo(file).absolutePath();
    if (file.open(QIODevice::WriteOnly)) {
        QTextStream fileStream (&file);
        m_module->getShownSimulation()->exportCalculationqs3DNoise_blade(fileStream);
    }
    file.close();
}

void NoiseMenu::onExportqs3DNoise_rotor() {
    QString fileName = m_module->getShownSimulation()->getName() + "-qs3D-rotor.csv";
    fileName.replace(' ', '_');
    fileName = QFileDialog::getSaveFileName(g_mainFrame, "Export quasi 3D Noise Rotor Simulation",
                                            g_mainFrame->m_ExportLastDirName + QDir::separator() + fileName,
                                            "Comma Separated Values (*.csv)");
    if (!fileName.endsWith(".csv")) {
        fileName.append(".csv");
    }

    QFile file (fileName);
    g_mainFrame->m_ExportLastDirName = QFileInfo(file).absolutePath();
    if (file.open(QIODevice::WriteOnly)) {
        QTextStream fileStream (&file);
        m_module->getShownSimulation()->exportCalculationqs3DNoise_rotor(fileStream);
    }
    file.close();
}
//Sara

void NoiseMenu::onModelValidityHint() {
//validation
    const QString message ("Airfoil TE noise model from Brooks, Pope & Marcolini, Airfoil Self-Noise and Prediction, "
                           "1989.\n\nThe original model was developed and validated for turbulent (tripped) flow up to "
                           "6 x 10⁵ Rec ≤ 1.5 × 10⁶, M < 0.21 and AOA up to 19.8°, for NACA 0012 airfoil.\n\nThe BPM directivity expression for High Frequency noise is not suitable for shallow angles (Θ → 180°). For details, see page 105 of (BPM, 1989).\n\nThe IAG Wind "
                           "Current validation range for the RDT-Modified LE Noise simulation showed good agreement in the range of 0.05588 ≤ M ≤ 0.18 and 1.3 x 10⁵ ≤ Rec ≥ 6 x 10⁵, with respect to experimental data (Bampanis et al., 2019, Juknevicius and Chong, 2018), for grid-generated turbulence structures."
                           "tunnel data (Herrig & Würz, 2008) showed good agreement with BPM prediction at Rec ~2.4 "
                           "× 10⁶ and M = 0.204, for peak Strouhal number and higher frequencies.");
    QMessageBox::information(g_mainFrame, "Model Validity Hint", message);
}
//Sara
