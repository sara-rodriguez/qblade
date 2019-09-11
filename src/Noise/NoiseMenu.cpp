#include "NoiseMenu.h"

#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>

#include "NoiseModule.h"
#include "NoiseSimulation.h"


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
    m_exportqs3DNoiseLog = new QAction("Export quasi 3D Noise Log", this);
    connect(m_exportqs3DNoiseLog, SIGNAL(triggered()), this, SLOT(onExportqs3DNoiseLog()));
    addAction(m_exportqs3DNoiseLog);

    m_exportqs3DNoiseComplete = new QAction("Export quasi 3D Noise Complete", this);
    connect(m_exportqs3DNoiseComplete, SIGNAL(triggered()), this, SLOT(onExportqs3DNoiseComplete()));
    addAction(m_exportqs3DNoiseComplete);

    m_exportqs3DNoise = new QAction("Export current quasi 3D Noise Simulation", this);
    connect(m_exportqs3DNoise, SIGNAL(triggered()), this, SLOT(onExportqs3DNoise()));
    addAction(m_exportqs3DNoise);
    //Sara
}

void NoiseMenu::onAboutToShow() {
	const bool simulationAvailable = (m_module->getShownSimulation() != NULL);
	m_exportNoise->setEnabled(simulationAvailable);
    
    //Sara
        m_exportqs3DNoise->setEnabled(simulationAvailable);
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
    QString fileName = m_module->getShownSimulation()->getName() + "-3D-log.csv";
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


void NoiseMenu::onExportqs3DNoise() {
    QString fileName = m_module->getShownSimulation()->getName() + "-3D.csv";
    fileName.replace(' ', '_');
    fileName = QFileDialog::getSaveFileName(g_mainFrame, "Export quasi 3D Noise Simulation",
                                            g_mainFrame->m_ExportLastDirName + QDir::separator() + fileName,
                                            "Comma Separated Values (*.csv)");
    if (!fileName.endsWith(".csv")) {
        fileName.append(".csv");
    }

    QFile file (fileName);
    g_mainFrame->m_ExportLastDirName = QFileInfo(file).absolutePath();
    if (file.open(QIODevice::WriteOnly)) {
        QTextStream fileStream (&file);
        m_module->getShownSimulation()->exportqs3DCalculation(fileStream);
    }
    file.close();
}

void NoiseMenu::onExportqs3DNoiseComplete() {
    QString fileName = m_module->getShownSimulation()->getName() + "-3D-complete.csv";
    fileName.replace(' ', '_');
    fileName = QFileDialog::getSaveFileName(g_mainFrame, "Export quasi 3D Noise Simulation Complete",
                                            g_mainFrame->m_ExportLastDirName + QDir::separator() + fileName,
                                            "Comma Separated Values (*.csv)");
    if (!fileName.endsWith(".csv")) {
        fileName.append(".csv");
    }

    QFile file (fileName);
    g_mainFrame->m_ExportLastDirName = QFileInfo(file).absolutePath();
    if (file.open(QIODevice::WriteOnly)) {
        QTextStream fileStream (&file);
        m_module->getShownSimulation()->exportqs3DCalculationComplete(fileStream);
    }
    file.close();
}
//Sara

void NoiseMenu::onModelValidityHint() {
	const QString message ("Airfoil TE noise model from Brooks, Pope & Marcolini, Airfoil Self-Noise and Prediction, "
						   "1989.\nThe original model was developed and validated for turbulent (tripped) flow up to "
						   "Rec ≤ 1.5 × 10^6, M < 0.21 and AOA up to 19.8°, for NACA 0012 airfoil.\nThe IAG Wind "
						   "tunnel data (Herrig & Würz, 2008) showed good agreement with BPM prediction at Rec ~2.4 "
						   "× 10^6 and M = 0.204, for peak Strouhal number and higher frequencies.");
	QMessageBox::information(g_mainFrame, "Model Validity Hint", message);
}
