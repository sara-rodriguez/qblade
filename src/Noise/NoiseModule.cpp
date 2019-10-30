#include "NoiseModule.h"

#include <QSettings>
#include <QMenuBar>

#include "../MainFrame.h"
#include "NoiseToolBar.h"
#include "NoiseDock.h"
#include "NoiseContextMenu.h"
#include "../TwoDGraphMenu.h"
#include "../TwoDWidget.h"
#include "../Store.h"
#include "NoiseMenu.h"
#include "NoiseCreatorDialog.h"//Sara


NoiseModule::NoiseModule(QMainWindow *mainWindow, QToolBar *toolbar)
{
    m_globalModuleIndentifier = NOISEMODULE;
    m_shownSimulation = nullptr;

    on3dGraphs(false);
	
//Sara
//	m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL (dB)", true, false});
//    m_graph[1] = new NewGraph ("NoiseGraphTwo",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_alpha", true, false});
//	m_graph[2] = new NewGraph ("NoiseGraphThree", this, {NewGraph::Noise, "Freq [Hz]", "SPL_S", true, false});
//	m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_P", true, false});
//    m_graph[4] = new NewGraph ("NoiseGraphFive",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_LE (dB)", true, false});//Alexandre MOD

    //QSettings settings(QSettings::NativeFormat, QSettings::UserScope, "QBLADE");
    QSettings settings("qblade.ini", QSettings::IniFormat);//Sara
    setGraphArrangement(static_cast<TwoDWidgetInterface::GraphArrangement>
                        (settings.value("modules/NoiseModule/graphArrangement", TwoDWidgetInterface::Quint).toInt()));//Sara

	m_menu = new NoiseMenu (mainWindow, this);
	
	m_toolBar = new NoiseToolBar (mainWindow, this);
	m_dock = new NoiseDock ("Noise Simulation", mainWindow, 0, this);
	registrateAtToolbar("PNoise", "Predict the noise generation", ":/images/NoiseIcon.png", toolbar);
	
	m_contextMenu = new NoiseContextMenu (m_twoDWidget, this);  // NM TODO move this up to TwoDInterface?
	setContextMenu(m_contextMenu);

    connect(&g_noiseSimulationStore, SIGNAL(objectListChanged(bool)), this, SLOT(reloadAllGraphs()));
}

NoiseModule::~NoiseModule() {
    delete m_graph[0];
    delete m_graph[1];
    delete m_graph[2];
    delete m_graph[3];
    delete m_graph[4];

    //QSettings settings(QSettings::NativeFormat, QSettings::UserScope,"QBLADE");
    QSettings settings("qblade.ini", QSettings::IniFormat);//Sara
	settings.setValue(QString("modules/NoiseModule/graphArrangement"), getGraphArrangement());
}

void NoiseModule::addMainMenuEntries() {
	g_mainFrame->menuBar()->addMenu(m_graphMenu);
	g_mainFrame->menuBar()->addMenu(m_menu);
}

//Sara experiment new
QList<NewCurve *> NoiseModule::prepareCurves(QString xAxis, QString yAxis, NewGraph::GraphType graphType,NewGraph::GraphType /*graphTypeMulti*/) {
	QList<NewCurve*> curves;

    NoiseSimulation *pNoiseSim = (NoiseSimulation *) g_mainFrame->m_pBEM;
	for (int simIndex = 0; simIndex < g_noiseSimulationStore.size(); ++simIndex) {
		NoiseSimulation *simulation = g_noiseSimulationStore.at(simIndex);
                if (simulation->getSelectFrom() == NoiseParameter::OriginalBpm) {
			NewCurve *curve = simulation->newCurve(xAxis, yAxis, graphType, 0);
			if (curve) {
//				curves.append(curve);
			}
		} else {
            for (int i = 0; i < simulation->getAnalyzedOpPoints().size(); ++i) {
                NewCurve* curve = simulation->newCurve(xAxis, yAxis, graphType, i);
                if (curve) {
                    curves.append(curve);
                }
			}
		}
	}
    return curves;
}

QStringList NoiseModule::getAvailableGraphVariables(bool /*xAxis*/) {
	return NoiseSimulation::getAvailableVariables();
}

QPair<ShowAsGraphInterface *, int> NoiseModule::getHighlightDot(NewGraph::GraphType /*graphType*/) {
	return QPair<ShowAsGraphInterface*,int> (NULL, -1);  // function not available
}

int NoiseModule::getHighlightIndex(NewGraph::GraphType /*graphTypeMulti*/) {
	// return which index is to be painted bold
	int count = 0;
	for (int i = 0; i < g_noiseSimulationStore.size(); ++i) {
		if (g_noiseSimulationStore.at(i) == m_shownSimulation) {
			return count + std::max(0, m_toolBar->getShownOpPointIndex());
		} else {
			count += std::max(1, g_noiseSimulationStore.at(i)->getAnalyzedOpPoints().size());
		}
	}
	return -1;
}

QStringList NoiseModule::prepareMissingObjectMessage() {
	return NoiseSimulation::prepareMissingObjectMessage();
}

bool NoiseModule::isColorByOpPoint() {
	return m_dock->isColorByOpPoint();
}

void NoiseModule::showAll() {
	g_noiseSimulationStore.showAllCurves(true);
	m_dock->adjustShowCheckBox();
}

void NoiseModule::hideAll() {
	g_noiseSimulationStore.showAllCurves(false, m_shownSimulation);
	m_dock->adjustShowCheckBox();
}

void NoiseModule::onActivationActionTriggered() {
	ModuleBase::onActivationActionTriggered();
	showModule();
	g_mainFrame->switchToTwoDWidget();
	m_dock->show();
	m_toolBar->show();
}

void NoiseModule::onModuleChanged() {
	if (g_mainFrame->getCurrentModule() == this) {
		ModuleBase::onModuleChanged();
		hideModule();
		m_dock->hide();
		m_toolBar->hide();
	}
}

void NoiseModule::onHideDocks(bool hide) {
	m_dock->setVisible(!hide);
}

//Sara
void NoiseModule::on3dGraphs(bool graphs){

//QMainWindow *mainWindow;
//QToolBar *toolbar;
m_globalModuleIndentifier = NOISEMODULE;
m_shownSimulation = nullptr;

if (graphs)	{
//    m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_LEX (dB)", true, false});
//    m_graph[1] = new NewGraph ("NoiseGraphTwo",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_PX", true, false});
//    m_graph[2] = new NewGraph ("NoiseGraphThree", this, {NewGraph::Noise, "Freq [Hz]", "SPL_SX", true, false});
//    m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_alphaX", true, false});
//    m_graph[4] = new NewGraph ("NoiseGraphFive",  this, {NewGraph::Noise, "Freq [Hz]", "SPL[3D]X (dB)", true, false});//Alexandre MOD

    m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL[3D]X (dB)", true, false});
    m_graph[1] = new NewGraph ("NoiseGraphTwo",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_alpha[3D]X", true, false});
    m_graph[2] = new NewGraph ("NoiseGraphThree", this, {NewGraph::Noise, "Freq [Hz]", "SPL_S[3D]X", true, false});
    m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_P[3D]X", true, false});
    m_graph[4] = new NewGraph ("NoiseGraphFive",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_LE[3D]X (dB)", true, false});//Alexandre MOD

} else {
    m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL (dB)", true, false});
    m_graph[1] = new NewGraph ("NoiseGraphTwo",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_alpha", true, false});
    m_graph[2] = new NewGraph ("NoiseGraphThree", this, {NewGraph::Noise, "Freq [Hz]", "SPL_S", true, false});
    m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_P", true, false});
    m_graph[4] = new NewGraph ("NoiseGraphFive",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_LE (dB)", true, false});//Alexandre MOD

//    m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL[3D] (dB)", true, false});
//    m_graph[1] = new NewGraph ("NoiseGraphTwo",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_alpha[3D]", true, false});
//    m_graph[2] = new NewGraph ("NoiseGraphThree", this, {NewGraph::Noise, "Freq [Hz]", "SPL_S[3D]", true, false});
//    m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_P[3D]", true, false});
//    m_graph[4] = new NewGraph ("NoiseGraphFive",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_LE[3D] (dB)", true, false});//Alexandre MOD
}

QSettings settings("qblade.ini", QSettings::IniFormat);//Sara

setGraphArrangement(static_cast<TwoDWidgetInterface::GraphArrangement>(settings.value("modules/NoiseModule/graphArrangement", TwoDWidgetInterface::Quint).toInt()));//Sara

//m_menu = new NoiseMenu (mainWindow, this);

//m_toolBar = new NoiseToolBar (mainWindow, this);
//m_dock = new NoiseDock ("Noise Simulation", mainWindow, 0, this);

m_contextMenu = new NoiseContextMenu (m_twoDWidget, this);  // NM TODO move this up to TwoDInterface?
setContextMenu(m_contextMenu);

connect(&g_noiseSimulationStore, SIGNAL(objectListChanged(bool)), this, SLOT(reloadAllGraphs()));
}
//Sara

void NoiseModule::setShownSimulation(NoiseSimulation *newSimulation, bool forceReload) {
	if (forceReload || m_shownSimulation != newSimulation) {
		m_shownSimulation = newSimulation;
		m_dock->setShownObject(m_shownSimulation);
		m_toolBar->setShownSimulation(m_shownSimulation);
		reloadForGraphType(NewGraph::Noise);
	}
}



NoiseModule *g_noiseModule;
