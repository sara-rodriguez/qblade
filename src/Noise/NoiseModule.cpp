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
#include "NoiseCalculation.h"//Sara

NoiseModule::NoiseModule(QMainWindow *mainWindow, QToolBar *toolbar)
{
    m_globalModuleIndentifier = NOISEMODULE;
    m_shownSimulation = nullptr;

//Sara
    onqs3dGraphs(false);//Sara
	
    //QSettings settings(QSettings::NativeFormat, QSettings::UserScope, "QBLADE");
    QSettings settings("qblade.ini", QSettings::IniFormat);//Sara
    setGraphArrangement(static_cast<TwoDWidgetInterface::GraphArrangement>
                        (settings.value("modules/NoiseModule/graphArrangement", TwoDWidgetInterface::Oct).toInt()));
//Sara

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
    delete m_graph[5];//Sara
    delete m_graph[6];//Sara
    delete m_graph[7];//Sara

    //QSettings settings(QSettings::NativeFormat, QSettings::UserScope,"QBLADE");
    QSettings settings("qblade.ini", QSettings::IniFormat);//Sara
	settings.setValue(QString("modules/NoiseModule/graphArrangement"), getGraphArrangement());
}

void NoiseModule::addMainMenuEntries() {
	g_mainFrame->menuBar()->addMenu(m_graphMenu);
	g_mainFrame->menuBar()->addMenu(m_menu);
}

QList<NewCurve *> NoiseModule::prepareCurves(QString xAxis, QString yAxis, NewGraph::GraphType graphType,NewGraph::GraphType /*graphTypeMulti*/) {
	QList<NewCurve*> curves;

	for (int simIndex = 0; simIndex < g_noiseSimulationStore.size(); ++simIndex) {
		NoiseSimulation *simulation = g_noiseSimulationStore.at(simIndex);
                if (simulation->getSelectFrom() == NoiseParameter::OriginalBpm) {
			NewCurve *curve = simulation->newCurve(xAxis, yAxis, graphType, 0);
			if (curve) {
                curves.append(curve);
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
    int index = 0;
    NoiseCalculation *pNoiseCalculation = (NoiseCalculation *) g_mainFrame->m_pBEM;
    int user_sel = pNoiseCalculation->user_sel;
    bool qs3d_check = pNoiseCalculation->user_qs3d_check;
    if(!qs3d_check){index = 0;}
    else {index = user_sel+1;}

    if(index==0){return NoiseSimulation::getAvailableVariables();}
    else if(index==1){return NoiseSimulation::getAvailableVariables_blade();}
    else if(index==2){return NoiseSimulation::getAvailableVariables_rotor();}
    else{return NoiseSimulation::getAvailableVariables();}
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
void NoiseModule::onqs3dGraphs(bool){
++index_qs3d;

NoiseCalculation *pNoiseCalculation = (NoiseCalculation *) g_mainFrame->m_pBEM;
    int user_sel = pNoiseCalculation->user_sel;
    bool qs3d_check = pNoiseCalculation->user_qs3d_check;
if (!qs3d_check){
onqs3dGraph2d();  index_qs3d=-1;
}
else if (user_sel==0){
if(index_qs3d==0){onqs3dGraph2d();}
if(index_qs3d==1){onqs3dGraphBlade(); index_qs3d=-1;}
}
else if (user_sel==1){
if(index_qs3d==0){onqs3dGraph2d();}
if(index_qs3d==1){onqs3dGraphBlade();}
if(index_qs3d==2){onqs3dGraphRotorLoops(); index_qs3d=-1;}
}

QSettings settings("qblade.ini", QSettings::IniFormat);

setGraphArrangement(static_cast<TwoDWidgetInterface::GraphArrangement>(settings.value("modules/NoiseModule/graphArrangement", TwoDWidgetInterface::Oct).toInt()));//Sara
//Sara

m_contextMenu = new NoiseContextMenu (m_twoDWidget, this);  // NM TODO move this up to TwoDInterface?
setContextMenu(m_contextMenu);

connect(&g_noiseSimulationStore, SIGNAL(objectListChanged(bool)), this, SLOT(reloadAllGraphs()));
}

void NoiseModule::onqs3dGraph2d(){
        m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_alpha", true, false});
        m_graph[1] = new NewGraph ("NoiseGraphTwo", this, {NewGraph::Noise, "Freq [Hz]", "SPL_S", true, false});
        m_graph[2] = new NewGraph ("NoiseGraphThree",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_P", true, false});
        m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_LE (dB)", true, false});//Alexandre MOD
        m_graph[4] = new NewGraph ("NoiseGraphFive",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_LBL_VS (dB)", true, false});//Sara
        m_graph[5] = new NewGraph ("NoiseGraphSix",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_blunt (dB)", true, false});//Sara
        m_graph[6] = new NewGraph ("NoiseGraphSeven",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_tipvortex (dB)", true, false});//Sara
        m_graph[7] = new NewGraph ("NoiseGraphEight",   this, {NewGraph::Noise, "Freq [Hz]", "SPL (dB)", true, false});
}

void NoiseModule::onqs3dGraphBlade(){
        m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_alpha_blade[qs3D]", true, false});
        m_graph[1] = new NewGraph ("NoiseGraphTwo", this, {NewGraph::Noise, "Freq [Hz]", "SPL_S_blade[qs3D]", true, false});
        m_graph[2] = new NewGraph ("NoiseGraphThree",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_P_blade[qs3D]", true, false});
        m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_LE_blade[qs3D] (dB)", true, false});//Alexandre MOD
        m_graph[4] = new NewGraph ("NoiseGraphFive",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_LBL_VS_blade[qs3D] (dB)", true, false});//Sara
        m_graph[5] = new NewGraph ("NoiseGraphSix",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_blunt_blade[qs3D] (dB)", true, false});//Sara
        m_graph[6] = new NewGraph ("NoiseGraphSeven",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_tipvortex_blade[qs3D] (dB)", true, false});//Sara
        m_graph[7] = new NewGraph ("NoiseGraphEight",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_blade[qs3D] (dB)", true, false});
}

void NoiseModule::onqs3dGraphRotor(){
        m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_alpha_rotor[qs3D]", true, false});
        m_graph[1] = new NewGraph ("NoiseGraphTwo", this, {NewGraph::Noise, "Freq [Hz]", "SPL_S_rotor[qs3D]", true, false});
        m_graph[2] = new NewGraph ("NoiseGraphThree",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_P_rotor[qs3D]", true, false});
        m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_LE_rotor[qs3D] (dB)", true, false});//Alexandre MOD
        m_graph[4] = new NewGraph ("NoiseGraphFive",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_LBL_VS_rotor[qs3D] (dB)", true, false});//Sara
        m_graph[5] = new NewGraph ("NoiseGraphSix",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_blunt_rotor[qs3D] (dB)", true, false});//Sara
        m_graph[6] = new NewGraph ("NoiseGraphSeven",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_tipvortex_rotor[qs3D] (dB)", true, false});//Sara
        m_graph[7] = new NewGraph ("NoiseGraphEight",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_rotor[qs3D] (dB)", true, false});
}

void NoiseModule::onqs3dGraphRotorLoops(){
        m_graph[0] = new NewGraph ("NoiseGraphOne",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_alpha_rotor[qs3D]", true, false});
        m_graph[1] = new NewGraph ("NoiseGraphTwo", this, {NewGraph::Noise, "Freq [Hz]", "SPL_S_rotor[qs3D]", true, false});
        m_graph[2] = new NewGraph ("NoiseGraphThree",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_P_rotor[qs3D]", true, false});
        m_graph[3] = new NewGraph ("NoiseGraphFour",  this, {NewGraph::Noise, "Freq [Hz]", "SPL_LE_rotor[qs3D] (dB)", true, false});//Alexandre MOD
        m_graph[4] = new NewGraph ("NoiseGraphFive",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_LBL_VS_rotor[qs3D] (dB)", true, false});//Sara
        m_graph[5] = new NewGraph ("NoiseGraphSix",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_blunt_rotor[qs3D] (dB)", true, false});//Sar
        m_graph[6] = new NewGraph ("NoiseGraphSeven",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_tipvortex_rotor[qs3D] (dB)", true, false});//Sara
        m_graph[7] = new NewGraph ("NoiseGraphEight",   this, {NewGraph::Noise, "Freq [Hz]", "SPL_rotor[qs3D] (dB)", true, false});
}
//Sara

void NoiseModule::setShownSimulation(NoiseSimulation *newSimulation, bool forceReload) {
    if (forceReload || m_shownSimulation != newSimulation) {
        onqs3dGraph2d();//Sara
        QSettings settings("qblade.ini", QSettings::IniFormat); //Sara
        setGraphArrangement(static_cast<TwoDWidgetInterface::GraphArrangement>(settings.value("modules/NoiseModule/graphArrangement", TwoDWidgetInterface::Oct).toInt()));//Sara


		m_shownSimulation = newSimulation;
		m_dock->setShownObject(m_shownSimulation);
		m_toolBar->setShownSimulation(m_shownSimulation);
		reloadForGraphType(NewGraph::Noise);
	}
}



NoiseModule *g_noiseModule;
