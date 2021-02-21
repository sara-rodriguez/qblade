#ifndef NOISEMODULE_H
#define NOISEMODULE_H

#include "../Module.h"
class NoiseDock;
class NoiseToolBar;
class NoiseSimulation;
class NoiseContextMenu;
class NoiseMenu;


class NoiseModule : public ModuleBase, public TwoDModule
{
	Q_OBJECT
	
public:
	NoiseModule(QMainWindow *mainWindow, QToolBar *toolbar);
	~NoiseModule();
	
	void addMainMenuEntries();
    QList<NewCurve*> prepareCurves (QString xAxis, QString yAxis, NewGraph::GraphType graphType,
									NewGraph::GraphType graphTypeMulti);
	QStringList getAvailableGraphVariables(bool xAxis);  // override from TwoDWidgetInterface
	QPair<ShowAsGraphInterface*,int> getHighlightDot(NewGraph::GraphType graphType);
	int getHighlightIndex(NewGraph::GraphType graphTypeMulti);
	QStringList prepareMissingObjectMessage();
	bool isColorByOpPoint();
	
private:
	void showAll();
	void hideAll();
	
	NoiseDock *m_dock;
	NoiseToolBar *m_toolBar;
	NoiseSimulation *m_shownSimulation;
	NoiseContextMenu *m_contextMenu;
	NoiseMenu *m_menu;	

    int index_qs3d=-1;//Sara
	
public slots:
	virtual void onActivationActionTriggered();  // override from ModuleBase
	virtual void onModuleChanged();  // override from ModuleBase
	void onHideDocks(bool hide);
    void onqs3dGraphs(bool graphs);//Sara
    void onqs3dGraph2d();//Sara
    void onqs3dGraphBlade();//Sara
    void onqs3dGraphRotor();//Sara
    void onqs3dGraphRotorLoops();//Sara
    void setShownSimulation(NoiseSimulation *newSimulation, bool forceReload = false);
	NoiseSimulation* getShownSimulation() { return m_shownSimulation; }

	void onNeedUpdate() { update(); }
    void reloadAllGraphs () { reloadAllGraphCurves(); }
};

extern NoiseModule *g_noiseModule;

#endif // NOISEMODULE_H
