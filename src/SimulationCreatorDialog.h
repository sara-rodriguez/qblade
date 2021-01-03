#ifndef SIMULATIONCREATORDIALOG_H
#define SIMULATIONCREATORDIALOG_H

class QGroupBox;

#include "CreatorDialog.h"
#include "ParameterViewer.h"


template <class ParameterGroup>
class SimulationCreatorDialog : public CreatorDialog, public ParameterViewer<ParameterGroup>
{
//Sara
public:
    double rho=1.225;
    double temperature=288.15;
    double relax_factor=0.35;
    double viscosity=0.00001647;
//Sara

protected:
	typedef ParameterGroup P;

	SimulationCreatorDialog();
	QGroupBox* constructParameterBox(QString defaultName);
	QGroupBox* constructCorrectionsBox();
	virtual void onRestoreDefaultClicked();
	virtual void loadValuesFromSettings();
	virtual void saveValuesToSettings();
	
	QVBoxLayout *m_leftContentVBox, *m_rightContentVBox;
};

#endif // SIMULATIONCREATORDIALOG_H
