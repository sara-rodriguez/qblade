/****************************************************************************

    SimuWidgetDMS Class
        Copyright (C) 2012 Juliane Wendler

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*****************************************************************************/

#ifndef SIMUWIDGETDMS_H
#define SIMUWIDGETDMS_H


#include <QtWidgets>
#include <QWidget>
#include "DMS.h"
#include "../Misc/LineCbBox.h"
#include "../Misc/LineButton.h"
#include "../Misc/LineDelegate.h"

class SimuWidgetDMS : public QWidget
{

    friend class QDMS;
    friend class MainFrame;

    Q_OBJECT

private slots:

    void OnStartDMS();
    void OnCreateDMS();
    void OnStartTurbineSimulation();
    void OnSetTurbineSimulationParameters();

    void OnShowCurve();
    void OnShowPoints();
    void OnShowOpPoint();
    void OnCurveWidth(int index);
    void OnCurveColor();
    void OnCurveStyle(int index);

	void OnCreateCharSim();
	void OnStartCharSim();
	void OnCheckDeltas();
	void OnFixCLicked();


public:
  
    SimuWidgetDMS(QWidget *parent);
    void SetupLayout();
    void Connect();

    void *m_pDMS;

    /////simulation params layout///
    QLabel *m_pctrlRho, *m_pctrlElements, *m_pctrlIteration, *m_pctrlEpsilon, *m_pctrlRelax, *m_pctrlVisc, *m_pctrlWindprofile;
    QLabel *m_pctrlRhoVal, *m_pctrlElementsVal, *m_pctrlIterationVal, *m_pctrlEpsilonVal, *m_pctrlRelaxVal, *m_pctrlViscVal, *m_pctrlWindprofileVal;
    QCheckBox *m_pctrlTipLoss, *m_pctrlAspectRatio, *m_pctrlVariable;
    QLabel *m_pctrlWind1Label, *m_pctrlWind2Label, *m_pctrlWindDeltaLabel;
    QStackedWidget *middleControls;

    //wing simu range layout
    QPushButton *m_pctrlCreateDMS, *m_pctrlStartDMS, *m_pctrlDefineTurbineSim, *m_pctrlStartTurbineSim, *m_pctrlDeleteDMS, *m_pctrlDeleteCDMS, *m_pctrlDeleteTDMS;
    QLabel *m_pctrlLSLabel, *m_pctrlLELabel, *m_pctrlLDLabel, *speed1, *speed2, *speed3;
    NumberEdit *m_pctrlLSLineEdit, *m_pctrlLDLineEdit, *m_pctrlLELineEdit, *m_pctrlWindspeed;
	NumberEdit *m_pctrlWind1, *m_pctrlWind2, *m_pctrlWindDelta;
	NumberEdit *WindStart, *WindEnd, *WindDelta, *RotStart, *RotEnd, *RotDelta, *PitchStart, *PitchEnd, *PitchDelta;
	QCheckBox *WindFixed, *RotFixed, *PitchFixed;
	QPushButton *CreateCharSim, *StartCharSim;
	QComboBox *Graph1Main, *Graph1Param, *Graph2Main, *Graph2Param, *Graph3Main, *Graph3Param, *Graph4Main, *Graph4Param;
	QLabel *WSpeed1, *WSpeed2, *WSpeed3;

    //Curve display layout///
    QCheckBox *m_pctrlShowCurve, *m_pctrlShowPoints, *m_pctrlShowOpPoint;
    LineCbBox *m_pctrlCurveStyle, *m_pctrlCurveWidth;
    LineButton *m_pctrlCurveColor;
    LineDelegate *m_pStyleDelegate, *m_pWidthDelegate;

private:

    void *m_pSimuWidgetDMS;

protected:




};


#endif // SIMUWIDGETDMS_H
