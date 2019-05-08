/****************************************************************************

	QXDirect Class
	Copyright (C) 2008-2009 Andre Deperrois adeperrois@xflr5.com

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

#ifndef QXDIRECT_H
#define QXDIRECT_H

#include <QPixmap>
#include <QSettings>
#include <QLabel>
#include <QRadioButton>
#include <QCheckBox>
#include <QStackedWidget>
#include "../Misc/LineButton.h"
#include "../Misc/LineCbBox.h"
#include "../Misc/LineDelegate.h"
#include "../Misc/NumberEdit.h"
#include "../Objects/Polar.h"
#include "../Objects/Foil.h"
#include "../Objects/OpPoint.h"
#include "../Graph/QGraph.h"
#include "XFoil.h"
#include "XFoilAnalysisDlg.h"
#include "../Graph/GraphDlg.h"
#include "BatchDlg.h"
#include "BatchThreadDlg.h"
#include "FoilPolarDlg.h"
#include "XFoilAdvancedDlg.h"
#include "../Misc/PolarFilterDlg.h"

#include "TwoDPanelDlg.h"
#include "InterpolateFoilsDlg.h"
#include "NacaFoilDlg.h"
#include "FoilCoordDlg.h"
#include "FoilGeomDlg.h"
#include "TEGapDlg.h"
#include "LEDlg.h"
#include "FlapDlg.h"
#include "CAddDlg.h"
#include "XDirectStyleDlg.h"
#include "ManageFoilsDlg.h"


class QXDirect : public QWidget
{
	friend class MainFrame;
	friend class TwoDWidget;
	friend class XFoilAnalysisDlg;
	friend class BatchDlg;
	friend class BatchThreadDlg;
	friend class FoilPolarDlg;
	friend class QFoilGeomDlg;
	friend class QTEGapDlg;
	friend class QCAddDlg;
	friend class EditPlrDlg;
	friend class XDirectStyleDlg;
	friend class DisplaySettingsDlg;
    friend class QAFoil;
        ////////new code DM///////////
        friend class QBEM;
        ///////end new code DM//////////

    Q_OBJECT

public:
	QXDirect(QWidget *parent = NULL);
	~QXDirect();

public slots:
	void UpdateView();


private slots:
	void OnXFoilAdvanced();
    void OnEditPolar();

	void OnCouplePolarGraphs();
	void OnAllPolarGraphs();
	void OnSinglePolarGraph();
	void OnAllPolarGraphsSetting();

	void OnCpGraphSettings();
	void OnPolarFilter();
	void OnInputChanged();
	void OnAnalyze();

	void OnBatchAnalysis();
	void OnMultiThreadedBatchAnalysis();
	void OnCpi();
	void OnCurOppOnly();
	void OnCurveStyle(int index);
	void OnCurveWidth(int index);
	void OnCurveColor();
	void OnDeleteCurFoil();
	void OnDelCurOpp();
	void OnDeleteCurPolar();
	void OnDeleteFoilPolars();
	void OnEditCurPolar();
	void OnExportCurFoil();
	void OnExportCurOpp();
	void OnExportCurPolar();
    void OnExportPolarNREL();
	void OnExportAllPolars();
    void OnExportAllPolarsNREL();
	void OnHideAllOpps();
	void OnHideAllPolars();
	void OnHideFoilPolars();
	void OnImportXFoilPolar();
	void OnInitBL();
	void OnOpPoints();
	void OnPolars();
	void OnResetAllPolarGraphsScales();
	void OnResetFoilScale();
	void OnResetXFoil();
	void OnResetCurPolar();
	void OnSequence();
	void OnShowFoilPolars();
	void OnShowBL();
	void OnShowCurve();
	void OnShowNeutralLine();
	void OnShowPoints();
	void OnShowPanels();
	void OnShowAllOpps();
	void OnShowAllPolars();
	void OnShowPressure();
    void OnNewPolar();
    void OnDeletePolar();
	void OnSpec();
	void OnViscous();
	void OnXDirectStyle();
	void OnGraphSettings();

	void OnShowPolarOpps();
	void OnHidePolarOpps();
	void OnDeletePolarOpps();
	void OnShowFoilOpps();
	void OnHideFoilOpps();
	void OnDeleteFoilOpps();
	void OnDerotateFoil();
	void OnNormalizeFoil();
	void OnCadd();
	void OnPanels();
	void OnFoilCoordinates();
	void OnFoilGeom();
	void OnResetGraphLegend();
	void OnSetTEGap();
	void OnSetLERadius();
	void OnSetFlap();
	void OnInterpolateFoils();
	void OnNacaFoils();
	void OnManageFoils();
	void OnRenamePolar();
	void OnHighlightOpp();

	void OnQGraph();
	void OnCpGraph();
	void OnExportCurXFoilResults();
	void OnCtPlot();
	void OnDbPlot();
	void OnDtPlot();
	void OnRtLPlot();
	void OnRtPlot();
	void OnNPlot();
	void OnCdPlot();
	void OnCfPlot();
	void OnUePlot();
	void OnHPlot();

	void OnOpPointProps();
	void OnPolarProps();
	void OnRenameFoil();


protected:
	void keyPressEvent(QKeyEvent *event);
	void keyReleaseEvent(QKeyEvent *event);
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void wheelEvent (QWheelEvent *event);
	void mouseDoubleClickEvent ( QMouseEvent * event );

	void AddOpData(OpPoint *pNewPoint);
	void SetControls();
	void Connect();
	void CreateOppCurves(OpPoint *pOpp= NULL);
	void CreatePolarCurves();
    void CreateSinglePolarCurve();
	void DeleteFoil(bool bAsk);
	void DeleteOpPoint(bool bCurrent);
	void PaintPolarLegend(QPoint place, int bottom,QPainter &painter);
	void FillComboBoxes(bool bEnable = true);
	void FillPolarCurve(CCurve *pCurve, CPolar *pPolar, int XVar, int YVar);
	void FillOppCurve(OpPoint *pOpp, Graph *pGraph, CCurve *pCurve, bool bInviscid=false);
	void InsertOpPoint(OpPoint *pNewPoint);
	void LoadSettings(QSettings *pSettings);

	void PaintBL(QPainter &painter, OpPoint* pOpPoint, double scale);
	void PaintPressure(QPainter &painter, OpPoint* pOpPoint, double scale);
	void PaintOpPoint(QPainter &painter);
    void PaintEmptyFoilDatabase(QPainter &painter);
    void PaintSingleGraph(QPainter &painter);
	void PaintCoupleGraphs(QPainter &painter);
	void PaintPolarGraphs(QPainter &painter);
	void PaintView(QPainter &painter);
	void ReadParams();
	void SaveSettings(QSettings *pSettings);
	void SetBufferFoil();
	void SetCurveParams();
	void SetFoilScale();
	void SetFoilScale(QRect CltRect);
	void SetHingeMoments(OpPoint *pOpPoint);
	void SetOpPointSequence();
	void SetAnalysisParams();
	void SetGraphTitles(Graph* pGraph, int iX, int iY);
	void SetGraphTitles(Graph* pGraph);
	void SetPolarLegendPos();
	void SetupLayout();	
	void UpdateCurve();


//	bool InitXFoil(CFoil *pFoil=NULL);
//	void InitXFoil2();

	void * GetVariable(CPolar *pPolar, int iVar);
	OpPoint* GetOpPoint(double Alpha);

	QGraph* GetGraph(QPoint &pt);

	CFoil* SetFoil(CFoil* pFoil=NULL);
	CFoil* SetFoil(QString FoilName);
	CPolar *SetPolar(CPolar *pPolar=NULL);
	CPolar *SetPolar(QString PlrName);
	OpPoint *SetOpp(double Alpha=-123456789.0);
    OpPoint *SetOpp(OpPoint *pOpp);


	OpPoint* AddOpPoint(CPolar *pPolar = NULL, bool bStoreOpp=false);

private:
	XFoilAdvancedDlg m_XFoilAdvancedDlg;
	GraphDlg m_GraphDlg;
	BatchDlg m_BatchDlg;
	BatchThreadDlg m_BatchThreadDlg;
	FoilPolarDlg m_FoilPolarDlg;

	TwoDPanelDlg        m_TwoDPanelDlg;
	InterpolateFoilsDlg m_InterpolateFoilsDlg;
	NacaFoilDlg         m_NacaFoilDlg;
	FoilCoordDlg        m_FoilCoordDlg;
	FoilGeomDlg         m_FoilGeomDlg;
	TEGapDlg            m_TEGapDlg;
	LEDlg               m_LEDlg;
	FlapDlg             m_FlapDlg;
	CAddDlg             m_CAddDlg;
	XDirectStyleDlg     m_XDirectStyleDlg;
	ManageFoilsDlg      m_ManageFoilsDlg;
	PolarFilterDlg      m_PolarFilterDlg;

	QStackedWidget *m_pctrlMiddleControls;

	QLabel *m_pctrlUnit1, *m_pctrlUnit2, *m_pctrlUnit3;

	QRadioButton *m_pctrlSpec1;
	QRadioButton *m_pctrlSpec2;
	QRadioButton *m_pctrlSpec3;

	QCheckBox *m_pctrlSequence;
	NumberEdit *m_pctrlAlphaMin  ;
	NumberEdit *m_pctrlAlphaMax ;
	NumberEdit *m_pctrlAlphaDelta;

	QCheckBox *m_pctrlViscous ;
	QCheckBox *m_pctrlInitBL;
    QPushButton *m_pctrlAnalyze, *m_pctrlNewPolar, *m_pctrlDeletePolar, *m_pctrlEditPolar;

	QCheckBox *m_pctrlShowBL, *m_pctrlShowPressure;

	QLabel *m_pctrlPolarProps;

	QCheckBox *m_pctrlShowCurve;
	QCheckBox *m_pctrlShowPoints;
	LineCbBox *m_pctrlCurveStyle;
	LineCbBox *m_pctrlCurveWidth;
	LineButton *m_pctrlCurveColor;

	LineDelegate *m_pStyleDelegate, *m_pWidthDelegate;

	int m_IconWidth, m_IconHeight; //in the comboboxes for curve style & color

	int m_CurveStyle, m_CurveWidth;
	QColor m_CurveColor;

	bool m_bAlpha;			//true if performing an analysis based on aoa, false if based on Cl
	bool m_bViscous;		// true if performing a viscous calculation, false if inviscid
	bool m_bInitBL;			// true if the boundary layer should be initialized for the next xfoil calculation
	bool m_bBL;				// true if the Boundary layer shoud be displayed
	bool m_bPressure;		// true if the pressure distirbution should be displayed
	bool m_bPolar;			// true if the polar view is selected, false if the operating point view is selected
	bool m_bShowUserGraph;	// true if the 5th polar graph should be displayed
	bool m_bShowPanels;		// true if the panels should be displayed on the foil surface
	bool m_bType1, m_bType2, m_bType3, m_bType4; // filter for polar diplay
	bool m_bAutoInitBL;		// true if the BL initialization is left to the code's decision
	bool m_bTrans;			// true if the user is dragging a view
	bool m_bTransGraph;		// true if the user is dragging a graph
	bool m_bFromList;		// true if the batch analysis is based on a list of Re values
	bool m_bFromZero;		// true if the batch analysis should start from Alpha=0
	bool m_bShowTextOutput;	// true if the batch analysis should display text result output
	bool m_bNeutralLine;	// true if the neutral line should be displayed
	bool m_bCurOppOnly;		// true if only the current operating point should be displayed
	bool m_bShowInviscid;	// true if the inviscid results should be displayed
	bool m_bCpGraph;		// true if the Cp graph should be displayed
	bool m_bSequence;		// true if a sequential analysis is to be performed
	bool m_bXPressed, m_bYPressed; //true if the corresponding key is pressed
	bool m_bHighlightOpp;

	int m_XFoilVar;			// defines the variable for current XFoil results
	int m_IterLim;			// max iteratio limit for XFoil

	int m_iPlrGraph;		// defines whch polar graph is selected if m_iPlrView=1
	int m_iPlrView;			// 0 is all, 1 is single, 2 is two !
	int m_FoilYPos;		// y position for the foil display, in pixels from the bottom of the screen

	double m_fFoilScale;	// foil display scale
	double m_ReList[30];	// for batch analysis
	double m_MachList[30];	// for batch analysis
	double m_NCritList[30];	// for batch analysis
	int m_NRe;				// number of Re values in the ReList

	double m_Alpha, m_AlphaMax, m_AlphaDelta;
	double m_Cl, m_ClMax, m_ClDelta;
	double m_Reynolds, m_ReynoldsMax, m_ReynoldsDelta;
	double m_Mach;
	double m_ASpec;
	double m_NCrit, m_XTopTr, m_XBotTr;
	enumPolarType m_PolarType;


	CFoil m_BufferFoil;		// used for screen drawing and temporary geometric mods
    CPolar *m_pCurPolar;	// pointer to the currently selected foil polar
	OpPoint * m_pCurOpp;	// pointer to the currently selected foil operating point

	QGraph* m_pCpGraph;		//pointers to the various graphs
	QGraph* m_pPolarGraph;
	QGraph* m_pCmGraph;
	QGraph* m_pCzGraph;
	QGraph* m_pTrGraph;
	QGraph* m_pUserGraph;
	QGraph* m_pCurGraph;

	QFile m_XFile;		//output file for the analysis


	void *m_pMainFrame;
	void *m_p2DWidget;


	QColor m_crFoilColor, m_crBLColor, m_crPressureColor, m_crNeutralColor; //foil display parameters
	int m_iFoilStyle, m_iFoilWidth;
	int m_iBLStyle, m_iBLWidth;
	int m_iPressureStyle, m_iPressureWidth;
	int m_iNeutralStyle, m_iNeutralWidth;

	QRect m_rCltRect;		// the client area
	QPoint m_PointDown;		// the client point for the last mouse left-click
	QPoint m_FoilOffset;		// the screen offset position for the foil display
	QPoint m_PolarLegendOffset;
//	CFoilAnalysisDlg m_Adlg;	// the polar definition dialog box
//	CUFOListDlg m_UFOdlg;		// the foil management dialog box
	XFoilAnalysisDlg m_XFdlg;	// the dialog class which manages the xfoil calculation and display output

	XFoil *m_pXFoil;		// a pointer to the XFoil object
};


#endif // QXDIRECT_H
