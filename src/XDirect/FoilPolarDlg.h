/****************************************************************************

	FoilPolarDlg Class
	Copyright (C) 2008 Andre Deperrois adeperrois@xflr5.com

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


#ifndef FOILPOLARDLG_H
#define FOILPOLARDLG_H

#include <QDialog>
#include <QRadioButton>
#include <QLabel>
#include <QString>
#include "../Misc/NumberEdit.h"
#include "../Objects/Polar.h"


class FoilPolarDlg : public QDialog
{
    friend class QOperWidget;

    Q_OBJECT

public:
	FoilPolarDlg(void *pParent=NULL);

	void ReadParams();
	void keyPressEvent(QKeyEvent *event);
	void InitDialog();
	void SetPlrName();
	void SetupLayout();
	void SetDensity();

	static void *s_pXDirect;
	static void *s_pMainFrame;

	QRadioButton *m_pctrlAuto1;
	QRadioButton *m_pctrlAuto2;

	QLabel *m_pctrlReLabel;
	QLabel *m_pctrlMachLabel;

	QLineEdit *m_pctrlAnalysisName;
	QRadioButton *m_rbtype1;
	QRadioButton *m_rbtype2;
	QRadioButton *m_rbtype3;
	QRadioButton *m_rbtype4;

	NumberEdit *m_pctrlReynolds;
	NumberEdit *m_pctrlMach;

	NumberEdit *m_pctrlChord, *m_pctrlMass, *m_pctrlSpan;
	QLabel *m_pctrlLengthUnit1, *m_pctrlLengthUnit2, *m_pctrlMassUnit;

	QRadioButton *m_pctrlUnit1, *m_pctrlUnit2;
	QLabel *m_pctrlRho, *m_pctrlNu, *m_pctrlViscosityUnit, *m_pctrlDensityUnit;
	NumberEdit *m_pctrlDensity, *m_pctrlViscosity;

	QPushButton *OKButton;
	QPushButton *CancelButton;

	NumberEdit *m_pctrlNCrit;
	NumberEdit *m_pctrlTopTrans;
	NumberEdit *m_pctrlBotTrans;

	bool  m_bAutoName;
	enumPolarType m_PolarType;
	int m_MaTypDef, m_ReTypDef;
	int m_UnitType;
	double m_Reynolds;
	double m_Viscosity, m_Density;
	double m_Chord, m_Span, m_Mass;
	double m_Mach;
	double m_ReDef;
	double m_ASpec;
	double m_XTopTr, m_XBotTr;
	double m_NCrit;
	QString m_FoilName;
	QString m_PlrName;

public slots:
	void OnAutoName();
	void OnOK();
	void OnPolarType();
	void OnNameChanged();
	void EditingFinished();
	void OnUnit();
	void OnEditingFinished();
};

#endif // FOILPOLARDLG_H
