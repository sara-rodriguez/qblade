/****************************************************************************

	UnitsDlg Class
	Copyright (C) 2009 Andre Deperrois adeperrois@xflr5.com

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


#ifndef UNITSDLG_H
#define UNITSDLG_H

#include <QComboBox>
#include <QDialog>
#include <QPushButton>
#include <QLabel>

class UnitsDlg : public QDialog
{
	Q_OBJECT
	friend class MainFrame;

public:
    UnitsDlg();

private slots:
	void OnSelChanged(const QString &strong);

private:
	QPushButton *OKButton, *CancelButton;
	QComboBox	*m_pctrlMoment;
	QComboBox	*m_pctrlSurface;
	QComboBox	*m_pctrlWeight;
	QComboBox	*m_pctrlSpeed;
	QComboBox	*m_pctrlLength;
	QComboBox	*m_pctrlForce;
	QLabel *m_pctrlForceFactor, *m_pctrlForceInvFactor;
	QLabel *m_pctrlLengthFactor, *m_pctrlLengthInvFactor;
	QLabel *m_pctrlSpeedFactor, *m_pctrlSpeedInvFactor;
	QLabel *m_pctrlSurfaceFactor, *m_pctrlSurfaceInvFactor;
	QLabel *m_pctrlWeightFactor, *m_pctrlWeightInvFactor;
	QLabel *m_pctrlMomentFactor, *m_pctrlMomentInvFactor;
	QLabel *m_pctrlQuestion;
        ///////////new code DM//////////
        QLabel *m_pctrlPowerFactor, *m_pctrlPowerInvFactor;
        QComboBox *m_pctrlPower;
        int m_Power;
        double m_WtoUnit;
        ///////end new code DM//////

private:
	void InitDialog();
	void SetupLayout();

	bool m_bLengthOnly;
	int m_Length, m_Area, m_Weight, m_Speed, m_Force, m_Moment;
	double m_mtoUnit;
	double m_kgtoUnit;
	double m_mstoUnit;
	double m_NtoUnit;
	double m_NmtoUnit;
	double m_m2toUnit;

	QString m_Question;
};

#endif // UNITSDLG_H
