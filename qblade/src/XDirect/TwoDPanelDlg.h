/****************************************************************************

	TwoDPanelDlg Class
	Copyright (C) 2008-2008 Andre Deperrois adeperrois@xflr5.com

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


#ifndef TWODPANELDLG_H
#define TWODPANELDLG_H

#include <QDialog>
#include "../Misc/NumberEdit.h"


class TwoDPanelDlg : public QDialog
{
	Q_OBJECT
	friend class QAFoil;
	friend class QXDirect;

private slots:
	void OnApply();
	void OnOK();
	void OnChanged();

public:
	TwoDPanelDlg();

	static void *s_pXFoil;

	void InitDialog();

private:
	void keyPressEvent(QKeyEvent *event);
	void SetupLayout();
	void ReadParams();

	QPushButton *OKButton, *CancelButton, *ApplyButton;

	QLineEdit  *m_pctrlNPanels;
	NumberEdit *m_pctrlCVpar,  *m_pctrlCTErat, *m_pctrlCTRrat;
	NumberEdit *m_pctrlXsRef1, *m_pctrlXsRef2, *m_pctrlXpRef1, *m_pctrlXpRef2;

	bool m_bApplied;
	bool m_bModified;

	int npan;
	double cvpar;
	double cterat;
	double ctrrat;
	double xsref1;
	double xsref2;
	double xpref1;
	double xpref2;

	void *m_pBufferFoil;
	void *m_pMemFoil;
	void *m_pXDirect;
	void *m_pAFoil;
};

#endif // TWODPANELDLG_H
