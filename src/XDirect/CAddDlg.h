/****************************************************************************

	Corner Add class
	Copyright (C) 2004-2009 Andre Deperrois adeperrois@xflr5.com

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

#ifndef CADDDLG_H
#define CADDDLG_H

#include <QDialog>
#include <QLabel>
#include <QRadioButton>


#include "../Objects/Foil.h"
#include "../Misc/NumberEdit.h"

class CAddDlg : public QDialog
{
	Q_OBJECT

	friend class MainFrame;
	friend class QXDirect;
	friend class QAFoil;

public:
	CAddDlg();
	void InitDialog();
	void SetupLayout();

private:
	void keyPressEvent(QKeyEvent *event);

private slots:
	void OnApply();
	void OnUniform();


private:
	static void* s_pXFoil;

	QPushButton	*ApplyButton, *OKButton, *CancelButton;
	QLabel *m_pctrlAtPanel;
	QLabel *m_pctrlTotal;
	QLabel *m_pctrlAdded;
	QLabel *m_pctrlMaxAngle;
	QRadioButton	*m_pctrlrb1;
	QRadioButton	*m_pctrlrb2;
	NumberEdit	*m_pctrlTo;
	NumberEdit	*m_pctrlFrom;
	NumberEdit	*m_pctrlAngTol;

	CFoil* m_pMemFoil;
	CFoil* m_pBufferFoil;
	void* m_pXDirect;
	void* m_pAFoil;
	double atol;
	int m_iSplineType;


};

#endif // CADDDLG_H
