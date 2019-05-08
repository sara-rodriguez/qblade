/****************************************************************************

	Naca Foil Dlg
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

#ifndef NACAFOILDLG_H
#define NACAFOILDLG_H

#include <QDialog>
#include <QLabel>
#include "../Objects/Foil.h"
#include "../Misc/NumberEdit.h"


class NacaFoilDlg : public QDialog
{
	Q_OBJECT
	friend class QXDirect;
	friend class QAFoil;
	friend class XFoil;


private slots:
	void EditingFinished();
	void OnOK();

public:
	NacaFoilDlg();
	void SetupLayout();
	void GenerateFoil();
	void keyPressEvent(QKeyEvent *event);

	static void *s_pXFoil;
	void *m_pXDirect;
	void *m_pAFoil;

	CFoil *m_pBufferFoil;
	QPushButton *OKButton, *CancelButton;
	NumberEdit *m_pctrlNumber;
	NumberEdit *m_pctrlPanels;
	QLabel * m_pctrlMessage;

	bool m_bGenerated;
	static int s_Digits;
	static int s_Panels;
};

#endif // NACAFOIL_H
