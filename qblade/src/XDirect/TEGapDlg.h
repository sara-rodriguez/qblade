/****************************************************************************

	TEGapDlg Class
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


#ifndef TEGAPDlg_H
#define TEGAPDlg_H

#include <QDialog>
#include "../Misc/NumberEdit.h"
#include "../Objects/Foil.h"


class TEGapDlg : public QDialog
{
	Q_OBJECT
	friend class QXDirect;
	friend class QAFoil;

public:
	TEGapDlg();
	void SetupLayout();
	void InitDialog();

private slots:
	void OnChanged();
	void OnOK();
	void OnApply();

private:
	void keyPressEvent(QKeyEvent *event);

public:
	QPushButton *OKButton, *CancelButton, *ApplyButton;
	NumberEdit *m_pctrlBlend, *m_pctrlGap;

	bool m_bApplied, m_bModified;
	double m_Gap, m_Blend;

	void *m_pXDirect;
	void *m_pAFoil;

	CFoil* m_pBufferFoil;
	CFoil* m_pMemFoil;
	static void* s_pXFoil;
};

#endif // TEGAPDLG_H
