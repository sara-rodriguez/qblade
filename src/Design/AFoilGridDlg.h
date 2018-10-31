/****************************************************************************

	AFoilGridDlg Class
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

#ifndef AFOILGRIDDLG_H
#define AFOILGRIDDLG_H

#include <QDialog>
#include <QPushButton>
#include <QCheckBox>
#include "../Misc/LineButton.h"
#include "../Misc/NumberEdit.h"


class AFoilGridDlg : public QDialog
{
	Q_OBJECT

	friend class QAFoil;

public:
	AFoilGridDlg(void *pParent = NULL);
	void InitDialog();

private slots:
	void OnApply();
	void OnOK();
	void OnScale();
	void OnNeutralStyle();
	void OnXMajStyle();
	void OnXMinStyle();
	void OnYMajStyle();
	void OnYMinStyle();
	void OnNeutralShow(bool bShow);
	void OnXMajShow(bool bShow);
	void OnYMajShow(bool bShow);
	void OnXMinShow(bool bShow);
	void OnYMinShow(bool bShow);

private:
	void SetupLayout();
	void keyPressEvent(QKeyEvent *event);

	void *m_pAFoil;
	QCheckBox  *m_pctrlNeutralShow, *m_pctrlScale, *m_pctrlXMajShow, *m_pctrlYMajShow, *m_pctrlXMinShow, *m_pctrlYMinShow;
	LineButton *m_pctrlNeutralStyle, *m_pctrlXMajStyle, *m_pctrlYMajStyle, *m_pctrlXMinStyle, *m_pctrlYMinStyle;
	NumberEdit *m_pctrlXUnit, *m_pctrlYUnit,*m_pctrlXMinUnit, *m_pctrlYMinUnit;
	QPushButton	*ApplyButton, *OKButton, *CancelButton;

	bool m_bNeutralLine, m_bScale;
	bool m_bXGrid,m_bYGrid;
	bool m_bXMinGrid, m_bYMinGrid;
	int m_XStyle, m_YStyle;
	int m_XWidth, m_YWidth;
	int m_XMinStyle, m_YMinStyle;
	int m_XMinWidth, m_YMinWidth;
	int m_NeutralStyle, m_NeutralWidth;
	double m_XUnit, m_YUnit;
	double m_XMinUnit, m_YMinUnit;
	QColor m_XColor,m_YColor;
	QColor m_XMinColor,m_YMinColor;
	QColor m_NeutralColor;

};

#endif // AFOILGRIDDLG_H
