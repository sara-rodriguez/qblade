/****************************************************************************

	ObjectPropsDlg Class
	Copyright (C) 2010 Andre Deperrois adeperrois@xflr5.com

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

#ifndef ObjectPropsDlg_H
#define ObjectPropsDlg_H
#include <QTextEdit>
#include <QDialog>
#include "../Objects/Polar.h"

class ObjectPropsDlg : public QDialog
{
	Q_OBJECT

	friend class MainFrame;
	friend class QXDirect;


public:
	ObjectPropsDlg();
	void InitDialog();

private:
	void SetupLayout();

	QTextEdit *m_pctrlDescription;
	CPolar *m_pPolar;
	OpPoint *m_pOpp;

	static void *s_pMainFrame;
	void *m_pXDirect;

};

#endif // ObjectPropsDlg_H
