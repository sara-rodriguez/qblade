/****************************************************************************

	ManageFoilsDlg Class
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


#ifndef MANAGEFOILSDLG_H
#define MANAGEFOILSDLG_H

#include <QDialog>

#include <QPushButton>
#include <QTableView>
#include <QStandardItemModel>
#include "../Design/FoilTableDelegate.h"
#include "../Objects/Foil.h"


class ManageFoilsDlg : public QDialog
{
	Q_OBJECT
	friend class QXDirect;
	friend class MainFrame;

public:
	ManageFoilsDlg();
	void InitDialog(QString FoilName);

private slots:
	void OnDelete();
	void OnRename();
	void OnExport();
	void OnFoilClicked(const QModelIndex& index);
	void OnDoubleClickTable(const QModelIndex &index);

private:
	void resizeEvent(QResizeEvent *event);
	void keyPressEvent(QKeyEvent *event);

	void FillFoilTable();
	void FillTableRow(int row);

	void SetupLayout();

private:
	QPushButton *CloseButton;
	QPushButton *m_pctrlRename, *m_pctrlDelete, *m_pctrlSelect, *m_pctrlExport;
	QTableView *m_pctrlFoilTable;
	QStandardItemModel *m_pFoilModel;
	FoilTableDelegate *m_pFoilDelegate;

	int m_iSelection;
	CFoil *m_pFoil;
	void *m_pMainFrame, *m_pXDirect;
};

#endif // MANAGEFOILSDLG_H
