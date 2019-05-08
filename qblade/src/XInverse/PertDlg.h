/****************************************************************************

	PertDlg class
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

#ifndef PERTDLG_H
#define PERTDLG_H

#include <QDialog>
#include <QCheckBox>
#include <QPushButton>
#include <QTableView>
#include <QStandardItem>
#include "../Misc/FloatEditDelegate.h"
#include "../Objects/Foil.h"
#include "../Misc/NumberEdit.h"
#include "../Params.h"

class PertDlg : public QDialog
{
	Q_OBJECT

public:
	PertDlg(void *pParent=NULL);
	friend class MainFrame;
	friend class QXInverse;

private slots:
	void OnCellChanged(QWidget *pWidget);
	void OnRestore();
	void OnApply();
	void OnOK();

private:
	void SetupLayout();
	void InitDialog();
	void FillCnModel() ;
	void ReadData();

private:

	QPushButton *OKButton, *CancelButton, *ApplyButton, *RestoreButton;

	QTableView *m_pctrlCnTable;
	QStandardItemModel *m_pCnModel;
	FloatEditDelegate *m_pFloatDelegate;

//	FloatEdit *m_pctrlCnr;
//	FloatEdit *m_pctrlCni;
//	QListBox *m_pctrlCnList;

protected:
	void keyPressEvent(QKeyEvent *event);



private:
	void * m_pXInverse;
	int   m_nc; 
	double m_cnr[IMX+1];
	double m_cni[IMX+1];
	double m_backr[IMX+1];
	double m_backi[IMX+1];
};
#endif

