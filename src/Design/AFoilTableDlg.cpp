/****************************************************************************

	AFoilTableDlg Class
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

#include "AFoil.h"
#include "AFoilTableDlg.h"
#include <QGridLayout>
#include <QVBoxLayout>
#include <QLabel>
#include "../Misc/LinePickerDlg.h"
#include <QKeyEvent>

AFoilTableDlg::AFoilTableDlg(void */*pParent*/)
{
	setWindowTitle(tr("Foil Table Columns"));

	SetupLayout();

	m_bFoilName = m_bPoints = true;
	m_bThickness = m_bThicknessAt = m_bCamber = m_bCamberAt = true;
	m_bTEFlapAngle = m_bTEXHinge = m_bTEYHinge = m_bLEFlapAngle = m_bLEXHinge = m_bLEYHinge = true;

}

void AFoilTableDlg::keyPressEvent(QKeyEvent *event)
{
	switch (event->key())
	{
		case Qt::Key_Escape:
		{
			done(0);
			break;
		}
		case Qt::Key_Return:
		{
			if(!OKButton->hasFocus() && !CancelButton->hasFocus())
			{
				OKButton->setFocus();
//				m_bApplied  = true;
			}
			else
			{
				QDialog::accept();
			}
			break;
		}
		default:
			event->ignore();
	}
}


void AFoilTableDlg::InitDialog()
{
	m_pctrlFoilName->setChecked(m_bFoilName);
	m_pctrlThickness->setChecked(m_bThickness);
	m_pctrlThicknessAt->setChecked(m_bThicknessAt);
	m_pctrlCamber->setChecked(m_bCamber);
	m_pctrlCamberAt->setChecked(m_bCamberAt);
	m_pctrlPoints->setChecked(m_bPoints);
	m_pctrlTEFlapAngle->setChecked(m_bTEFlapAngle);
	m_pctrlTEXHinge->setChecked(m_bTEXHinge);
	m_pctrlTEYHinge->setChecked(m_bTEYHinge);
	m_pctrlLEFlapAngle->setChecked(m_bLEFlapAngle);
	m_pctrlLEXHinge->setChecked(m_bLEXHinge);
	m_pctrlLEYHinge->setChecked(m_bLEYHinge);
}


void AFoilTableDlg::SetupLayout()
{
	QVBoxLayout *ColumnsLayout = new QVBoxLayout;
	m_pctrlFoilName    = new QCheckBox(tr("Foil Name"));
	m_pctrlThickness   = new QCheckBox(tr("Thickness"));
	m_pctrlThicknessAt = new QCheckBox(tr("Thickness max. position"));
	m_pctrlCamber      = new QCheckBox(tr("Camber"));
	m_pctrlCamberAt    = new QCheckBox(tr("Camber max. position"));
	m_pctrlPoints      = new QCheckBox(tr("Number of points"));
	m_pctrlTEFlapAngle = new QCheckBox(tr("Trailing edge flap angle"));
	m_pctrlTEXHinge    = new QCheckBox(tr("Trailing edge hinge x-position"));
	m_pctrlTEYHinge    = new QCheckBox(tr("Trailing edge hinge y-position"));
	m_pctrlLEFlapAngle = new QCheckBox(tr("Leading edge flap angle"));
	m_pctrlLEXHinge    = new QCheckBox(tr("Leading edge hinge x-position"));
	m_pctrlLEYHinge    = new QCheckBox(tr("Leading edge hinge y-position"));

	ColumnsLayout->addWidget(m_pctrlFoilName);
	ColumnsLayout->addWidget(m_pctrlThickness);
	ColumnsLayout->addWidget(m_pctrlThicknessAt);
	ColumnsLayout->addWidget(m_pctrlCamber);
	ColumnsLayout->addWidget(m_pctrlCamberAt);
	ColumnsLayout->addWidget(m_pctrlPoints);
	ColumnsLayout->addWidget(m_pctrlTEFlapAngle);
	ColumnsLayout->addWidget(m_pctrlTEXHinge);
	ColumnsLayout->addWidget(m_pctrlTEYHinge);
	ColumnsLayout->addWidget(m_pctrlLEFlapAngle);
	ColumnsLayout->addWidget(m_pctrlLEXHinge);
	ColumnsLayout->addWidget(m_pctrlLEYHinge);

	QHBoxLayout *CommandButtons = new QHBoxLayout;
	OKButton      = new QPushButton(tr("OK"));
	CancelButton  = new QPushButton(tr("Cancel"));
	CommandButtons->addWidget(OKButton);
	CommandButtons->addWidget(CancelButton);

	QVBoxLayout *MainLayout = new QVBoxLayout;
	MainLayout->addLayout(ColumnsLayout);
	MainLayout->addStretch(1);
	MainLayout->addLayout(CommandButtons);
	setLayout(MainLayout);

	connect(OKButton, SIGNAL(clicked()),this, SLOT(OnOK()));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}

void AFoilTableDlg::OnOK()
{
	m_bFoilName    = m_pctrlFoilName->isChecked();
	m_bThickness   = m_pctrlThickness->isChecked();
	m_bThicknessAt = m_pctrlThicknessAt->isChecked();
	m_bCamber      = m_pctrlCamber->isChecked();
	m_bCamberAt    = m_pctrlCamberAt->isChecked();
	m_bPoints      = m_pctrlPoints->isChecked();
	m_bTEFlapAngle = m_pctrlTEFlapAngle->isChecked();
	m_bTEXHinge    = m_pctrlTEXHinge->isChecked();
	m_bTEYHinge    = m_pctrlTEYHinge->isChecked();
	m_bLEFlapAngle = m_pctrlLEFlapAngle->isChecked();
	m_bLEXHinge    = m_pctrlLEXHinge->isChecked();
	m_bLEYHinge    = m_pctrlLEYHinge->isChecked();

	accept();
}
