/****************************************************************************

	FlapDlg class
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

#include <QGridLayout>
#include <QHBoxLayout>
#include <QMessageBox>
#include "FlapDlg.h"
#include "XDirect.h"
#include "../Design/AFoil.h"


FlapDlg::FlapDlg(void */*pParent*/)
{
	setWindowTitle(tr("Flap Dlg"));

	m_pAFoil   = NULL;
	m_pXDirect = NULL;
	m_pMemFoil    = NULL;
	m_pBufferFoil = NULL;
	m_bTEFlap     = false;
	m_TEFlapAngle = 0.0;
	m_TEXHinge    = 80.0;//%
	m_TEYHinge    = 50.0;//%
	m_bLEFlap     = false;
	m_LEFlapAngle = 0.0;
	m_LEXHinge    = 20.0;//%
	m_LEYHinge    = 50.0;//%
	m_pXFoil    = NULL;

	m_bModified   = false;
	m_bApplied    = true;

	SetupLayout();

	connect(ApplyButton, SIGNAL(clicked()),this, SLOT(OnApply()));
	connect(OKButton, SIGNAL(clicked()),this, SLOT(OnOK()));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(reject()));

	connect(m_pctrlLEFlapCheck, SIGNAL(stateChanged(int)), this, SLOT(OnLEFlapCheck(int)));
	connect(m_pctrlTEFlapCheck, SIGNAL(stateChanged(int)), this, SLOT(OnTEFlapCheck(int)));

	connect(m_pctrlLEXHinge, SIGNAL(editingFinished()), this, SLOT(OnChanged()));
	connect(m_pctrlLEYHinge, SIGNAL(editingFinished()), this, SLOT(OnChanged()));
	connect(m_pctrlTEXHinge, SIGNAL(editingFinished()), this, SLOT(OnChanged()));
	connect(m_pctrlTEYHinge, SIGNAL(editingFinished()), this, SLOT(OnChanged()));
	connect(m_pctrlLEFlapAngle, SIGNAL(editingFinished()), this, SLOT(OnChanged()));
	connect(m_pctrlTEFlapAngle, SIGNAL(editingFinished()), this, SLOT(OnChanged()));

}


void FlapDlg::SetupLayout()
{
	QGridLayout *FlapData = new QGridLayout;
	m_pctrlLEFlapCheck = new QCheckBox(tr("L.E. Flap"));
	m_pctrlTEFlapCheck = new QCheckBox(tr("T.E. Flap"));
	m_pctrlLEXHinge    = new NumberEdit;
	m_pctrlLEYHinge    = new NumberEdit;
	m_pctrlTEXHinge    = new NumberEdit;
	m_pctrlTEYHinge    = new NumberEdit;
	m_pctrlTEFlapAngle = new NumberEdit;
	m_pctrlLEFlapAngle = new NumberEdit;

	QLabel *lab1 = new QLabel(tr("Flap Angle"));
	QLabel *lab2 = new QLabel(QString::fromUtf8("deg (")+tr("+ is down") +")");
	QLabel *lab3 = new QLabel(tr("Hinge X Position"));
	QLabel *lab4 = new QLabel(tr("% Chord"));
	QLabel *lab5 = new QLabel(tr("Hinge Y Position"));
	QLabel *lab6 = new QLabel(tr("% Thickness"));

	FlapData->addWidget(m_pctrlLEFlapCheck, 1, 2);
	FlapData->addWidget(m_pctrlTEFlapCheck, 1, 3);
	FlapData->addWidget(lab1, 2, 1);
	FlapData->addWidget(m_pctrlLEFlapAngle, 2, 2);
	FlapData->addWidget(m_pctrlTEFlapAngle, 2, 3);
	FlapData->addWidget(lab2, 2, 4);
	FlapData->addWidget(lab3, 3, 1);
	FlapData->addWidget(m_pctrlLEXHinge, 3, 2);
	FlapData->addWidget(m_pctrlTEXHinge, 3, 3);
	FlapData->addWidget(lab4, 3, 4);
	FlapData->addWidget(lab5, 4, 1);
	FlapData->addWidget(m_pctrlLEYHinge, 4, 2);
	FlapData->addWidget(m_pctrlTEYHinge, 4, 3);
	FlapData->addWidget(lab6, 4, 4);

	QHBoxLayout *CommandButtons = new QHBoxLayout;
	OKButton      = new QPushButton(tr("OK"));
	CancelButton  = new QPushButton(tr("Cancel"));
	ApplyButton  = new QPushButton(tr("Apply"));

	CommandButtons->addStretch(1);
	CommandButtons->addWidget(ApplyButton);
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(OKButton);
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(CancelButton);
	CommandButtons->addStretch(1);

	QVBoxLayout *MainLayout = new QVBoxLayout;
	MainLayout->addLayout(FlapData);
	MainLayout->addLayout(CommandButtons);
	setLayout(MainLayout);
}


void FlapDlg::InitDialog()
{
	m_pctrlTEFlapCheck->setChecked(m_pMemFoil->m_bTEFlap);

	EnableTEFlap(m_pMemFoil->m_bTEFlap);
	m_pctrlTEFlapAngle->setValue(m_pMemFoil->m_TEFlapAngle);
	m_pctrlTEXHinge->setValue(m_pMemFoil->m_TEXHinge);
	m_pctrlTEYHinge->setValue(m_pMemFoil->m_TEYHinge);

	m_pctrlLEFlapCheck->setChecked(m_pMemFoil->m_bLEFlap);
	EnableLEFlap(m_pMemFoil->m_bLEFlap);
	m_pctrlLEFlapAngle->setValue(m_pMemFoil->m_LEFlapAngle);
	m_pctrlLEXHinge->setValue(m_pMemFoil->m_LEXHinge);
	m_pctrlLEYHinge->setValue(m_pMemFoil->m_LEYHinge);
}


void FlapDlg::ReadParams()
{
	m_bLEFlap = m_pctrlLEFlapCheck->isChecked();
	m_LEFlapAngle = m_pctrlLEFlapAngle->getValue();
	m_LEXHinge    = m_pctrlLEXHinge->getValue();
	m_LEYHinge    = m_pctrlLEYHinge->getValue();

	m_bTEFlap = m_pctrlTEFlapCheck->isChecked();
	m_TEFlapAngle = m_pctrlTEFlapAngle->getValue();
	m_TEXHinge    = m_pctrlTEXHinge->getValue();
	m_TEYHinge    = m_pctrlTEYHinge->getValue();

	if(m_LEXHinge>=m_TEXHinge && m_bLEFlap && m_bTEFlap)
	{
		QMessageBox::information(window(), tr("Warning"), tr("The trailing edge hinge must be downstream of the leading edge hinge"));
		m_pctrlLEXHinge->setFocus();
		m_pctrlLEXHinge->selectAll();
	}
}


void FlapDlg::OnApply()
{
	QXDirect *pXDirect = (QXDirect*)m_pXDirect;
	QAFoil *pAFoil = (QAFoil*)m_pAFoil;

	if(m_bApplied) return;
	//reset everything and retry

	ReadParams();

	m_pBufferFoil->SetTEFlapData(m_bTEFlap, m_TEXHinge, m_TEYHinge, m_TEFlapAngle);
	m_pBufferFoil->SetLEFlapData(m_bLEFlap, m_LEXHinge, m_LEYHinge, m_LEFlapAngle);
	m_pBufferFoil->SetFlap();

	m_bApplied = true;
	m_bModified = true;
	if(pXDirect)    pXDirect->UpdateView();
	else if(pAFoil) pAFoil->UpdateView();
}


void FlapDlg::keyPressEvent(QKeyEvent *event)
{
	switch (event->key())
	{
		case Qt::Key_Escape:
		{
			done(0);
			return;
		}
		case Qt::Key_Return:
		{
			if(!OKButton->hasFocus() && !CancelButton->hasFocus())
			{
				OnApply();
				OKButton->setFocus();
				m_bApplied  = true;
			}
			else
			{
				QDialog::accept();
			}
			break;
		}
		default:
			event->ignore();
			break;
	}
}

void FlapDlg::EnableLEFlap(bool bEnable)
{
	m_pctrlLEFlapAngle->setEnabled(bEnable);
	m_pctrlLEXHinge->setEnabled(bEnable);
	m_pctrlLEYHinge->setEnabled(bEnable);
}

void FlapDlg::EnableTEFlap(bool bEnable)
{
	m_pctrlTEFlapAngle->setEnabled(bEnable);
	m_pctrlTEXHinge->setEnabled(bEnable);
	m_pctrlTEYHinge->setEnabled(bEnable);
}

void FlapDlg::OnTEFlapCheck(int /*state*/)
{
	if(m_pctrlTEFlapCheck->isChecked())
	{
		EnableTEFlap(true);
		m_pctrlTEFlapAngle->setFocus();
	}
	else
		EnableTEFlap(false);
	m_bApplied = false;
	OnApply();
}

void FlapDlg::OnLEFlapCheck(int /*state*/)
{
	if(m_pctrlLEFlapCheck->isChecked())
	{
		EnableLEFlap(true);
		m_pctrlLEFlapAngle->setFocus();
	}
	else
		EnableLEFlap(false);
	m_bApplied = false;
	OnApply();
}

void FlapDlg::OnChanged()
{
	m_bApplied = false;
}

void FlapDlg::OnOK()
{
	if(!m_bApplied)
	{
		OnApply();
		accept();
	}
	else done(1);
}
