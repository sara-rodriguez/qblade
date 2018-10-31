/****************************************************************************

	Naca Foil Dlg
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


#include "NacaFoilDlg.h"
#include "XDirect.h"
#include "../Design/AFoil.h"



void *NacaFoilDlg::s_pXFoil;
int NacaFoilDlg::s_Digits = 0;
int NacaFoilDlg::s_Panels = 100;

NacaFoilDlg::NacaFoilDlg()
{
	setWindowTitle(tr("NACA Foils"));
	m_pAFoil = NULL;
	m_pXDirect = NULL;
	m_bGenerated = false;
	m_pBufferFoil = NULL;

	SetupLayout();
	m_pctrlNumber->setValue(s_Digits);
	m_pctrlPanels->setValue(s_Panels);
}


void NacaFoilDlg::SetupLayout()
{
	QGridLayout *MainGrid = new QGridLayout;
	QLabel *NacaNumber   = new QLabel(tr("4 or 5 digits"));
	QLabel *PanelNumber  = new QLabel(tr("Number of Panels"));

	m_pctrlNumber = new NumberEdit(NumberEdit::Standard, 0);
	m_pctrlPanels = new NumberEdit(NumberEdit::Standard, 0);
	m_pctrlNumber->setValue(100);
	m_pctrlMessage = new QLabel();
	m_pctrlMessage->setMinimumWidth(120);

//	QValidator *PanelValid = new QIntValidator(0, IQX, this);
//	m_pctrlPanels->setValidator(PanelValid);

//	QValidator *NacaValid = new QIntValidator(0, 100000000, this);
//	m_pctrlNumber->setValidator(NacaValid);

	m_pctrlNumber->setAlignment(Qt::AlignRight);
	m_pctrlPanels->setAlignment(Qt::AlignRight);
	MainGrid->addWidget(NacaNumber,     1,1, 1,1, Qt::AlignRight);
	MainGrid->addWidget(m_pctrlNumber,  1,2, 1,1, Qt::AlignRight);
	MainGrid->addWidget(m_pctrlMessage, 2,1, 1,2, Qt::AlignRight);
	MainGrid->addWidget(PanelNumber,    3,1, 1,1, Qt::AlignRight);
	MainGrid->addWidget(m_pctrlPanels,  3,2, 1,1, Qt::AlignRight);


	QHBoxLayout *CommandButtons = new QHBoxLayout;
	OKButton = new QPushButton(tr("OK"));
	OKButton->setAutoDefault(false);
	CancelButton = new QPushButton(tr("Cancel"));
	CancelButton->setAutoDefault(false);
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(OKButton);
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(CancelButton);
	CommandButtons->addStretch(1);
	connect(OKButton, SIGNAL(clicked()),this, SLOT(OnOK()));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(reject()));

	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addLayout(MainGrid);
	mainLayout->addStretch(1);
	mainLayout->addSpacing(30);
	mainLayout->addLayout(CommandButtons);

	setLayout(mainLayout);

	connect(m_pctrlNumber, SIGNAL(editingFinished()), this, SLOT(EditingFinished()));
	connect(m_pctrlPanels, SIGNAL(editingFinished()), this, SLOT(EditingFinished()));
}


void NacaFoilDlg::EditingFinished()
{
	s_Digits = m_pctrlNumber->text().toInt();
	s_Panels = m_pctrlPanels->getValue();

	GenerateFoil();
	OKButton->setFocus();
}


void NacaFoilDlg::GenerateFoil()
{
	int itype;
	QAFoil *pAFoil = (QAFoil*)m_pAFoil;
	QXDirect *pXDirect = (QXDirect*)m_pXDirect;
	XFoil *pXFoil = (XFoil*)s_pXFoil;
	pXFoil->lflap = false;
	pXFoil->lbflap = false;

	if(s_Digits<=25099) itype = 5;
	if(s_Digits<=9999 ) itype = 4;

	if(itype==4) pXFoil->naca4(s_Digits, (int)(s_Panels/2));
	else if(itype==5)
	{
		int three  = s_Digits/100;
		if(three!=210 && three !=220 && three !=230 && three !=240 && three !=250)
		{
			m_pctrlNumber->selectAll();
			m_pctrlMessage->setText(tr("Illegal NACA Number"));
			m_bGenerated = false;
			return;
		}
		if(!pXFoil->naca5(s_Digits, s_Panels))
		{
			m_bGenerated = false;
			m_pctrlMessage->setText(tr("Illegal NACA Number"));
			return;
		}
	}
	else
	{
		m_pctrlNumber->selectAll();
		m_pctrlMessage->setText(tr("Illegal NACA Number"));
		m_bGenerated = false;
		return;
	}
	m_pctrlMessage->setText(" ");

	for (int j=0; j< pXFoil->nb; j++)
	{
		m_pBufferFoil->xb[j] = pXFoil->xb[j+1];
		m_pBufferFoil->yb[j] = pXFoil->yb[j+1];
		m_pBufferFoil->x[j]  = pXFoil->xb[j+1];
		m_pBufferFoil->y[j]  = pXFoil->yb[j+1];
	}
	m_pBufferFoil->nb = pXFoil->nb;
	m_pBufferFoil->n = pXFoil->nb;
	m_pBufferFoil->InitFoil();

	if(pXDirect) pXDirect->UpdateView();
	else if(pAFoil) pAFoil->UpdateView();
	m_bGenerated = true;
}



void NacaFoilDlg::keyPressEvent(QKeyEvent *event)
{
	// Prevent Return Key from closing App
	// Generate the foil instead
	switch (event->key())
	{
		case Qt::Key_Return:
		{
			if(!OKButton->hasFocus() && !CancelButton->hasFocus())
			{
				GenerateFoil();
				if(m_bGenerated) OKButton->setFocus();
				else
				{
					m_pctrlNumber->selectAll();
				}
			}
			else if (OKButton->hasFocus())
			{
				OnOK();
			}
			return;
		}
		case Qt::Key_Escape:
		{
			reject();
			return;
		}
	}
}


void NacaFoilDlg::OnOK()
{
	GenerateFoil();
	accept();
}
