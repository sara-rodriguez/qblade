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

#include <QMessageBox>

#include "CAddDlg.h"
#include "XDirect.h"
#include "../Design/AFoil.h"
#include "XFoil.h"



void *CAddDlg::s_pXFoil;

CAddDlg::CAddDlg()
{
	setWindowTitle(tr("Local Panel Refinement"));
	m_pAFoil      = NULL;
	m_pXDirect    = NULL;
	m_pBufferFoil = NULL;
	m_pMemFoil    = NULL;

	m_iSplineType = 0;
	SetupLayout();
	connect(ApplyButton, SIGNAL(clicked()),this, SLOT(OnApply()));
	connect(OKButton, SIGNAL(clicked()),this, SLOT(accept()));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}

void CAddDlg::SetupLayout()
{
	QGridLayout *RefineGrid =new QGridLayout;
	QLabel *lab1 = new QLabel(tr("Angle Criterion ")+QString::fromUtf8("(deg)"));
	QLabel *lab2 = new QLabel(tr("Type of Spline"));
	QLabel *lab3 = new QLabel(tr("Refinement X Limits"));
	QLabel *lab4 = new QLabel(tr("From"));
	QLabel *lab5 = new QLabel(tr("To"));
	lab4->setAlignment(Qt::AlignCenter);
	lab5->setAlignment(Qt::AlignCenter);
	m_pctrlAngTol = new NumberEdit;
	m_pctrlFrom   = new NumberEdit;
	m_pctrlTo     = new NumberEdit;

	m_pctrlrb1 = new QRadioButton(tr("Uniform"));
	m_pctrlrb2 = new QRadioButton(tr("Arc Length"));

	RefineGrid->addWidget(lab1,1,1);
	RefineGrid->addWidget(lab2,2,1);
	RefineGrid->addWidget(lab3,5,1);
	RefineGrid->addWidget(m_pctrlAngTol,1,2);
	RefineGrid->addWidget(m_pctrlrb1,2,2);
	RefineGrid->addWidget(m_pctrlrb2,2,3);
	RefineGrid->addWidget(lab4, 4, 2);
	RefineGrid->addWidget(lab5, 4, 3);
	RefineGrid->addWidget(m_pctrlFrom,5,2);
	RefineGrid->addWidget(m_pctrlTo,5,3);

	m_pctrlTotal    = new QLabel(tr("Total"));
	m_pctrlAdded    = new QLabel(tr("Added"));
	m_pctrlMaxAngle = new QLabel(tr("MaxAngle"));
	m_pctrlAtPanel  = new QLabel(tr("At Panel"));

	QHBoxLayout *CommandButtons = new QHBoxLayout;
	OKButton      = new QPushButton(tr("Accept"));
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
	MainLayout->addLayout(RefineGrid);
	MainLayout->addStretch(1);
	MainLayout->addWidget(m_pctrlTotal);
	MainLayout->addWidget(m_pctrlAdded);
	MainLayout->addWidget(m_pctrlMaxAngle);
	MainLayout->addWidget(m_pctrlAtPanel);
	MainLayout->addStretch(1);
	MainLayout->addLayout(CommandButtons);
	setLayout(MainLayout);
	setMinimumHeight(300);
}


void CAddDlg::OnApply()
{
	XFoil *pXFoil = (XFoil*)s_pXFoil;
	QAFoil * pAFoil = (QAFoil*)m_pAFoil;
	QXDirect *pXDirect = (QXDirect*)m_pXDirect;
	int i;

	for (i=0; i< m_pMemFoil->nb; i++)
	{
		pXFoil->xb[i+1] = m_pMemFoil->xb[i] ;
		pXFoil->yb[i+1] = m_pMemFoil->yb[i];
	}
	pXFoil->nb = m_pMemFoil->nb;

	pXFoil->lflap = false;
	pXFoil->lbflap = false;

	if(pXFoil->Preprocess())
	{
		pXFoil->CheckAngles();
/*		for (int k=0; k<pXFoil->n;k++){
			m_pMemFoil->nx[k] = pXFoil->nx[k+1];
			m_pMemFoil->ny[k] = pXFoil->ny[k+1];
		}
		m_pMemFoil->n = pXFoil->n;*/
	}
	else
	{
		QMessageBox::information(window(), tr("Warning"), tr("Unrecognized foil format"));
		return;
	}

	if (m_pctrlrb1->isChecked())
		m_iSplineType = 1;
	else
		m_iSplineType = 2;

	int added = pXFoil->cadd(m_iSplineType, m_pctrlAngTol->getValue(),
						m_pctrlFrom->getValue(), m_pctrlTo->getValue());
	pXFoil->abcopy();

	QString strong;
	strong  =QString(tr("Total number of points is %1")).arg(pXFoil->n);
	m_pctrlTotal->setText(strong);
	strong = QString(tr("(added %1 points to original foil)")).arg(added);
	m_pctrlAdded->setText(strong);

	for (i=0; i< pXFoil->n; i++)
	{
		m_pBufferFoil->xb[i] = pXFoil->x[i+1];
		m_pBufferFoil->yb[i] = pXFoil->y[i+1];
	}
	m_pBufferFoil->nb = pXFoil->n;
//		pXFoil->SetFoilFlap(m_pBufferFoil);
	m_pBufferFoil->InitFoil();
	m_pBufferFoil->SetFlap();

	pXFoil->CheckAngles();
	strong = QString(tr("Maximum panel angle is %1")).arg( pXFoil->amax,0,'f',1);
	m_pctrlMaxAngle->setText(strong);
	strong = QString(tr("at panel position %1")).arg(pXFoil->imax);
	m_pctrlAtPanel->setText(strong);

	if(pXDirect)    pXDirect->UpdateView();
	else if(pAFoil) pAFoil->UpdateView();
}


void CAddDlg::OnUniform()
{
	if(m_pctrlrb1->isChecked()) m_iSplineType = 1;
	else                      m_iSplineType = 2;
}

void CAddDlg::InitDialog()
{
	XFoil *pXFoil = (XFoil*)s_pXFoil;

	double xbmin, xbmax, xrf1, xrf2;
//	int ispl;

	xbmin = pXFoil->xb[1];
	xbmax = pXFoil->xb[1];

	for( int i=1; i<= pXFoil->nb; i++)
	{
		xbmin = qMin(xbmin, pXFoil->xb[i]);
		xbmax = qMax(xbmax, pXFoil->xb[i]);
	}

	//----- default inputs
	atol = 0.5 * pXFoil->amax;
//	ispl = 1;
	xrf1 = xbmin - 0.1*(xbmax-xbmin);
	xrf2 = xbmax + 0.1*(xbmax-xbmin);

	m_pctrlrb1->setChecked(1);
	m_iSplineType = 1;
	m_pctrlFrom->setValue(xrf1);
	m_pctrlTo->setValue(xrf2);
	m_pctrlAngTol->setValue(atol);

	QString strong;
	pXFoil->CheckAngles();
	strong = QString(tr("Maximum panel angle is %1 deg")).arg(pXFoil->amax,0,'f',1);
	m_pctrlMaxAngle->setText(strong);
	strong = QString(tr("at panel position %1")).arg(pXFoil->imax);
	m_pctrlAtPanel->setText(strong);
	m_pctrlAdded->setText("  ");
	strong = QString(tr("Total number of points is %1")).arg(pXFoil->n);
	m_pctrlTotal->setText(strong);

}


void CAddDlg::keyPressEvent(QKeyEvent *event)
{
	// Prevent Return Key from closing dialog
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

