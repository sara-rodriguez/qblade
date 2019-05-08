/****************************************************************************

	InterpolateFoilsDlg Class
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

#include <QGroupBox>

#include "InterpolateFoilsDlg.h"
#include "../MainFrame.h"
#include "XFoil.h"
#include "XDirect.h"
#include "../Design/AFoil.h"
#include "../Store.h"

void *InterpolateFoilsDlg::s_pXFoil;


InterpolateFoilsDlg::InterpolateFoilsDlg()
{
	setWindowTitle(tr("Interpolate Foils"));
	m_pXDirect = NULL;
	m_pAFoil = NULL;
	m_pMainFrame = NULL;
	m_pBufferFoil = NULL;

	SetupLayout();

	connect(m_pctrlFoil1,  SIGNAL(activated(int)),    this, SLOT(OnSelChangeFoil1(int)));
	connect(m_pctrlFoil2,  SIGNAL(activated(int)),    this, SLOT(OnSelChangeFoil2(int)));
	connect(m_pctrlFrac,   SIGNAL(editingFinished()), this, SLOT(OnFrac()));
	connect(m_pctrlSlider, SIGNAL(sliderMoved(int)),  this, SLOT(OnVScroll(int)));
}


void InterpolateFoilsDlg::SetupLayout()
{
	QVBoxLayout *LeftSide = new QVBoxLayout;
	m_pctrlFoil1 = new QComboBox;
	m_pctrlFoil2 = new QComboBox;
	m_pctrlCamb1 = new QLabel(tr("Camb1"));
	m_pctrlCamb2 = new QLabel(tr("Camb2"));
	m_pctrlCamb3 = new QLabel(tr("Camb3"));
	m_pctrlThick1 = new QLabel(tr("Thick1"));
	m_pctrlThick2 = new QLabel(tr("Thick2"));
	m_pctrlThick3 = new QLabel(tr("Thick3"));
	m_pctrlCamb1->setMinimumWidth(250);
	m_pctrlCamb2->setMinimumWidth(250);
	m_pctrlCamb3->setMinimumWidth(250);
	m_pctrlThick1->setMinimumWidth(250);
	m_pctrlThick2->setMinimumWidth(250);
	m_pctrlThick3->setMinimumWidth(250);
	LeftSide->addWidget(m_pctrlFoil1);
	LeftSide->addWidget(m_pctrlCamb1);
	LeftSide->addWidget(m_pctrlThick1);
	LeftSide->addStretch(1);
	LeftSide->addWidget(m_pctrlFoil2);
	LeftSide->addWidget(m_pctrlCamb2);
	LeftSide->addWidget(m_pctrlThick2);
	LeftSide->addStretch(2);
	QVBoxLayout *Foil3 = new QVBoxLayout;
	m_pctrlNewFoilName = new QLineEdit(tr("New Foil Name"));
	Foil3->addWidget(m_pctrlNewFoilName);
	Foil3->addWidget(m_pctrlCamb3);
	Foil3->addWidget(m_pctrlThick3);
	QGroupBox *Foil3Group = new QGroupBox(tr("Interpolated Foil"));
	Foil3Group->setLayout(Foil3);
	LeftSide->addWidget(Foil3Group);

	QVBoxLayout *RightSide = new QVBoxLayout;
	m_pctrlSlider = new QSlider;
	m_pctrlFrac = new NumberEdit;
	m_pctrlSlider->setMinimumHeight(300);
	m_pctrlSlider->setMinimum(0);
	m_pctrlSlider->setMaximum(100);
	m_pctrlSlider->setSliderPosition(100);
	m_pctrlSlider->setTickInterval(10);
	m_pctrlSlider->setTickPosition(QSlider::TicksLeft);
	RightSide->addWidget(m_pctrlSlider);
	RightSide->addWidget(m_pctrlFrac);
	RightSide->addStretch(1);

	QHBoxLayout *CommandButtons = new QHBoxLayout;
	OKButton = new QPushButton(tr("OK"));
	CancelButton = new QPushButton(tr("Cancel"));
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(OKButton);
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(CancelButton);
	CommandButtons->addStretch(1);
	connect(OKButton, SIGNAL(clicked()),this, SLOT(OnOK()));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(reject()));
	LeftSide->addStretch(1);
	LeftSide->addLayout(CommandButtons);
	LeftSide->addStretch(1);

	QHBoxLayout *MainLayout = new QHBoxLayout;
	MainLayout->addLayout(LeftSide);
	MainLayout->addStretch(1);
	MainLayout->addLayout(RightSide);
	setLayout(MainLayout);
	setMinimumWidth(400);
	setMinimumHeight(400);
}


void InterpolateFoilsDlg::InitDialog()
{

	int i;
	CFoil* pFoil;
	m_pctrlFoil1->clear();
	m_pctrlFoil2->clear();
	for (i=0; i<g_foilStore.size(); i++)
	{
		pFoil = g_foilStore.at(i);
		if(pFoil)
		{
            m_pctrlFoil1->addItem(pFoil->getName());
            m_pctrlFoil2->addItem(pFoil->getName());
		}
	}
	m_pctrlFoil1->setCurrentIndex(0);
	m_pctrlFoil2->setCurrentIndex(1);

	m_Frac = 0.0;
	m_pctrlFrac->setValue(100.0-m_Frac);

	OnSelChangeFoil1(0);
	OnSelChangeFoil2(1);
}


void InterpolateFoilsDlg::keyPressEvent(QKeyEvent *event)
{
	// Prevent Return Key from closing App
	// Generate the foil instead
	switch (event->key())
	{
		case Qt::Key_Return:
		{
			if(!OKButton->hasFocus() && !CancelButton->hasFocus())
			{
				Update();
				OKButton->setFocus();
			}
			else if (OKButton->hasFocus())
			{
				OnOK();
			}
			break;
		}
		case Qt::Key_Escape:
		{
			reject();
			return;
		}
		default:
			event->ignore();
	}
}


void InterpolateFoilsDlg::OnSelChangeFoil1(int /*i*/)
{
	MainFrame * pMainFrame = (MainFrame*)m_pMainFrame;
	QString strong  = m_pctrlFoil1->currentText();

//	i=0;
	CFoil* pFoil = pMainFrame->GetFoil(strong);

	if(pFoil)
	{
		QString str;
		str = QString(tr("Camb.=%1")).arg(pFoil->m_fCamber*100,5,'f',2);
		str += "%";
		strong = QString(tr(" at x=%1")).arg(pFoil->m_fXCamber*100,5,'f',1);
		strong += "%";
		str+=strong;
		m_pctrlCamb1->setText(str);

		str = QString(tr("Thick.=%1")).arg(pFoil->m_fThickness*100,5,'f',2);
		str += "%";
		strong = QString(tr(" at x=%1")).arg(pFoil->m_fXThickness*100,5,'f',1);
		strong += "%";
		str+=strong;
		m_pctrlThick1->setText(str);
	}
	Update();
}


void InterpolateFoilsDlg::OnSelChangeFoil2(int /*i*/)
{
//	i=0;
	MainFrame * pMainFrame = (MainFrame*)m_pMainFrame;
	QString strong  = m_pctrlFoil2->currentText();

	CFoil* pFoil = pMainFrame->GetFoil(strong);

	if(pFoil)
	{
		QString str;
		str = QString(tr("Camb.=%1")).arg(pFoil->m_fCamber*100,5,'f',2);
		str += "%";
		strong = QString(tr(" at x=%1")).arg(pFoil->m_fXCamber*100,5,'f',1);
		strong += "%";
		str+=strong;
		m_pctrlCamb2->setText(str);

		str = QString(tr("Thick.=%1")).arg(pFoil->m_fThickness*100,5,'f',2);
		str += "%";
		strong = QString(tr(" at x=%1")).arg(pFoil->m_fXThickness*100,5,'f',1);
		strong += "%";
		str+=strong;
		m_pctrlThick2->setText(str);
	}
	Update();
}


void InterpolateFoilsDlg::Update()
{
	MainFrame * pMainFrame = (MainFrame*)m_pMainFrame;
	QAFoil *pAFoil = (QAFoil*)m_pAFoil;
	QXDirect *pXDirect = (QXDirect*)m_pXDirect;
	XFoil *pXFoil = (XFoil*)s_pXFoil;
	QString strong;

	strong = m_pctrlFoil1->currentText();
	CFoil* pFoil1 = pMainFrame->GetFoil(strong);

	strong = m_pctrlFoil2->currentText();
	CFoil* pFoil2 = pMainFrame->GetFoil(strong);

	if(!pFoil1 || !pFoil2) return;

	pXFoil->Interpolate(pFoil1->x, pFoil1->y, pFoil1->n,
						pFoil2->x, pFoil2->y, pFoil2->n,
						m_Frac/100.0);

	for (int j=0; j< pFoil1->n; j++)
	{
		m_pBufferFoil->x[j]  = pXFoil->xb[j+1];
		m_pBufferFoil->y[j]  = pXFoil->yb[j+1];
		m_pBufferFoil->xb[j] = pXFoil->xb[j+1];
		m_pBufferFoil->yb[j] = pXFoil->yb[j+1];
	}
	m_pBufferFoil->n  = pFoil1->n;
	m_pBufferFoil->nb = pFoil1->n;

	m_pBufferFoil->InitFoil();

	QString str;
	str = QString(tr("Camb.=%1")).arg(m_pBufferFoil->m_fCamber*100,5,'f',2);
	str += "%";
	strong = QString(tr(" at x=%1")).arg(m_pBufferFoil->m_fXCamber*100,5,'f',1);
	strong += "%";
	str+=strong;
	m_pctrlCamb3->setText(str);

	str = QString(tr("Thick.=%1")).arg(m_pBufferFoil->m_fThickness*100,5,'f',2);
	str += "%";
	strong = QString(tr(" at x=%1")).arg(m_pBufferFoil->m_fXThickness*100,5,'f',1);
	strong += "%";
	str+=strong;
	m_pctrlThick3->setText(str);

	if(pXDirect) pXDirect->UpdateView();
	else if(pAFoil) pAFoil->UpdateView();
}


void InterpolateFoilsDlg::OnFrac()
{
	if(m_pctrlFrac->getValue()>100.0) m_pctrlFrac->setValue(100.0);
	if(m_pctrlFrac->getValue()<0.0)   m_pctrlFrac->setValue(0.0);

        m_Frac = m_pctrlFrac->getValue();
	m_pctrlSlider->setSliderPosition((int)m_Frac);
        m_Frac = 100.0 - m_Frac;
	Update();
}


void InterpolateFoilsDlg::OnOK()
{
	m_NewFoilName = m_pctrlNewFoilName->text();
    m_pBufferFoil->getName() = m_NewFoilName;

	QDialog::accept();
}


void InterpolateFoilsDlg::OnVScroll(int val)
{
	val = m_pctrlSlider->sliderPosition();
	m_Frac = 100.0 - (double)val;
	m_pctrlFrac->setValue(val);
	Update();
}










