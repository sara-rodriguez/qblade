/****************************************************************************

	BatchThreadDlg Class
        Copyright (C) 2003-2010 Andre Deperrois adeperrois@xflr5.com

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

#include "BatchThreadDlg.h"
#include "XDirect.h"
#include "ReListDlg.h"
#include "../MainFrame.h"
#include "../XInverse/FoilSelectionDlg.h"
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QDir>
#include <QDateTime>
#include <QCoreApplication>
#include <QThread>
#include <QThreadPool>
#include <QDebug>
#include "../Store.h"


//bool BatchThreadDlg::s_bStoreOpp = false;
bool BatchThreadDlg::s_bCurrentFoil=true;
bool BatchThreadDlg::s_bUpdatePolarView = false;
void * BatchThreadDlg::s_pXFoil;
void * BatchThreadDlg::s_pMainFrame;
void * BatchThreadDlg::s_pXDirect;


BatchThreadDlg::BatchThreadDlg(void */*pParent*/)
{
	QString str = tr("Multi-threaded batch analysis");
	setWindowTitle(str);

	m_FoilList.clear();
	
	m_ReInc =  50000.0;
	m_ReMax = 300000.0;
	m_ReMin = 100000.0;
	m_Mach  = 0.0;

	m_bAlpha          = true;
	m_bFromList       = false;
	m_bFromZero       = false;
	m_bInitBL         = false;
	m_bIsRunning      = false;
	m_bCancel         = false;

	XFoil::s_bCancel = false;

	SetupLayout();
	connect(m_pctrlFoil1, SIGNAL(clicked()), this, SLOT(OnFoilSelectionType()));
	connect(m_pctrlFoil2, SIGNAL(clicked()), this, SLOT(OnFoilSelectionType()));
	connect(m_pctrlFoilList, SIGNAL(clicked()), this, SLOT(OnFoilList()));
	connect(m_pctrlClose, SIGNAL(clicked()), this, SLOT(OnClose()));
	connect(m_pctrlAnalyze, SIGNAL(clicked()), this, SLOT(OnAnalyze()));
	connect(m_pctrlAlpha, SIGNAL(toggled(bool)), this, SLOT(OnAcl()));
	connect(m_pctrlCl, SIGNAL(toggled(bool)), this, SLOT(OnAcl()));
	connect(m_rbRange1, SIGNAL(toggled(bool)), this, SLOT(OnRange()));
	connect(m_rbRange2, SIGNAL(toggled(bool)), this, SLOT(OnRange()));
	connect(m_pctrlEditList, SIGNAL(clicked()), this, SLOT(OnEditReList()));
	connect(m_pctrlFromZero, SIGNAL(stateChanged(int)), this, SLOT(OnFromZero(int)));
	connect(m_pctrlSpecMin, SIGNAL(editingFinished()), this, SLOT(OnSpecChanged()));
	connect(m_pctrlSpecMax, SIGNAL(editingFinished()), this, SLOT(OnSpecChanged()));
	connect(m_pctrlSpecDelta, SIGNAL(editingFinished()), this, SLOT(OnSpecChanged()));
}



void BatchThreadDlg::SetupLayout()
{
	QSizePolicy szPolicyExpanding;
	szPolicyExpanding.setHorizontalPolicy(QSizePolicy::Expanding);
	szPolicyExpanding.setVerticalPolicy(QSizePolicy::Expanding);

	QSizePolicy szPolicyMinimum;
	szPolicyMinimum.setHorizontalPolicy(QSizePolicy::Minimum);
	szPolicyMinimum.setVerticalPolicy(QSizePolicy::MinimumExpanding);

	QSizePolicy szPolicyMaximum;
	szPolicyMaximum.setHorizontalPolicy(QSizePolicy::Maximum);
	szPolicyMaximum.setVerticalPolicy(QSizePolicy::Maximum);

	QGroupBox *FoilBox = new QGroupBox(tr("Foil Selection"));
	{
		QHBoxLayout *FoilLayout = new QHBoxLayout;
		m_pctrlFoil1 = new QRadioButton(tr("Current foil only"));
		m_pctrlFoil2 = new QRadioButton(tr("Foil list"));
		m_pctrlFoilList = new QPushButton(tr("Foil list"));
		FoilLayout->addWidget(m_pctrlFoil1);
		FoilLayout->addWidget(m_pctrlFoil2);
		FoilLayout->addStretch(1);
		FoilLayout->addWidget(m_pctrlFoilList);
		FoilBox->setLayout(FoilLayout);
	}

	QGroupBox *BatchVarsGroup = new QGroupBox(tr("Batch Variables"));
	{
		QGridLayout *BatchVars = new QGridLayout;
		m_rbRange1 = new QRadioButton(tr("Range"));
		m_rbRange2 = new QRadioButton(tr("Re List"));
		m_pctrlEditList = new QPushButton(tr("Edit List"));
		QHBoxLayout *RangeSpecLayout = new QHBoxLayout;
		RangeSpecLayout->addWidget(m_rbRange1);
		RangeSpecLayout->addWidget(m_rbRange2);
		RangeSpecLayout->addStretch(1);
		RangeSpecLayout->addWidget(m_pctrlEditList);

		QLabel *MinVal   = new QLabel(tr("Min"));
		QLabel *MaxVal   = new QLabel(tr("Max"));
		QLabel *DeltaVal = new QLabel(tr("Increment"));
		MinVal->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
		MaxVal->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
		DeltaVal->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);

		m_pctrlReType  = new QLabel("Reynolds=");
		m_pctrlMaType  = new QLabel("Mach=");
		QLabel *NCritLabel = new QLabel(tr("NCrit="));
		m_pctrlReType->setAlignment(Qt::AlignVCenter|Qt::AlignRight);
		NCritLabel->setAlignment(Qt::AlignVCenter|Qt::AlignRight);
		m_pctrlMaType->setAlignment(Qt::AlignVCenter|Qt::AlignRight);
		m_pctrlNCrit   = new NumberEdit();
		m_pctrlNCrit->setValue(9);
		m_pctrlReMin   = new NumberEdit(NumberEdit::Standard, 0);
		m_pctrlReMin->setValue(100000);
		m_pctrlReMax   = new NumberEdit(NumberEdit::Standard, 0);
		m_pctrlReMax->setValue(150000);
		m_pctrlReDelta = new NumberEdit(NumberEdit::Standard, 0);
		m_pctrlReDelta->setValue(50000);
		m_pctrlMach    = new NumberEdit();

		BatchVars->addWidget(MinVal, 2, 2);
		BatchVars->addWidget(MaxVal, 2, 3);
		BatchVars->addWidget(DeltaVal, 2, 4);
		BatchVars->addWidget(m_pctrlReType, 3, 1);
		BatchVars->addWidget(m_pctrlReMin, 3, 2);
		BatchVars->addWidget(m_pctrlReMax, 3, 3);
		BatchVars->addWidget(m_pctrlReDelta, 3, 4);
		BatchVars->addWidget(m_pctrlMaType, 4, 1);
		BatchVars->addWidget(m_pctrlMach, 4, 2);
		BatchVars->addWidget(NCritLabel, 5,1);
		BatchVars->addWidget(m_pctrlNCrit, 5, 2);

		QVBoxLayout *BatchVarsGroupLayout = new QVBoxLayout;
		BatchVarsGroupLayout->addLayout(RangeSpecLayout);
		BatchVarsGroupLayout->addLayout(BatchVars);
		BatchVarsGroup->setLayout(BatchVarsGroupLayout);
	}

	QGroupBox *RangeVarsGroup = new QGroupBox(tr("Analysis Range"));
	{
		QLabel *Spec = new QLabel(tr("Specify:"));
		m_pctrlAlpha = new QRadioButton(tr("Alpha"));
		m_pctrlCl = new QRadioButton(tr("Cl"));
		m_pctrlFromZero   = new QCheckBox(tr("From Zero"));
		QHBoxLayout *RangeSpecLayout = new QHBoxLayout;
		RangeSpecLayout->addWidget(Spec);
		RangeSpecLayout->addWidget(m_pctrlAlpha);
//		RangeSpecLayout->addWidget(m_pctrlCl);
		RangeSpecLayout->addStretch(1);
		RangeSpecLayout->addWidget(m_pctrlFromZero);

		QLabel *SpecMin   = new QLabel(tr("Min"));
		QLabel *SpecMax   = new QLabel(tr("Max"));
		QLabel *SpecDelta = new QLabel(tr("Increment"));
		SpecMin->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
		SpecMax->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
		SpecDelta->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
		m_pctrlSpecVar    = new QLabel(tr("Spec ="));
		m_pctrlSpecMin    = new NumberEdit();
		m_pctrlSpecMax    = new NumberEdit();
		m_pctrlSpecMax->setValue(1);
		m_pctrlSpecDelta  = new NumberEdit();
		m_pctrlSpecDelta->setValue(0.5);

		QGridLayout *RangeVarsLayout = new QGridLayout;
		RangeVarsLayout->addWidget(SpecMin, 1, 2);
		RangeVarsLayout->addWidget(SpecMax, 1, 3);
		RangeVarsLayout->addWidget(SpecDelta, 1, 4);
		RangeVarsLayout->addWidget(m_pctrlSpecVar, 2, 1);
		RangeVarsLayout->addWidget(m_pctrlSpecMin, 2, 2);
		RangeVarsLayout->addWidget(m_pctrlSpecMax, 2, 3);
		RangeVarsLayout->addWidget(m_pctrlSpecDelta, 2, 4);

		QVBoxLayout *RangeVarsGroupLayout = new QVBoxLayout;
		RangeVarsGroupLayout->addLayout(RangeSpecLayout);
		RangeVarsGroupLayout->addLayout(RangeVarsLayout);
		RangeVarsGroup->setLayout(RangeVarsGroupLayout);
	}

	QGroupBox *TransVarsGroup = new QGroupBox(tr("Forced Transitions"));
	{
		QGridLayout *TransVars = new QGridLayout;
		TransVars->setColumnStretch(0,4);
		TransVars->setColumnStretch(1,1);
		QLabel *TopTransLabel = new QLabel(tr("Top transition location (x/c)"));
		QLabel *BotTransLabel = new QLabel(tr("Bottom transition location (x/c)"));
		TopTransLabel->setAlignment(Qt::AlignVCenter|Qt::AlignRight);
		BotTransLabel->setAlignment(Qt::AlignVCenter|Qt::AlignRight);
		m_pctrlXTopTr = new NumberEdit();
		m_pctrlXTopTr->setValue(1);
		m_pctrlXBotTr = new NumberEdit();
		m_pctrlXBotTr->setValue(1);

		TransVars->addWidget(TopTransLabel, 2, 1);
		TransVars->addWidget(m_pctrlXTopTr, 2, 2);
		TransVars->addWidget(BotTransLabel, 3, 1);
		TransVars->addWidget(m_pctrlXBotTr, 3, 2);
		TransVarsGroup->setLayout(TransVars);
	}


	QHBoxLayout *OptionsLayout = new QHBoxLayout;
	m_pctrlInitBL          = new QCheckBox(tr("Initialize BLs between polars"));
	m_pctrlUpdatePolarView = new QCheckBox(tr("Update polar view"));
	m_pctrlUpdatePolarView->setToolTip(tr("Update the polar graphs after the completion of each foil/polar pair"));
	OptionsLayout->addWidget(m_pctrlInitBL);
	OptionsLayout->addStretch(1);
	OptionsLayout->addWidget(m_pctrlUpdatePolarView);

	//_*_*_*_*_*_*_**_*_*_**_*_*_*_
	m_pctrlTextOutput = new QTextEdit;
	m_pctrlTextOutput->setReadOnly(true);
	m_pctrlTextOutput->setLineWrapMode(QTextEdit::NoWrap);
	m_pctrlTextOutput->setWordWrapMode(QTextOption::NoWrap);
	m_pctrlTextOutput->setSizePolicy(szPolicyMinimum);
	m_pctrlTextOutput->setMinimumWidth(350);

	QHBoxLayout *CommandButtons = new QHBoxLayout;
	m_pctrlClose     = new QPushButton(tr("Close"));
	m_pctrlAnalyze   = new QPushButton(tr("Analyze"))	;
	m_pctrlAnalyze->setAutoDefault(true);

	CommandButtons->addStretch(1);
	CommandButtons->addWidget(m_pctrlAnalyze);
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(m_pctrlClose);
	CommandButtons->addStretch(1);

	QVBoxLayout *LeftSide = new QVBoxLayout;
	LeftSide->addWidget(FoilBox);
	LeftSide->addWidget(BatchVarsGroup);
	LeftSide->addWidget(TransVarsGroup);
	LeftSide->addWidget(RangeVarsGroup);
	LeftSide->addStretch(1);
	LeftSide->addSpacing(20);
	LeftSide->addLayout(CommandButtons);

	QVBoxLayout *RightSide = new QVBoxLayout;
//	RightSide->addStretch(1);
	RightSide->addLayout(OptionsLayout);
	RightSide->addWidget(m_pctrlTextOutput,1);
//	RightSide->addStretch(1);


	QHBoxLayout *BoxesLayout = new QHBoxLayout;
	BoxesLayout->setSpacing(10);
	BoxesLayout->addLayout(LeftSide);
	BoxesLayout->addLayout(RightSide);


	setLayout(BoxesLayout);
}


void BatchThreadDlg::AddOpPoint(CPolar *pPolar)
{
	XFoil *pXFoil = (XFoil*)s_pXFoil;

	pPolar->AddData(pXFoil);
}



void BatchThreadDlg::CleanUp()
{
	if(m_pXFile->isOpen())
	{
		QTextStream out(m_pXFile);
		out<<m_pctrlTextOutput->toPlainText();
		m_pXFile->close();
	}
	m_pctrlClose->setEnabled(true);
	m_pctrlAnalyze->setText(tr("Analyze"));
	m_bIsRunning = false;
	m_bCancel    = false;
	XFoil::s_bCancel = false;
	m_pctrlClose->setFocus();
}


CPolar * BatchThreadDlg::CreatePolar(CFoil *pFoil, double Spec, double Mach, double NCrit)
{
	if(!pFoil)
		return NULL;

	MainFrame *pMainFrame = (MainFrame*)s_pMainFrame;
	CPolar *pNewPolar = new CPolar;
	pNewPolar->setSingleParent(pFoil);
	pNewPolar->m_bIsVisible = true;

	pNewPolar->m_PolarType = m_PolarType;

	switch (pNewPolar->m_PolarType)
	{
		case 1:
			pNewPolar->m_MaType = 1;
			pNewPolar->m_ReType = 1;
			break;
		case 2:
			pNewPolar->m_MaType = 2;
			pNewPolar->m_ReType = 2;
			break;
		case 3:
			pNewPolar->m_MaType = 1;
			pNewPolar->m_ReType = 3;
			break;
		case 4:
			pNewPolar->m_MaType = 1;
			pNewPolar->m_ReType = 1;
			break;
		default:
			pNewPolar->m_ReType = 1;
			pNewPolar->m_MaType = 1;
			break;
	}
	if(m_PolarType!=FIXEDAOAPOLAR)  pNewPolar->m_Reynolds = Spec;
	else            pNewPolar->m_ASpec    = Spec;

	pNewPolar->m_Mach  = Mach;
	pNewPolar->m_ACrit = NCrit;
	pNewPolar->m_XTop  = m_XTopTr;
	pNewPolar->m_XBot  = m_XBotTr;
	pNewPolar->m_Color = pMainFrame->GetColor(1);

	SetPlrName(pNewPolar);
	
	CPolar *pOldPolar = g_polarStore.getObjectByName(pNewPolar->getName(), pFoil);
	if(pOldPolar) {
        delete pNewPolar;
        pNewPolar = pOldPolar;
	} else {
		pNewPolar = (g_polarStore.add(pNewPolar) ? pNewPolar : NULL);
    }
	return pNewPolar;
}


void BatchThreadDlg::keyPressEvent(QKeyEvent *event)
{
	// Prevent Return Key from closing App
	switch (event->key())
	{
		case Qt::Key_Return:
		{
			if(m_pctrlClose->hasFocus())	     done(1);
			else if(m_pctrlAnalyze->hasFocus())  OnAnalyze();
			else                                 m_pctrlAnalyze->setFocus();
			break;
		}
		case Qt::Key_Escape:
		{
			if(m_bIsRunning)
			{
				OnAnalyze();//will cancel the threads
			}
			else
			{
				OnClose(); // will close the dialog box
			}
			break;
		}
		default:
			event->ignore();
	}
}


void BatchThreadDlg::InitDialog()
{
	if(!m_pCurFoil) return;

	m_PolarType = FIXEDSPEEDPOLAR; //no choice...

	m_pctrlFoil1->setChecked(s_bCurrentFoil);
	m_pctrlFoil2->setChecked(!s_bCurrentFoil);
	m_pctrlFoilList->setEnabled(!s_bCurrentFoil);
	
	m_pctrlMach->setAutomaticPrecision(2);


	m_pctrlReMin->setAutomaticPrecision(0);
	m_pctrlReMax->setAutomaticPrecision(0);
	m_pctrlReDelta->setAutomaticPrecision(0);

	m_pctrlSpecMin->setAutomaticPrecision(1);
	m_pctrlSpecMax->setAutomaticPrecision(1);
	m_pctrlSpecDelta->setAutomaticPrecision(1);

	if(m_ReMin<=0.0) m_ReMin = fabs(m_ReInc);
	m_pctrlReMin->setValue(m_ReMin);
	m_pctrlReMax->setValue(m_ReMax);
	m_pctrlReDelta->setValue(m_ReInc);
	m_pctrlSpecMin->setValue(m_AlphaMin);
	m_pctrlSpecMax->setValue(m_AlphaMax);
	m_pctrlSpecDelta->setValue(m_AlphaInc);


	m_pctrlMach->setValue(m_Mach);
	m_pctrlNCrit->setValue(m_NCrit);
	m_pctrlXTopTr->setValue(m_XTopTr);
	m_pctrlXBotTr->setValue(m_XBotTr);

	if(m_bAlpha) m_pctrlAlpha->setChecked(true);
	else         m_pctrlCl->setChecked(m_bAlpha);
	OnAcl();

	if(!m_bFromList)  m_rbRange1->setChecked(true);
	else              m_rbRange2->setChecked(true);
	OnRange();

	if(m_bFromZero)  m_pctrlFromZero->setChecked(true);
	else             m_pctrlFromZero->setChecked(false);

	m_pctrlInitBL->setChecked(true);
	m_pctrlUpdatePolarView->setChecked(s_bUpdatePolarView);
}



void BatchThreadDlg::OnAcl()
{
	if(m_PolarType==FIXEDAOAPOLAR) return;
	m_bAlpha = m_pctrlAlpha->isChecked();
	if(m_bAlpha)
	{
		m_pctrlSpecVar->setText(tr("Alpha"));
		m_pctrlSpecMin->setValue(m_AlphaMin);
		m_pctrlSpecMax->setValue(m_AlphaMax);
		m_pctrlSpecDelta->setValue(m_AlphaInc);
		m_pctrlFromZero->setEnabled(true);
	}
	else
	{
		m_pctrlSpecVar->setText(tr("CL"));
		m_pctrlSpecMin->setValue(m_ClMin);
		m_pctrlSpecMax->setValue(m_ClMax);
		m_pctrlSpecDelta->setValue(m_ClInc);
		m_bFromZero = false;
		m_pctrlFromZero->setChecked(false);
		m_pctrlFromZero->setEnabled(false);
	}
}


void BatchThreadDlg::OnSpecChanged()
{
	ReadParams();
}



void BatchThreadDlg::OnAnalyze()
{
	if(m_bIsRunning)
	{
		m_bCancel = true;
		XFoilTask::s_bCancel = true;
		XFoil::s_bCancel = true;
		return;
	}
	m_bCancel    = false;
	m_bIsRunning = true;

	XFoil *pXFoil = (XFoil*)s_pXFoil;

	m_pctrlTextOutput->setText("");

	m_pctrlClose->setEnabled(false);

	QString FileName = QDir::tempPath() + "/XFLR5.log";
	m_pXFile = new QFile(FileName);
	if (!m_pXFile->open(QIODevice::WriteOnly | QIODevice::Text)) m_pXFile = NULL;

	pXFoil->pXFile = m_pXFile;
	pXFoil->lvisc = true;

	ReadParams();

	SetFileHeader();
	m_bInitBL = m_pctrlInitBL->isChecked();

	m_pctrlAnalyze->setFocus();
	StartAnalysis();
}



void BatchThreadDlg::OnClose()
{
	s_bUpdatePolarView = m_pctrlUpdatePolarView->isChecked();
	m_bCancel = true;
	XFoilTask::s_bCancel = true;
	QThreadPool::globalInstance()->waitForDone();
	ReadParams();
	accept();
}



void BatchThreadDlg::OnEditReList()
{
	ReListDlg dlg;

	for (int i=0; i<m_NRe; i++)
	{
		dlg.m_ReList[i]    = m_ReList[i];
		dlg.m_MachList[i]  = m_MachList[i];
		dlg.m_NCritList[i] = m_NCritList[i];
	}
	dlg.m_NRe = m_NRe;
	dlg.InitDialog();

	if(QDialog::Accepted == dlg.exec())
	{
		for (int i=0; i<dlg.m_NRe; i++)
		{
			m_ReList[i]    = dlg.m_ReList[i];
			m_MachList[i]  = dlg.m_MachList[i];
			m_NCritList[i] = dlg.m_NCritList[i];
		}
		m_NRe = dlg.m_NRe;
	}
}


void BatchThreadDlg::OnFoilList()
{
	FoilSelectionDlg dlg;
	dlg.SetSelectionMode(true);

	dlg.m_FoilList.clear();
	for(int i=0; i<m_FoilList.size(); i++)
	{
		dlg.m_FoilList.append(m_FoilList.at(i));
	}
	dlg.InitDialog();

	m_FoilList.clear();
	if(QDialog::Accepted == dlg.exec())
	{
		for(int i=0; i<dlg.m_FoilList.count();i++)
		{
			m_FoilList.append(dlg.m_FoilList.at(i));
		}
	}
}



void BatchThreadDlg::OnFoilSelectionType()
{
	s_bCurrentFoil = m_pctrlFoil1->isChecked();
	m_pctrlFoilList->setEnabled(!s_bCurrentFoil);
}


void BatchThreadDlg::OnFromZero(int /*state*/)
{
//	state = 0;
	if(m_pctrlFromZero->isChecked()) m_bFromZero = true;
	else                             m_bFromZero = false;
}

void BatchThreadDlg::OnInitBL(int /*state*/)
{
//	state = 0;
	if (m_pctrlInitBL->isChecked()) m_bInitBL = true;
	else                            m_bInitBL = false;
}


void BatchThreadDlg::OnRange()
{
	if(m_rbRange1->isChecked())
		m_bFromList = false;
	else
		m_bFromList = true;

	m_pctrlEditList->setEnabled(m_bFromList);
	m_pctrlReMin->setEnabled(!m_bFromList);
	m_pctrlReMax->setEnabled(!m_bFromList);
	m_pctrlReDelta->setEnabled(!m_bFromList);
	m_pctrlMach->setEnabled(!m_bFromList);
	m_pctrlNCrit->setEnabled(!m_bFromList);
}



void BatchThreadDlg::ReadParams()
{
	m_bAlpha = m_pctrlAlpha->isChecked();

	if(m_PolarType!=FIXEDAOAPOLAR)
	{
		m_ReInc = m_pctrlReDelta->getValue();
		m_ReMax = m_pctrlReMax->getValue();
		m_ReMin = m_pctrlReMin->getValue();

		if(m_bAlpha)
		{
			m_AlphaInc = fabs(m_pctrlSpecDelta->getValue());
			m_AlphaMax = m_pctrlSpecMax->getValue();
			m_AlphaMin = m_pctrlSpecMin->getValue();
		}
		else
		{
			m_ClInc = fabs(m_pctrlSpecDelta->getValue());
			m_ClMin = m_pctrlSpecMin->getValue();
			m_ClMax = m_pctrlSpecMax->getValue();
		}
	}

	if(m_ReMin<=0.0) m_ReMin = fabs(m_ReInc);
	if(m_ReMax<=0.0) m_ReMax = fabs(m_ReMax);


	m_Mach     = m_pctrlMach->getValue();
	if(m_Mach<=0.0) m_Mach = 0.0;
	m_NCrit    = m_pctrlNCrit->getValue();
	m_XTopTr   = m_pctrlXTopTr->getValue();
	m_XBotTr   = m_pctrlXBotTr->getValue();
	
}




void BatchThreadDlg::SetFileHeader()
{
	QXDirect *pXDirect = (QXDirect*)s_pXDirect;
	MainFrame *pMainFrame = (MainFrame*)s_pMainFrame;

	QTextStream out(m_pXFile);

	out << "\n";
	out << pMainFrame->m_VersionName;
	out << "\n";
    out << m_pCurFoil->getName();
	out << "\n";
	if(pXDirect && pXDirect->m_pCurPolar)
	{
//	out << pXDirect->m_pCurPolar->m_PlrName;
//	out << "\n";
	}

	QDateTime dt = QDateTime::currentDateTime();
	QString str = dt.toString("dd.MM.yyyy  hh:mm:ss");

	out << str;
	out << "\n___________________________________\n\n";

}



void BatchThreadDlg::SetPlrName(CPolar *pNewPolar)
{
	if(m_PolarType!=FIXEDAOAPOLAR)
	{
		double R = pNewPolar->m_Reynolds/1000000.;
		pNewPolar->setName(QString("T%1_Re%2_M%3")
								 .arg(pNewPolar->m_PolarType)
								 .arg(R,0,'f',3)
								 .arg( pNewPolar->m_Mach,0,'f',2));
	}
	else
	{
		pNewPolar->setName(QString("T%1_Al%2_M%3")
								 .arg(pNewPolar->m_PolarType)
								 .arg(pNewPolar->m_ASpec,5,'f',2)
								 .arg(pNewPolar->m_Mach,0,'f',2));
	}
	QString str;
	str = QString("_N%1").arg(pNewPolar->m_ACrit,0,'f',1);
	pNewPolar->setName(pNewPolar->getName() + str);
}


void BatchThreadDlg::UpdateOutput(QString &str)
{
	m_pctrlTextOutput->insertPlainText(str);
	m_pctrlTextOutput->ensureCursorVisible();
}




void BatchThreadDlg::WriteString(QString &strong)
{
	if(!m_pXFile) return;
	if(!m_pXFile->isOpen()) return;
	QTextStream ds(m_pXFile);
	ds << strong;
}



void BatchThreadDlg::StartAnalysis()
{
	QXDirect *pXDirect = (QXDirect*)s_pXDirect;
	MainFrame *pMainFrame = (MainFrame*)s_pMainFrame;
	CFoil *pFoil;
	CPolar *pPolar;
	QString strong;
	int iRe, nRe;

	m_pctrlAnalyze->setText(tr("Cancel"));
	strong ="Launching multi-threaded batch analysis\n\n";
	UpdateOutput(strong);

	if(s_bCurrentFoil)
	{
		m_FoilList.clear();
        m_FoilList.append(m_pCurFoil->getName());
	}

	if(!m_bFromList) nRe = (int)fabs((m_ReMax-m_ReMin)/m_ReInc);
	else             nRe = m_NRe-1;

//	QThreadPool::globalInstance()->setExpiryTimeout(60000);//ms

	//build an array of all analysis pairs to run
	m_nAnalysis = 0;
	Analysis *pAnalysis=NULL;
	for(int i=0; i<m_FoilList.count(); i++)
	{
		pFoil = pMainFrame->GetFoil(m_FoilList.at(i));
		if(pFoil)
		{
			for (iRe=0; iRe<=nRe; iRe++)
			{
				pAnalysis = new Analysis;
				m_AnalysisPair.append(pAnalysis);
				pAnalysis->pFoil = pFoil;

				if(!m_bFromList) pPolar = CreatePolar(pFoil, m_ReMin + iRe *m_ReInc, m_Mach, m_NCrit);
				else             pPolar = CreatePolar(pFoil, m_ReList[iRe], m_MachList[iRe], m_NCritList[iRe]);
				pAnalysis->pPolar=pPolar;

				m_nAnalysis++;
			}
		}
	}
	strong = QString("Found %1 foil/polar pairs to analyze\n").arg(m_nAnalysis);
	UpdateOutput(strong);

	//Start as many threads as the system will support
	m_nThreads = QThread::idealThreadCount();
	m_pXFoilTask = new XFoilTask[m_nThreads];
	XFoilTask::s_bAutoInitBL= pXDirect->m_bAutoInitBL;
	XFoilTask::s_pBatchThreadDlg = this;
	XFoilTask::s_bCancel = false;
	strong = QString("Using %1 threads\n\n").arg(m_nThreads);
	UpdateOutput(strong);

	for(int it=0; it<m_nThreads; it++)	StartThread();

}


void BatchThreadDlg::timerEvent(QTimerEvent *event)
{
	// used out of purpose as a convenient way to receive a notification
	// of the analysis pair which has ended its analysis
	qApp->processEvents();

	if(m_pXFoilTask)
	{
        QString str = "   ...Finished "+m_pXFoilTask[event->timerId()].m_pFoil->getName()+" / "
								 +m_pXFoilTask[event->timerId()].m_pPolar->getName()+"\n";

		UpdateOutput(str);
		m_pXFoilTask[event->timerId()].m_bIsFinished = true;
	}

	//time to launch another analysis, if any is left
	QString strong;
	QXDirect *pXDirect = (QXDirect*)s_pXDirect;
	if(m_pctrlUpdatePolarView->isChecked())
	{
		pXDirect->CreatePolarCurves();
		pXDirect->UpdateView();
	}

	if(m_AnalysisPair.size()<=0 || m_bCancel)
	{
		//nothing left to launch... just wait and enjoy the show
		bool bAllFinished = true;
		for (int it=0; it<m_nThreads; it++)
		{
			if(m_pXFoilTask && !m_pXFoilTask[it].m_bIsFinished)
			{
				bAllFinished = false;
				break;
			}
		}

		if(bAllFinished)
		{
			QThreadPool::globalInstance()->waitForDone();

			if(m_pXFoilTask)
			{
				delete [] m_pXFoilTask;
				m_pXFoilTask = NULL;

				if(m_bCancel) strong = "\n_____Analysis cancelled_____\n";
				else          strong = "\n_____Analysis completed_____\n";
				UpdateOutput(strong);

				CleanUp();

				if(pXDirect->m_bPolar)
				{
					pXDirect->CreatePolarCurves();
					pXDirect->UpdateView();
				}
			}
		}
	}
	else  StartThread(); // analyze a new pair
}



void BatchThreadDlg::StartThread()
{
	Analysis *pAnalysis;
	QString strong;
	//  browse through the array until we find an available thread
	if(QThreadPool::globalInstance()->activeThreadCount()<m_nThreads)
	{
		//re-use a runnable object, don't recreate one, to save memory allocation time
		//on the other hand, this loop will run permanently;
		for (int it=0; it<m_nThreads; it++)
		{
            if(m_pXFoilTask && m_pXFoilTask[it].m_bIsFinished && m_AnalysisPair.size())
			{
				m_pXFoilTask[it].m_bIsFinished = false;
				m_pXFoilTask[it].m_Id = it;

				//take the last analysis in the array
				pAnalysis = (Analysis*)m_AnalysisPair.at(m_AnalysisPair.size()-1);

				//initiate the task
				m_pXFoilTask[it].Init(pAnalysis->pFoil, pAnalysis->pPolar);

				//launch it
				strong = "Starting "+pAnalysis->pFoil->getName()+" / "+pAnalysis->pPolar->getName()+"\n";
				UpdateOutput(strong);
				QThreadPool::globalInstance()->start(m_pXFoilTask+it);

				//remove it from the todo array
				m_AnalysisPair.removeLast();
				break;
			}
		}
	}
}








