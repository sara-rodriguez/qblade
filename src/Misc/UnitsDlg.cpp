/****************************************************************************

	UnitsDlg Class
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

#include "../Globals.h"
#include "UnitsDlg.h"
#include <QGridLayout>
#include <QVBoxLayout>

UnitsDlg::UnitsDlg()
{
	m_Length  = 0;
	m_Area    = 1;
	m_Weight  = 1;
	m_Speed   = 0;
	m_Force   = 0;
        //////////new code DM////////
        m_Power   = 2;
        ///////////end new code DM///////////
	m_bLengthOnly = false;
	m_Question = tr("Select units for this project :");
	setWindowTitle(tr("Units Dialog"));
	SetupLayout();
}

void UnitsDlg::SetupLayout()
{
	QGridLayout *UnitsLayout = new QGridLayout;

	QLabel *lab1 = new QLabel(tr("Length"));
	QLabel *lab2 = new QLabel(tr("Area"));
	QLabel *lab3 = new QLabel(tr("Speed"));
	QLabel *lab4 = new QLabel(tr("Mass"));
	QLabel *lab5 = new QLabel(tr("Force"));
	QLabel *lab6 = new QLabel(tr("Moment"));
	UnitsLayout->addWidget(lab1, 1,1);
	UnitsLayout->addWidget(lab2, 2,1);
	UnitsLayout->addWidget(lab3, 3,1);
	UnitsLayout->addWidget(lab4, 4,1);
	UnitsLayout->addWidget(lab5, 5,1);
	UnitsLayout->addWidget(lab6, 6,1);

	m_pctrlQuestion = new QLabel(tr("Define the project units"));

	m_pctrlLengthFactor = new QLabel(" ");
	m_pctrlSurfaceFactor = new QLabel(" ");
	m_pctrlWeightFactor = new QLabel(" ");
	m_pctrlSpeedFactor = new QLabel(" ");
	m_pctrlForceFactor = new QLabel(" ");
	m_pctrlMomentFactor = new QLabel(" ");
	m_pctrlLengthFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlSurfaceFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlWeightFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlSpeedFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlForceFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlMomentFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);


	UnitsLayout->addWidget(m_pctrlLengthFactor, 1,2);
	UnitsLayout->addWidget(m_pctrlSurfaceFactor, 2,2);
	UnitsLayout->addWidget(m_pctrlSpeedFactor, 3,2);
	UnitsLayout->addWidget(m_pctrlWeightFactor, 4,2);
	UnitsLayout->addWidget(m_pctrlForceFactor, 5,2);
	UnitsLayout->addWidget(m_pctrlMomentFactor, 6,2);

	m_pctrlLength  = new QComboBox;
	m_pctrlSurface = new QComboBox;
	m_pctrlSpeed   = new QComboBox;
	m_pctrlWeight  = new QComboBox;
	m_pctrlForce   = new QComboBox;
	m_pctrlMoment  = new QComboBox;
	UnitsLayout->addWidget(m_pctrlLength,  1,3);
	UnitsLayout->addWidget(m_pctrlSurface, 2,3);
	UnitsLayout->addWidget(m_pctrlSpeed,   3,3);
	UnitsLayout->addWidget(m_pctrlWeight,  4,3);
	UnitsLayout->addWidget(m_pctrlForce,   5,3);
	UnitsLayout->addWidget(m_pctrlMoment,  6,3);


	m_pctrlLengthInvFactor = new QLabel(" ");
	m_pctrlSurfaceInvFactor = new QLabel(" ");
	m_pctrlWeightInvFactor = new QLabel(" ");
	m_pctrlSpeedInvFactor = new QLabel(" ");
	m_pctrlForceInvFactor = new QLabel(" ");
	m_pctrlMomentInvFactor = new QLabel(" ");
	m_pctrlLengthInvFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlSurfaceInvFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlWeightInvFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlSpeedInvFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlForceInvFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	m_pctrlMomentInvFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
	UnitsLayout->addWidget(m_pctrlLengthInvFactor, 1,4);
	UnitsLayout->addWidget(m_pctrlSurfaceInvFactor, 2,4);
	UnitsLayout->addWidget(m_pctrlSpeedInvFactor, 3,4);
	UnitsLayout->addWidget(m_pctrlWeightInvFactor, 4,4);
	UnitsLayout->addWidget(m_pctrlForceInvFactor, 5,4);
	UnitsLayout->addWidget(m_pctrlMomentInvFactor, 6,4);
	UnitsLayout->setColumnStretch(4,2);
	UnitsLayout->setColumnMinimumWidth(4,220);


	QHBoxLayout *CommandButtons = new QHBoxLayout;
	OKButton      = new QPushButton(tr("OK"));
	CancelButton  = new QPushButton(tr("Cancel"));
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(OKButton);
	CommandButtons->addStretch(1);
	CommandButtons->addWidget(CancelButton);
	CommandButtons->addStretch(1);

	QVBoxLayout *MainLayout = new QVBoxLayout;
	MainLayout->addWidget(m_pctrlQuestion);
	MainLayout->addLayout(UnitsLayout);
	MainLayout->addStretch(1);
	MainLayout->addSpacing(20);
	MainLayout->addLayout(CommandButtons);
	MainLayout->addStretch(1);

	setLayout(MainLayout);

	connect(OKButton, SIGNAL(clicked()),this, SLOT(accept()));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(reject()));

	connect(m_pctrlLength, SIGNAL(activated(const QString &)),this, SLOT(OnSelChanged(const QString &)));
	connect(m_pctrlSurface, SIGNAL(activated(const QString &)),this, SLOT(OnSelChanged(const QString &)));
	connect(m_pctrlSpeed, SIGNAL(activated(const QString &)),this, SLOT(OnSelChanged(const QString &)));
	connect(m_pctrlWeight, SIGNAL(activated(const QString &)),this, SLOT(OnSelChanged(const QString &)));
	connect(m_pctrlForce, SIGNAL(activated(const QString &)),this, SLOT(OnSelChanged(const QString &)));
	connect(m_pctrlMoment, SIGNAL(activated(const QString &)),this, SLOT(OnSelChanged(const QString &)));

        ////////////////////new code DM///////////////
        QLabel *lab7 = new QLabel(tr("Power"));
        UnitsLayout->addWidget(lab7, 7,1);
        m_pctrlPowerFactor = new QLabel(" ");
        m_pctrlPowerInvFactor = new QLabel(" ");

        m_pctrlPowerFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);
        m_pctrlPowerInvFactor->setAlignment(Qt::AlignRight | Qt::AlignCenter);

        UnitsLayout->addWidget(m_pctrlPowerFactor, 7,2);
        UnitsLayout->addWidget(m_pctrlMomentInvFactor, 7,4);

        m_pctrlPower = new QComboBox;
        UnitsLayout->addWidget(m_pctrlPower, 7,3);
        connect(m_pctrlPower, SIGNAL(activated(const QString &)),this, SLOT(OnSelChanged(const QString &)));
        //////////end new code DM//////////

}

void UnitsDlg::InitDialog()
{
	QStringList list;
	list <<"mm" << "cm"<<"dm"<<"m"<<"in"<<"ft";
	m_pctrlLength->clear();
	m_pctrlLength->addItems(list);		//5

	m_pctrlSurface->clear();
	m_pctrlSurface->addItem(QString::fromUtf8("mm²"));		//0
	m_pctrlSurface->addItem(QString::fromUtf8("cm²"));		//1
	m_pctrlSurface->addItem(QString::fromUtf8("dm²"));		//2
	m_pctrlSurface->addItem(QString::fromUtf8("m²"));		//3
	m_pctrlSurface->addItem(QString::fromUtf8("in²"));		//4
	m_pctrlSurface->addItem(QString::fromUtf8("ft²"));		//5

	m_pctrlSpeed->clear();
	m_pctrlSpeed->addItem("m/s");		//0
	m_pctrlSpeed->addItem("km/h");		//1
	m_pctrlSpeed->addItem("ft/s");		//2
	m_pctrlSpeed->addItem("kt (int.)");	//3
	m_pctrlSpeed->addItem("mph");		//4

	m_pctrlWeight->clear();
	m_pctrlWeight->addItem("g");		//0
	m_pctrlWeight->addItem("kg");		//1
	m_pctrlWeight->addItem("oz");		//2
	m_pctrlWeight->addItem("lb");		//3

	m_pctrlForce->clear();
	m_pctrlForce->addItem("N");		//0
	m_pctrlForce->addItem("lbf");		//1

	m_pctrlMoment->clear();
	m_pctrlMoment->addItem("N.m");	//0
	m_pctrlMoment->addItem("lbf.in");	//1
	m_pctrlMoment->addItem("lbf.ft");	//2

        //////////new code DM///////////
        m_pctrlPower->clear();
        m_pctrlPower->addItem("MW");
        m_pctrlPower->addItem("kW");
        m_pctrlPower->addItem("W");
        m_pctrlPower->setCurrentIndex(m_Power);
        if(m_bLengthOnly)
        {
        m_pctrlPower->setEnabled(false);
        }
        ///////////end new code DM//////////

	m_pctrlLength->setCurrentIndex(m_Length);
	m_pctrlWeight->setCurrentIndex(m_Weight);
	m_pctrlSurface->setCurrentIndex(m_Area);
	m_pctrlSpeed->setCurrentIndex(m_Speed);
	m_pctrlForce->setCurrentIndex(m_Force);
	m_pctrlMoment->setCurrentIndex(m_Moment);

	m_pctrlLength->setFocus();
	OnSelChanged(" ");

	if(m_bLengthOnly)
	{
		m_pctrlSpeed->setEnabled(false);
		m_pctrlSurface->setEnabled(false);
		m_pctrlWeight->setEnabled(false);
		m_pctrlForce->setEnabled(false);
		m_pctrlMoment->setEnabled(false);
	}
	m_pctrlQuestion->setText(m_Question);
}


void UnitsDlg::OnSelChanged(const QString &/*strong*/)
{
	m_Length  = m_pctrlLength->currentIndex();
	m_Area    = m_pctrlSurface->currentIndex();
	m_Weight  = m_pctrlWeight->currentIndex();
	m_Speed   = m_pctrlSpeed->currentIndex();
	m_Force   = m_pctrlForce->currentIndex();
	m_Moment  = m_pctrlMoment->currentIndex();
        ////new code DM//////////
        m_Power = m_pctrlPower->currentIndex();
        //////end new code DM//////////

        SetUnits(m_Length, m_Area, m_Speed, m_Weight, m_Force, m_Moment, m_Power,
                         m_mtoUnit, m_m2toUnit, m_mstoUnit, m_kgtoUnit, m_NtoUnit, m_NmtoUnit, m_WtoUnit);

	QString str, strange;

	GetLengthUnit(str,m_Length);
	strange= QString("     1 m = %1").arg(m_mtoUnit,15,'f',5);
	m_pctrlLengthFactor->setText(strange);
	strange= "1 "+str+" = " +QString("%1 m").arg(1./m_mtoUnit,15,'f',5);
	m_pctrlLengthInvFactor->setText(strange);


	GetAreaUnit(str,m_Area);
	strange= QString(QString::fromUtf8("     1 m² = %1")).arg(m_m2toUnit,15,'f',5);
	m_pctrlSurfaceFactor->setText(strange);
	strange= "1 "+str+" = " +QString("%1 m").arg(1./m_m2toUnit,15,'f',5);
	strange += QString::fromUtf8("²");
	m_pctrlSurfaceInvFactor->setText(strange);

	GetWeightUnit(str,m_Weight);
	strange= QString("     1 kg = %1").arg(m_kgtoUnit,15,'f',5);
	m_pctrlWeightFactor->setText(strange);
	strange= "1 "+str+" = " +QString("%1 kg").arg(1./m_kgtoUnit,15,'f',5);
	m_pctrlWeightInvFactor->setText(strange);

	GetSpeedUnit(str,m_Speed);
	strange= QString("     1 m/s = %1").arg(m_mstoUnit,15,'f',5);
	m_pctrlSpeedFactor->setText(strange);
	strange= "1 "+str+" = " +QString("%1 m/s").arg(1./m_mstoUnit,15,'f',5);
	m_pctrlSpeedInvFactor->setText(strange);

	GetForceUnit(str,m_Force);
	strange= QString("     1 N = %1").arg(m_NtoUnit,15,'f',5);
	m_pctrlForceFactor->setText(strange);
	strange= "1 "+str+" = " +QString("%1 N").arg(1./m_NtoUnit,15,'f',5);
	m_pctrlForceInvFactor->setText(strange);

	GetMomentUnit(str,m_Moment);
	strange= QString("     1 N.m = %1").arg(m_NmtoUnit,15,'f',5);
	m_pctrlMomentFactor->setText(strange);
	strange= "1 "+str+" = " +QString("%1 N.m").arg(1./m_NmtoUnit,15,'f',5);
	m_pctrlMomentInvFactor->setText(strange);

        /////////////new code DM//////////////
        GetPowerUnit(str,m_Power);
        strange= QString("     1 W = %1").arg(m_WtoUnit,15,'f',5);
        m_pctrlPowerFactor->setText(strange);
        strange= "1 "+str+" = " +QString("%1 W").arg(1./m_WtoUnit,15,'f',5);
		m_pctrlPowerInvFactor->setText(strange);
		///////////end new code DM////////////

}
