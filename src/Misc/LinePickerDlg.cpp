/****************************************************************************

	LinePicker Class
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


#include <QKeyEvent>
#include <QColorDialog>
#include <QLabel>
#include <QGridLayout>

#include "LinePickerDlg.h"
#include "../Globals.h"



LinePickerDlg::LinePickerDlg(void *pParent)
{
	setWindowTitle(tr("Line Picker"));
	m_pParent = pParent;
	m_Style  = 0;
	m_Width  = 1;
	m_Color  = QColor(0,255,0);
	SetupLayout();

	m_pStyleDelegate = new LineDelegate(this);//will intercept painting operations
	m_pWidthDelegate = new LineDelegate(this);//will intercept painting operations

	m_pctrlStyle->setItemDelegate(m_pStyleDelegate);
	m_pctrlWidth->setItemDelegate(m_pWidthDelegate);

	connect(m_pctrlColor, SIGNAL(clicked()), this, SLOT(OnColor()));
	connect(m_pctrlStyle, SIGNAL(activated(int)), this, SLOT(OnStyle(int)));
	connect(m_pctrlWidth, SIGNAL(activated(int)), this, SLOT(OnWidth(int)));

	connect(OKButton, SIGNAL(clicked()),this, SLOT(accept()));
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}


void LinePickerDlg::FillBoxes()
{
	m_pctrlStyle->SetLine(m_Style, m_Width, m_Color);
	m_pctrlWidth->SetLine(m_Style, m_Width, m_Color);

	m_pctrlColor->SetStyle(m_Style);
	m_pctrlColor->SetWidth(m_Width);

	m_pctrlColor->SetColor(m_Color);
	m_pStyleDelegate->SetLineColor(m_Color);
	m_pWidthDelegate->SetLineColor(m_Color);

	m_pctrlStyle->setCurrentIndex(m_Style);
	m_pctrlWidth->setCurrentIndex(m_Width-1);

	int LineWidth[5];
	for (int i=0; i<5;i++) LineWidth[i] = m_Width;
	m_pStyleDelegate->SetLineWidth(LineWidth); // the same selected width for all styles


	int LineStyle[5];
	for (int i=0; i<5;i++) LineStyle[i] = m_Style;
	m_pWidthDelegate->SetLineStyle(LineStyle); //the same selected style for all widths
}


void LinePickerDlg::InitDialog(int style, int width, QColor color)
{
	m_Color = color;
	m_Width = width;
	m_Style = style;
	QString str;
	for (int i=0; i<5; i++)
	{
		str = QString("%1").arg(i);
		m_pctrlWidth->addItem(str);
	}
	m_pctrlStyle->addItem("solid");//string doesn't matter
	m_pctrlStyle->addItem("dash");
	m_pctrlStyle->addItem("dot");
	m_pctrlStyle->addItem("dashdot");
	m_pctrlStyle->addItem("dashdotdot");

	FillBoxes();
}

void LinePickerDlg::disableLineStyle() {
	m_pctrlStyle->setEnabled(false);
}



void LinePickerDlg::InitDialog()
{
	QString str;
	for (int i=0; i<5; i++)
	{
		str = QString("%1").arg(i);
		m_pctrlWidth->addItem(str);
	}
	m_pctrlStyle->addItem("solid");//string doesn't matter
	m_pctrlStyle->addItem("dash");
	m_pctrlStyle->addItem("dot");
	m_pctrlStyle->addItem("dashdot");
	m_pctrlStyle->addItem("dashdotdot");

	FillBoxes();
}



void LinePickerDlg::keyPressEvent(QKeyEvent *event)
{
	// Prevent Return Key from closing App
	switch (event->key())
	{
		case Qt::Key_Return:
		{
			if(!OKButton->hasFocus())
			{
				OKButton->setFocus();
				return;
			}
			else
			{
				accept();
				return;
			}
			break;
		}
		case Qt::Key_Escape:
		{
			reject();
			return;
		}
	}
}



void LinePickerDlg::OnWidth(int val)
{
	m_Width = val+1;
	FillBoxes();
	repaint();
	OKButton->setFocus();
}


void LinePickerDlg::OnStyle(int val)
{
	m_Style = val;
	FillBoxes();
	repaint();
	OKButton->setFocus();
}


void LinePickerDlg::OnColor()
{
	QColor Color = QColorDialog::getColor(m_Color,
								   this, "Color Selection", QColorDialog::ShowAlphaChannel);
	if(Color.isValid()) m_Color = Color;

	FillBoxes();
	repaint();
	OKButton->setFocus();
}



int & LinePickerDlg::GetWidth()
{
	return m_Width;
}


int & LinePickerDlg::GetStyle()
{
	return m_Style;
}


QColor & LinePickerDlg::GetColor()
{
	return m_Color;
}

void LinePickerDlg::SetColor(QColor color)
{
	m_Color = color;
	FillBoxes();
	repaint();
}





void LinePickerDlg::SetStyle(int style)
{
	m_Style = style;
	FillBoxes();
	repaint();
}


void LinePickerDlg::SetupLayout()
{
	QGridLayout *StyleLayout = new QGridLayout;
	QLabel *lab1 = new QLabel(tr("Style"));
	QLabel *lab2 = new QLabel(tr("Width"));
	QLabel *lab3 = new QLabel(tr("Color"));
	lab1->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
	lab2->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
	lab3->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
	lab1->setMinimumWidth(60);
	lab2->setMinimumWidth(60);
	lab3->setMinimumWidth(60);
	m_pctrlColor = new LineButton;
	m_pctrlStyle = new LineCbBox;
	m_pctrlWidth = new LineCbBox;
	m_pctrlStyle->setMinimumWidth(140);
	m_pctrlWidth->setMinimumWidth(140);
	m_pctrlColor->setMinimumWidth(140);

	StyleLayout->addWidget(lab1,1,1);
	StyleLayout->addWidget(lab2,2,1);
	StyleLayout->addWidget(lab3,3,1);
	StyleLayout->addWidget(m_pctrlStyle,1,2);
	StyleLayout->addWidget(m_pctrlWidth,2,2);
	StyleLayout->addWidget(m_pctrlColor,3,2);

	QHBoxLayout *CommandButtons = new QHBoxLayout;
	OKButton = new QPushButton(tr("OK"));
	CancelButton = new QPushButton(tr("Cancel"));
	CommandButtons->addStretch();
	CommandButtons->addWidget(CancelButton);
	CommandButtons->addWidget(OKButton);

	QVBoxLayout *MainLayout = new QVBoxLayout;
	MainLayout->addStretch(1);
	MainLayout->addLayout(StyleLayout);
	MainLayout->addStretch(1);
	MainLayout->addLayout(CommandButtons);
	MainLayout->addStretch(1);

	setMinimumHeight(170);
	setLayout(MainLayout);
}



void LinePickerDlg::SetWidth(int width)
{
	m_Width = width;
	FillBoxes();
	repaint();
}













