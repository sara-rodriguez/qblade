/****************************************************************************

    Spline Foil Class
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



#include "math.h"
#include "Sf.h"
#include "../Globals.h"



CSF::CSF()
{
	m_FoilStyle = 0;
	m_FoilWidth = 1;
	m_FoilColor   = QColor(119, 183, 83);
	m_OutPoints    = 0;
	m_bVisible     = true;
	m_bOutPoints   = false;
	m_bModified    = false;
	m_bCenterLine  = false;
	m_Intrados.SetSplineParams(m_FoilStyle, m_FoilWidth, m_FoilColor);
	m_Extrados.SetSplineParams(m_FoilStyle, m_FoilWidth, m_FoilColor);
	m_bSymetric = false;
}


CSF::~CSF()
{
}


void CSF::SetCurveParams(int style, int width, QColor color)
{
	m_FoilStyle = style;
	m_FoilWidth = width;
	m_FoilColor = color;
	m_Intrados.SetSplineParams(style, width, color);
	m_Extrados.SetSplineParams(style, width, color);
}


bool CSF::InitSplineFoil()
{
	m_bModified   = false;
	m_strFoilName = QObject::tr("Spline Foil");
	
	m_Extrados.m_Color = m_FoilColor;
	m_Extrados.m_Style = m_FoilStyle;
	m_Extrados.m_Width = m_FoilWidth;
	m_Extrados.m_iCtrlPoints = 0;
	m_Extrados.m_iRes = 30;
	m_Extrados.InsertPoint(0.0 , 0.0);
	m_Extrados.InsertPoint(0.0 , 0.00774066);
	m_Extrados.InsertPoint(0.0306026, 0.0343829);
	m_Extrados.InsertPoint(0.289036 , 0.0504014);
	m_Extrados.InsertPoint(0.576000,  0.0350933);
	m_Extrados.InsertPoint(0.736139 , 0.0269428);
	m_Extrados.InsertPoint(1. , 0.);
	m_Extrados.SplineKnots();
	m_Extrados.SplineCurve();

	m_Intrados.m_Color = m_FoilColor;
	m_Intrados.m_Style = m_FoilStyle;
	m_Intrados.m_Width = m_FoilWidth;
	m_Intrados.m_iCtrlPoints = 0;
	m_Intrados.m_iRes = 30;
	m_Intrados.InsertPoint(0. , 0.);
	m_Intrados.InsertPoint(0. , -0.00774066);
	m_Intrados.InsertPoint(0.0306026, -0.0343829);
	m_Intrados.InsertPoint(0.289036 , -0.0504014);
	m_Intrados.InsertPoint(0.576000,  -0.0350933);
	m_Intrados.InsertPoint(0.736139 , -0.0269428);
	m_Intrados.InsertPoint(1. , 0.);
	m_Intrados.SplineKnots();
	m_Intrados.SplineCurve();

	CompMidLine();
	m_OutPoints = m_Extrados.m_iRes + m_Intrados.m_iRes;

	return true;
}


bool CSF::CompMidLine()
{
	double x, yex, yin;
	m_fThickness = 0.0;
	m_fCamber    = 0.0;
	m_fxCambMax  = 0.0;
	m_fxThickMax = 0.0;

	m_rpMid[0].x   = 0.0;
	m_rpMid[0].y   = 0.0;
	m_rpMid[100].x = 1.0;
	m_rpMid[100].y = 0.0;
	// use 0.01 step;
	double step = 0.001;
	for (int k=0; k<=1000; k++)
	{
		x = k*step;
		yex = m_Extrados.GetY(x);
		yin = m_Intrados.GetY(x);
		m_rpMid[k].x = x;
		m_rpMid[k].y = (yex+yin)/2.0;
		if(fabs(yex-yin)>m_fThickness)
		{
			m_fThickness = fabs(yex-yin);
			m_fxThickMax = x;
		}
		if(fabs(m_rpMid[k].y)>fabs(m_fCamber))
		{
			m_fCamber = m_rpMid[k].y;
			m_fxCambMax = x;
		}	
	}
	return true;
}







void CSF::Copy(CSF* pSF)
{
	m_FoilColor = pSF->m_FoilColor;
	m_FoilStyle = pSF->m_FoilStyle;
	m_FoilWidth = pSF->m_FoilWidth;
	m_Extrados.Copy(&pSF->m_Extrados);
	m_Intrados.Copy(&pSF->m_Intrados);
	m_OutPoints  = pSF->m_OutPoints;
	m_fCamber    = pSF->m_fCamber;
	m_fThickness = pSF->m_fThickness;
	m_fxCambMax  = pSF->m_fxCambMax; 
	m_fxThickMax = pSF->m_fxThickMax;
	m_bSymetric  = pSF->m_bSymetric;
}




void CSF::ExportToBuffer(CFoil *pFoil)
{
	int i;
	int j = m_Extrados.m_iRes;
	int k = m_Intrados.m_iRes;
	for (i=0; i<m_Extrados.m_iRes; i++)
	{
		pFoil->x[j-i-1] = m_Extrados.m_Output[i].x;
		pFoil->y[j-i-1] = m_Extrados.m_Output[i].y;
	}
	for (i=1; i<m_Intrados.m_iRes;i++)
	{
		pFoil->x[i+j-1] = m_Intrados.m_Output[i].x;
		pFoil->y[i+j-1] = m_Intrados.m_Output[i].y;
	}
	pFoil->n = j+k-1;
	memcpy(pFoil->xb, pFoil->x, sizeof(pFoil->x));
	memcpy(pFoil->yb, pFoil->y, sizeof(pFoil->y));
	pFoil->nb = pFoil->n;
    pFoil->setName(m_strFoilName);
}


void CSF::ExportToFile(QTextStream &out)
{
	m_Extrados.Export(out, true);
	m_Intrados.Export(out, false);
}

void CSF::SetViewRect(QRect &rc)
{
	m_Intrados.m_rViewRect = rc;
	m_Extrados.m_rViewRect = rc;
}


bool CSF::Serialize(QDataStream &ar, bool bIsStoring)
{
	float f;
	int ArchiveFormat,k ;

	if(bIsStoring)
	{
		ar << 100306;
		// 100306 : changed to C convention	the number of control points

		WriteCString(ar, m_strFoilName);
		WriteCOLORREF(ar, m_FoilColor);
		ar << m_FoilStyle << m_FoilWidth;

		ar << m_Extrados.m_iCtrlPoints;
		ar << m_Extrados.m_iDegree;
		for (k=0; k<m_Extrados.m_iCtrlPoints;k++)
		{
			ar << (float)m_Extrados.m_Input[k].x << (float)m_Extrados.m_Input[k].y;
		}

		ar << m_Intrados.m_iCtrlPoints;
		ar << m_Intrados.m_iDegree;
		for (k=0; k<m_Intrados.m_iCtrlPoints;k++)
		{
			ar << (float)m_Intrados.m_Input[k].x << (float)m_Intrados.m_Input[k].y;
		}

		ar << m_Extrados.m_iRes << m_Intrados.m_iRes;
		if(m_bVisible)    ar<<1; else ar <<0;
		if(m_bOutPoints)  ar<<1; else ar <<0;
		if(m_bCenterLine) ar<<1; else ar <<0;

		m_bModified = false;
		return true;
	}
	else
	{
		ar >> ArchiveFormat;
		if(ArchiveFormat < 100000 || ArchiveFormat > 110000) return false;
		ReadCString(ar, m_strFoilName);
		m_strFoilName = QObject::tr("Spline Foil");
		ReadCOLORREF(ar, m_FoilColor);
		ar >>m_FoilStyle >> m_FoilWidth;

		ar >> m_Extrados.m_iCtrlPoints;
		ar >> m_Extrados.m_iDegree;
		for (k=0; k<m_Extrados.m_iCtrlPoints;k++)
		{
			ar >> f; m_Extrados.m_Input[k].x =f;
			ar >> f; m_Extrados.m_Input[k].y =f;
		}
		if(ArchiveFormat<100306)
		{
			ar >> f; m_Extrados.m_Input[k].x =f;
			ar >> f; m_Extrados.m_Input[k].y =f;
		}

		ar >> m_Intrados.m_iCtrlPoints;
		ar >> m_Intrados.m_iDegree;
		for (k=0; k<m_Intrados.m_iCtrlPoints;k++)
		{
			ar >> f; m_Intrados.m_Input[k].x =f;
			ar >> f; m_Intrados.m_Input[k].y =f;
		}
		if(ArchiveFormat<100306)
		{
			ar >> f; m_Extrados.m_Input[k].x =f;
			ar >> f; m_Extrados.m_Input[k].y =f;
		}
		ar >> m_Extrados.m_iRes >> m_Intrados.m_iRes;
		ar >> k;
		if(k !=0 && k !=1)			return false;

		if(k) m_bVisible = true; else m_bVisible = false;

		ar >> k;
		if(k !=0 && k !=1)			return false;

		if(k) m_bOutPoints = true; else m_bOutPoints = false;

		ar >> k;
		if(k !=0 && k !=1)			return false;

		if(k) m_bCenterLine = true; else m_bCenterLine = false;


		m_Extrados.SplineKnots();
		m_Intrados.SplineKnots();

		Update(true);
		Update(false);
		m_bModified = false;
		return true;
	}
}


void CSF::Update(bool bExtrados)
{
	if(bExtrados)
	{
		m_Extrados.SplineCurve();
	}
	else
	{
		m_Intrados.SplineCurve();
	}
	CompMidLine();
	m_OutPoints = m_Extrados.m_iRes + m_Intrados.m_iRes;
}


void CSF::UpdateKnots()
{
	m_Extrados.SplineKnots();
	m_Intrados.SplineKnots();
}


void CSF::UpdateSelected(double x, double y)
{
	int i;
	for (i=0; i<m_Extrados.m_iCtrlPoints;i++)
	{
		if (m_Extrados.m_iSelect == i)
		{
			m_Extrados.m_Input[i].x = x;
			m_Extrados.m_Input[i].y = y;
			Update(true);
			return;
		}
	}
	for (i=0; i<m_Intrados.m_iCtrlPoints;i++)
	{
		if (m_Intrados.m_iSelect == i)
		{
			m_Intrados.m_Input[i].x = x;
			m_Intrados.m_Input[i].y = y;
			Update(false);
			return;
		}
	}
}



void CSF::DrawCtrlPoints(QPainter &painter, double scalex, double scaley, QPoint Offset)
{
	m_Extrados.DrawCtrlPoints(painter, scalex, scaley, Offset);
	m_Intrados.DrawCtrlPoints(painter, scalex, scaley, Offset);
}


void CSF::DrawOutPoints(QPainter & painter, double scalex, double scaley, QPoint Offset)
{
	m_Extrados.DrawOutputPoints(painter, scalex, scaley, Offset);
	m_Intrados.DrawOutputPoints(painter, scalex, scaley, Offset);
}


void CSF::DrawFoil(QPainter &painter, double scalex, double scaley, QPoint Offset)
{
	m_Extrados.DrawSpline(painter, scalex, scaley, Offset);
	m_Intrados.DrawSpline(painter, scalex, scaley, Offset);
}


void CSF::DrawMidLine(QPainter &painter, double scalex, double scaley, QPoint Offset)
{
	painter.save();
	int k;
	QPoint From, To;
	
	QPen MidPen(m_FoilColor);
	MidPen.setStyle(Qt::DashLine);
	MidPen.setWidth(1);
	painter.setPen(MidPen);

	From = QPoint((int)(m_rpMid[0].x*scalex) + Offset.x(), (int)(-m_rpMid[0].y*scaley) + Offset.y());

	for (k=1; k<=100; k++)
	{
		To.rx() = (int)( m_rpMid[k*10].x*scalex) + Offset.x();
		To.ry() = (int)(-m_rpMid[k*10].y*scaley) + Offset.y();
		painter.drawLine(From, To);
		From = To;
	}
	painter.restore();
}
