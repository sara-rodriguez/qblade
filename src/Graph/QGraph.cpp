/****************************************************************************

	QGraph Classes
	Copyright (C) 2008-2010 Andre Deperrois adeperrois@xflr5.com

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

#include "QGraph.h"

#include <QPainter>
#include <QTextStream>
#include <QSettings>
#include <QDebug>
#include <QRect>
#include <limits>
#include <QDate>
#include <QTime>
#include "../MainFrame.h"

#include "../Globals.h"


QGraph::QGraph()
{
	m_bHighlightPoint = false;

	m_iMargin = 10;

	m_rCltRect.setRect(0,0, 200, 300);

	SetDefaults();
}


QGraph::~QGraph()
{
	DeleteCurves();
}


void QGraph::DrawGraph(QRect const &rect, QPainter &painter)
{
	m_rCltRect = rect;
	DrawGraph(painter);
}


void QGraph::DrawGraph(QPainter &painter)
{
	static QColor color;
	painter.save();

//	Paint background
//	QBrush bg(m_BkColor);
//	painter.setBackground(bg);

//	Draw Border
	if(m_bBorder) color = m_BorderColor;
	else          color = m_BkColor;
	QPen BorderPen(color);
	BorderPen.setStyle(GetStyle(m_BorderStyle));
	BorderPen.setWidth(m_BorderWidth);

	painter.setPen(BorderPen);
	painter.fillRect(m_rCltRect, m_BkColor);
	painter.drawRect(m_rCltRect);
	Init();

	painter.setClipRect(m_rCltRect);

	painter.setBackgroundMode(Qt::TransparentMode);

	if(m_bXMinGrid) DrawXMinGrid(painter);
	if(m_bYMinGrid) DrawYMinGrid(painter);
	if(m_bXMajGrid) DrawXMajGrid(painter);
	if(m_bYMajGrid) DrawYMajGrid(painter);

	DrawAxes(painter);

	DrawXTicks(painter);

	DrawYTicks(painter);

	for (int nc=0; nc < m_oaCurves.size(); nc++)	DrawCurve(nc,painter);

	DrawTitles(painter);

	painter.setClipping(false);
	painter.restore();
}

void QGraph::DrawCurve(int nIndex, QPainter &painter)
{
	painter.save();
	static double scaley;
	static int i, ptside;
	static QPoint From, To, Min, Max;
	static QRect rViewRect;

	ptside = 2;
	CCurve* pCurve = GetCurve(nIndex);

	scaley = m_scaley;

	QBrush FillBrush(m_BkColor);
	painter.setBrush(FillBrush);

	QPen CurvePen(pCurve->GetColor());
	CurvePen.setStyle(GetStyle(pCurve->GetStyle()));
	CurvePen.setWidth((int)pCurve->GetWidth());
	painter.setPen(CurvePen);

	Min.setX(int(xmin/m_scalex) +m_ptoffset.x());
	Min.setY(int(ymin/scaley) +m_ptoffset.y());
	Max.setX(int(xmax/m_scalex) +m_ptoffset.x());
	Max.setY(int(ymax/scaley) +m_ptoffset.y());
	rViewRect.setTopLeft(Min);
	rViewRect.setBottomRight(Max);

	if(pCurve->n>=1)
	{
		From.setX(int(pCurve->x[0]/m_scalex+m_ptoffset.x()));
		From.setY(int(pCurve->y[0]/scaley  +m_ptoffset.y()));

		if(pCurve->IsVisible())
		{
			for (i=1; i<pCurve->n;i++)
			{
				To.setX(int(pCurve->x[i]/m_scalex+m_ptoffset.x()));
				To.setY(int(pCurve->y[i]/scaley  +m_ptoffset.y()));
				painter.drawLine(From, To);
	
				From = To;
			}
		}

		if(pCurve->PointsVisible())
		{
			for (i=0; i<pCurve->n;i++)
			{
				if(pCurve->GetSelected() !=i)
					painter.drawRect(int(pCurve->x[i]/m_scalex+m_ptoffset.x())-ptside,
									 int(pCurve->y[i]/  scaley+m_ptoffset.y())-ptside,
									 2*ptside,2*ptside);
			}
		}
	}

	if(m_bHighlightPoint)
	{
		int point = pCurve->GetSelected();
		if(point>=0)
		{
			//highlight
			QColor HighColor(0,40, 150);
			CurvePen.setWidth((int)pCurve->GetWidth());
			CurvePen.setColor(HighColor);
			painter.setPen(CurvePen);
			To.setX(int(pCurve->x[point]/m_scalex+m_ptoffset.x()));
			To.setY(int(pCurve->y[point]/scaley  +m_ptoffset.y()));
			painter.drawRect(To.x()-ptside-1,To.y()-ptside-1, 2*(ptside+1),2*(ptside+1));
		}
	}
	painter.restore();
}


void QGraph::DrawAxes(QPainter &painter)
{	
	static double xp, yp, scaley;
	static QPen AxesPen;
	scaley = m_scaley;
	painter.save();

	AxesPen.setColor(m_AxisColor);
	AxesPen.setStyle(GetStyle(m_AxisStyle));
	AxesPen.setWidth(m_AxisWidth);
	painter.setPen(AxesPen);

	//vertical axis
	if(xo>=xmin && xo<=xmax) xp = xo;
	else if(xo>xmax)         xp = xmax;
	else                     xp = xmin;

	painter.drawLine((int)(xp/m_scalex) + m_ptoffset.x(), (int)(ymin/scaley) + m_ptoffset.y(),
					 (int)(xp/m_scalex) + m_ptoffset.x(), (int)(ymax/scaley) + m_ptoffset.y());

	//horizontal axis
	if(yo>=ymin && yo<=ymax)	yp = yo;
	else if(yo>ymax)		yp = ymax;
	else				yp = ymin;


	painter.drawLine((int)(xmin/m_scalex) +m_ptoffset.x(), (int)(yp/scaley) + m_ptoffset.y(),
					 (int)(xmax/m_scalex) +m_ptoffset.x(), (int)( yp/scaley) + m_ptoffset.y());

	painter.restore();
}


void QGraph::DrawTitles(QPainter &painter)
{
	//draws the x & y axis name
	static double scaley;
	static int XPosXTitle, YPosXTitle, XPosYTitle, YPosYTitle;
	static double xp, yp;

	scaley = m_scaley;
	painter.save();
	XPosXTitle = -35;
	YPosXTitle = -10;
	XPosYTitle = -5;
	YPosYTitle =  15;

	if(xo>=xmin && xo<=xmax) xp = xo;
	else if(xo>xmax)         xp = xmax;
	else					 xp = xmin;

	if(yo>=ymin && yo<=ymax) yp = yo;
	else if(yo>ymax)         yp = ymax;
	else                     yp = ymin;

	painter.setFont(m_TitleLogFont);

	QPen TitlePen(m_TitleColor);
	painter.setPen(TitlePen);

	painter.drawText(  (int)(xmax/m_scalex) + m_ptoffset.x() + XPosXTitle,
					   (int)(yp  /scaley)   + m_ptoffset.y() + YPosXTitle, m_XTitle);

	painter.drawText(  m_ptoffset.x() + (int)(xp/m_scalex)   + XPosYTitle,
					   m_rCltRect.top() + m_iMargin - YPosYTitle, m_YTitle);
	painter.restore();
}


void QGraph::DrawXTicks(QPainter &painter)
{
	static double main, scaley, xt, yp;
	static int exp, TickSize, height, yExpOff/*, xMainOff*/, nx;

	if(fabs(xunit)<0.00000001) return;
	if(fabs(xmax-xmin)/xunit>30.0) return;

	scaley = m_scaley;
	painter.save();
	QString strLabel, strLabelExp;

	QFontMetrics fm(m_LabelLogFont);

	painter.setFont(m_LabelLogFont);


	TickSize = 5;
	height  = fm.height()/2;
	yExpOff = height/2;

	QPen LabelPen(m_AxisColor);

	LabelPen.setStyle(GetStyle(m_AxisStyle));
	LabelPen.setWidth(m_AxisWidth);
	painter.setPen(LabelPen);
	xt = xo-(xo-xmin);//one tick at the origin
	nx = (int)((xo-xmin)/xunit);
	xt = xo - nx*xunit;


	if(yo>=ymin && yo<=ymax) yp = yo;
	else if(yo>ymax)         yp = ymax;
	else                     yp = ymin;

	while(xt<=xmax*1.0001)
	{
		//Draw ticks
		if(xt>=xmin)
		{
			painter.setPen(LabelPen);
			painter.drawLine(int(xt/m_scalex) + m_ptoffset.x(),int(yp/scaley) +TickSize + m_ptoffset.y(),
							 int(xt/m_scalex) + m_ptoffset.x(),int(yp/scaley)           + m_ptoffset.y());
			painter.setPen(m_LabelColor);

			if(exp_x>=4 || exp_x<=-4)
			{
				main = xt;
				ExpFormat(main, exp);

				strLabel = QString("%1 10").arg(main,5,'f',2);
				painter.drawText(int(xt/m_scalex) - fm.width(strLabel)/2  +m_ptoffset.x(),
								 int(yp/scaley)   + TickSize*2 +height    +m_ptoffset.y(),
								 strLabel);
				strLabelExp = QString("%1").arg(exp);

				painter.drawText(int(xt/m_scalex) + fm.width(strLabel)/2       +m_ptoffset.x(),
								 int(yp/scaley)   + TickSize*2 +height-yExpOff +m_ptoffset.y(),
								 strLabelExp);
			}
			else
			{
				if(exp_x>0)        strLabel = QString("%1").arg(xt,0,'f',0);
				else if (exp_x>=-1) strLabel = QString("%1").arg(xt,6,'f',1);
				else if (exp_x>=-2) strLabel = QString("%1").arg(xt,6,'f',2);
				else if (exp_x>=-3) strLabel = QString("%1").arg(xt,6,'f',3);
				painter.drawText((int)(xt/m_scalex) - fm.width(strLabel)/2 + m_ptoffset.x(),
								 (int)(yp/scaley)   + TickSize*2 +height   + m_ptoffset.y(),
								 strLabel);
			}
		}
		xt += xunit ;
	}
	painter.restore();
}




void QGraph::DrawYTicks(QPainter &painter)
{
	static double scaley, xp, main, yt;
	static int TickSize, fmheight, fmheight4, exp;
	if(fabs(xunit)<0.00000001) return;
	if(fabs(ymax-ymin)/yunit>30.0) return;
	scaley = m_scaley;
	painter.save();
	QString strLabel, strLabelExp;
	exp = 0.0;
	QFontMetrics fm(m_LabelLogFont);
	painter.setFont(m_LabelLogFont);

	fmheight  = fm.height();
	fmheight4 = (int)(fmheight/4);

	TickSize = 5;

	QPen LabelPen(m_AxisColor);
	LabelPen.setStyle(GetStyle(m_AxisStyle));
	LabelPen.setWidth(m_AxisWidth);


	if(xo>=xmin && xo<=xmax) xp = xo;
	else if(xo>xmax)         xp = xmax;
	else                     xp = xmin;

	yt = yo-int((yo-ymin)*1.0001/yunit)*yunit;//one tick at the origin
	
		while(yt<=ymax*1.0001)
		{
			//Draw ticks
			if(yt>=ymin)
			{
				painter.setPen(LabelPen);
				painter.drawLine((int)(xp/m_scalex)          +m_ptoffset.x(), (int)(yt/scaley) +m_ptoffset.y(),
								 (int)(xp/m_scalex)-TickSize +m_ptoffset.x(), (int)(yt/scaley) +m_ptoffset.y());
	
				painter.setPen(m_LabelColor);
	
				if(exp_y>=3 || exp_y<=-3)
				{
					main = yt;
					ExpFormat (main, exp);
	
					strLabel    = QString("%1 10").arg(main,5,'f',2);
					strLabelExp = QString("%1").arg(exp);
	
					painter.drawText((int)(xp/m_scalex) - fm.width(strLabel)-TickSize*3 + m_ptoffset.x(),
									 (int)(yt/scaley)   + fmheight4                     + m_ptoffset.y(),
									 strLabel);
	
					if(exp_y>=3)
					{
						painter.drawText(int(xp/m_scalex) -TickSize*3 + m_ptoffset.x(),
										 int(yt/scaley)               + m_ptoffset.y(),
										 strLabelExp);
					}
					else
					{
						painter.drawText(int(xp/m_scalex) -TickSize*3 + 2 + m_ptoffset.x(),
										 int(yt/scaley)                   + m_ptoffset.y(),
										 strLabelExp);
					}
				}
				else
				{
					if(exp_y>=0)        strLabel = QString("%1").arg(yt,6,'f',0);
					else if (exp_y>=-1) strLabel = QString("%1").arg(yt,6,'f',1);
					else if (exp_y>=-2) strLabel = QString("%1").arg(yt,6,'f',2);
					else if (exp_y>=-3) strLabel = QString("%1").arg(yt,6,'f',3);
					
					painter.drawText((int)(xp/m_scalex) - fm.width(strLabel)-TickSize*2 +m_ptoffset.x(),
									 (int)(yt/scaley)   + fmheight4 +m_ptoffset.y(),
									 strLabel);
				}
			}
			yt += yunit ;
		}
	
	painter.restore();
}


void QGraph::DrawXMajGrid(QPainter &painter)
{
	double scaley = m_scaley;
	if(fabs(xunit)<0.00000001)     return;
	if(fabs(xmax-xmin)/xunit>30.0) return;

	painter.save();
	int YMin, YMax;

	QPen GridPen(m_XMajClr);

	GridPen.setStyle(GetStyle(m_XMajStyle));
	GridPen.setWidth(m_XMajWidth);
	painter.setPen(GridPen);

	YMin = (int)(ymin/scaley) + m_ptoffset.y();
	YMax = (int)(ymax/scaley) + m_ptoffset.y();


	double xt = xo-int((xo-xmin)*1.0001/xunit)*xunit;//one tick at the origin
	while(xt<=xmax*1.001)
	{
		if(xt>=xmin)
		{
			painter.drawLine(int(xt/m_scalex) + m_ptoffset.x(), YMin, int(xt/m_scalex) + m_ptoffset.x(), YMax);
		}
		xt += xunit ;
	}
	painter.restore();
}


void QGraph::DrawYMajGrid(QPainter &painter)
{
	double scaley = m_scaley;
	if(fabs(yunit)<0.00000001) return;
	if(fabs(ymax-ymin)/yunit>30.0) return;

	painter.save();
	int width;
	if(m_YMajWidth<=1) width = 1;

	QPen GridPen(m_YMajClr);

	GridPen.setStyle(GetStyle(m_YMajStyle));
	GridPen.setWidth(width);
	painter.setPen(GridPen);

	double yt = yo-int((yo-ymin)*1.0001/yunit)*yunit;//one tick at the origin

	int XMin = qMax((int)(xmin/m_scalex + m_ptoffset.x()), m_rCltRect.left());
	int XMax = qMin((int)(xmax/m_scalex + m_ptoffset.x()), m_rCltRect.right());

	while(yt<=ymax*1.0001)
	{
		if(yt>=ymin)
		{
			painter.drawLine(XMin, (int)(yt/scaley)   + m_ptoffset.y(), XMax, (int)(yt/scaley)   + m_ptoffset.y());
		}
		yt += yunit ;
	}
	painter.restore();
}

void QGraph::DrawXMinGrid(QPainter &painter)
{
	double scaley = m_scaley;
	if(fabs(xunit)<0.00000001) return;
	if(fabs(m_XMinorUnit)<0.00000001) return;
	if(fabs(xmax-xmin)/xunit>30.0) return;
	if(fabs(xmax-xmin)/m_XMinorUnit>100.0) return;
	int YMin, YMax;

	painter.save();
	QPen GridPen(m_XMinClr);
	GridPen.setStyle(GetStyle(m_XMinStyle));
	GridPen.setWidth(m_XMinWidth);
	painter.setPen(GridPen);


	YMin = (int)(ymin/scaley)+ m_ptoffset.y();
	YMax = (int)(ymax/scaley)+ m_ptoffset.y();

	double xDelta = m_XMinorUnit;
	double xt = xo-int((xo-xmin)*1.0001/xDelta)*xDelta;//one tick at the origin


	while(xt<=xmax*1.001)
	{
		if(xt>=xmin)
		{
			painter.drawLine(int(xt/m_scalex) + m_ptoffset.x(), YMin, int(xt/m_scalex) + m_ptoffset.x(), YMax);
		}
		xt += xDelta;
	}
	painter.restore();
}

void QGraph::DrawYMinGrid(QPainter &painter)
{
	double scaley = m_scaley;
	if(fabs(yunit)<0.00000001) return;
	if(fabs(m_YMinorUnit)<0.00000001) return;
	if(fabs(ymax-ymin)/yunit>30.0) return;
	if(fabs(ymax-ymin)/m_YMinorUnit>100.0) return;

	painter.save();
	QPen GridPen(m_YMinClr);
	GridPen.setStyle(GetStyle(m_YMinStyle));
	GridPen.setWidth(m_YMinWidth);
	painter.setPen(GridPen);

	double yDelta = m_YMinorUnit;
	double yt = yo-int((yo-ymin)*1.0001/yDelta)*yDelta;//one tick at the origin
	int XMin = qMax((int)(xmin/m_scalex + m_ptoffset.x()), m_rCltRect.left());
	int XMax = qMin((int)(xmax/m_scalex + m_ptoffset.x()), m_rCltRect.right());

	while(yt<=ymax*1.0001)
	{
		if(yt>=ymin)
		{
			painter.drawLine(XMin, (int)(yt/scaley)   + m_ptoffset.y(), XMax, (int)(yt/scaley)   + m_ptoffset.y());
		}
		yt += yDelta ;
	}
	painter.restore();
}


void QGraph::DrawLegend(QPainter &painter, QPoint &Place, QFont &LegendFont, QColor &LegendColor)
{
	painter.save();
	int LegendSize, ypos;
	QString strong;

	LegendSize = 30;
	ypos = 12;

	painter.setFont(LegendFont);

	CCurve* pCurve;

	QPen TextPen(LegendColor);
	QPen LegendPen(Qt::gray);

	int npos = 0;
	for (int nc=0; nc< m_oaCurves.size(); nc++)
	{
		pCurve = (CCurve*) m_oaCurves[nc];
		if(pCurve->IsVisible())
		{
			pCurve->GetTitle(strong);
			if(pCurve->n>0 && strong.length())//is there anything to draw ?
			{

				LegendPen.setColor(pCurve->GetColor());
				LegendPen.setStyle(GetStyle(pCurve->GetStyle()));
				LegendPen.setWidth(pCurve->GetWidth());

				painter.setPen(LegendPen);

				painter.drawLine(Place.x(),                     Place.y() + ypos*npos,
								 Place.x() + (int)(LegendSize), Place.y() + ypos*npos);

				painter.setPen(TextPen);
				painter.drawText(Place.x() + (int)(1.5*LegendSize),    Place.y()  + ypos*npos+(int)(ypos/2),
								 strong);

				npos++;
			}
		}
	}

	painter.restore();

}

void QGraph::ExpFormat(double &f, int &exp) {
	if (fabs(f) < 0.00000000001) {
		f = 0;
		exp = 0;
	} else {
		double f1 = fabs(f);
		if (f1 < 1)
			exp = (int)log10(f1)-1;
		else
			exp = (int)log10(f1);
		
		f = f/pow(10.0,exp);
	}
}


void QGraph::ExportToFile(QFile &XFile, int FileType)
{
	int i,j, maxpoints;
	CCurve *pCurve;
	QString strong, xTitle;
	QTextStream out(&XFile);
    QDate date = QDate::currentDate();
    QTime time = QTime::currentTime();


	maxpoints = 0;

    out << "Export File Created with "<< g_mainFrame->m_VersionName<<" on "<<date.toString("dd.MM.yyyy") << " at " << time.toString("hh:mm:ss") << "\n" ;


	for(i=0; i<m_oaCurves.size(); i++)
	{
		pCurve = GetCurve(i);
        if(pCurve && pCurve->n > 2)
		{
			maxpoints = qMax(maxpoints,pCurve->n); 
			pCurve->GetTitle(strong);
            out << QString("  %1").arg(strong, 30);
		}
    }

    out<<"\n";

    for(i=0; i<m_oaCurves.size(); i++)
    {
        pCurve = GetCurve(i);
        if(pCurve && pCurve->n > 2)
        {
            maxpoints = qMax(maxpoints,pCurve->n);

            if(FileType==1)	out <<" "<<m_XTitle<<" "<< m_YTitle;
            else            out << m_XTitle<<";"<< m_YTitle << ";";

        }
	}
	out<<"\n"; //end of title line

	for(j=0; j<maxpoints; j++)
	{
		for(i=0; i<m_oaCurves.size(); i++)
		{
			pCurve = GetCurve(i);
            if(pCurve && j<pCurve->n && pCurve->n > 2)
			{
                if(FileType==1)	strong= QString(" %1 %2")
                                                .arg(pCurve->x[j],15,'e',5).arg(pCurve->y[j],15,'e',5);
                else            strong= QString(" %1; %2;")
                                                .arg(pCurve->x[j],15,'e',5).arg(pCurve->y[j],15,'e',5).replace(".",",");
			}
            else if(pCurve->n > 2)
			{
                if(FileType==1)	strong=QString().fill(' ', 32);
                else            strong= ";;";
            }
            if (pCurve->n > 2) out << strong;
		}
		out<<"\n"; //end of data line
	}
	out<<"\n"; //end of file
	XFile.close();
}



QPoint QGraph::GetOffset()
{
	return m_ptoffset;
}

void QGraph::Highlight(QPainter &painter, CCurve *pCurve, int ref)
{
	painter.save();
	int x = int(pCurve->x[ref]/m_scalex)  +m_ptoffset.x();
	int y = int(pCurve->y[ref]/m_scaley)  +m_ptoffset.y();

	QPen HighlightPen(QColor(255,100,100));
	HighlightPen.setWidth(2);
	painter.setPen(HighlightPen);
	QRect r(x-3,y-3,6,6);;
	painter.drawRect(r);
	painter.restore();
}


void QGraph::SaveSettings(QSettings *pSettings)
{
	QFont lgft;
	QColor clr;
	int k,s,w;
	bool ba, bs;
	double f;

	pSettings->beginGroup(m_GraphName);
	{
		//read variables
		clr = GetAxisColor();
		pSettings->setValue("AxisColorRed", clr.red());
		pSettings->setValue("AxisColorGreen",clr.green());
		pSettings->setValue("AxisColorBlue", clr.blue());

		k = GetAxisStyle();
		pSettings->setValue("AxisStyle", k);
		k = GetAxisWidth();
		pSettings->setValue("AxisWidth", k);

		clr = GetTitleColor();
		pSettings->setValue("TitleColorRed", clr.red());
		pSettings->setValue("TitleColorGreen",clr.green());
		pSettings->setValue("TitleColorBlue", clr.blue());

		clr = GetLabelColor();
		pSettings->setValue("LabelColorRed", clr.red());
		pSettings->setValue("LabelColorGreen",clr.green());
		pSettings->setValue("LabelColorBlue", clr.blue());

		GetTitleLogFont(&lgft);
		pSettings->setValue("TitleFontName", lgft.family());
		pSettings->setValue("TitleFontSize", lgft.pointSize());

		GetLabelLogFont(&lgft);
		pSettings->setValue("LabelFontName", lgft.family());
		pSettings->setValue("LabelFontSize", lgft.pointSize());

		GetXMajGrid(bs,clr,s,w);
		pSettings->setValue("XMajGridColorRed", clr.red());
		pSettings->setValue("XMajGridColorGreen",clr.green());
		pSettings->setValue("XMajGridColorBlue", clr.blue());
		pSettings->setValue("XMajGridShow",bs);
		pSettings->setValue("XMajGridStyle",s);
		pSettings->setValue("XMajGridWidth",w);

		GetYMajGrid(bs,clr,s,w);
		pSettings->setValue("YMajGridColorRed", clr.red());
		pSettings->setValue("YMajGridColorGreen",clr.green());
		pSettings->setValue("YMajGridColorBlue", clr.blue());
		pSettings->setValue("YMajGridShow",bs);
		pSettings->setValue("YMajGridStyle",s);
		pSettings->setValue("YMajGridWidth",w);

		GetXMinGrid(bs,ba,clr,s,w,f);
		pSettings->setValue("XMinGridColorRed", clr.red());
		pSettings->setValue("XMinGridColorGreen",clr.green());
		pSettings->setValue("XMinGridColorBlue", clr.blue());
		pSettings->setValue("XMinGridAuto",ba);
		pSettings->setValue("XMinGridShow",bs);
		pSettings->setValue("XMinGridStyle",s);
		pSettings->setValue("XMinGridWidth",w);
		pSettings->setValue("XMinGridUnit",f);

		GetYMinGrid(bs,ba,clr,s,w,f);
		pSettings->setValue("YMinGridColorRed", clr.red());
		pSettings->setValue("YMinGridColorGreen",clr.green());
		pSettings->setValue("YMinGridColorBlue", clr.blue());
		pSettings->setValue("YMinGridAuto",ba);
		pSettings->setValue("YMinGridShow",bs);
		pSettings->setValue("YMinGridStyle",s);
		pSettings->setValue("YMinGridWidth",w);
		pSettings->setValue("YMinGridUnit",f);

		clr = GetBorderColor();
		s   = GetBorderStyle();
		w   = GetBorderWidth();
		pSettings->setValue("BorderColorRed", clr.red());
		pSettings->setValue("BorderColorGreen", clr.green());
		pSettings->setValue("BorderColorBlue", clr.blue());
		pSettings->setValue("BorderStyle", s);
		pSettings->setValue("BorderWidth", w);
		pSettings->setValue("BorderShow", m_bBorder);

		clr = GetBackColor();
		pSettings->setValue("BackColorRed", clr.red());
		pSettings->setValue("BackColorGreen", clr.green());
		pSettings->setValue("BackColorBlue", clr.blue());

		pSettings->setValue("Inverted", m_bYInverted);

		pSettings->setValue("XVariable", m_X);
		pSettings->setValue("YVariable", m_Y);

        ///new code DM///
        pSettings->setValue("MainVar", m_MainVar);
        pSettings->setValue("Param", m_Param);
        ///end new code DM///

        ///new code JW///
        pSettings->setValue("Type", m_Type);
        ///end new code JW///


	}
	pSettings->endGroup();
}


void QGraph::LoadSettings(QSettings *pSettings)
{
	int k;

	QString FontName;
	QFont lgft;
	bool bs, ba;
	int s,w;
	int r,g,b;
	double f;

	pSettings->beginGroup(m_GraphName);
	{
		//read variables
        r = pSettings->value("AxisColorRed",200).toInt();
        g = pSettings->value("AxisColorGreen",200).toInt();
        b = pSettings->value("AxisColorBlue",200).toInt();
		SetAxisColor(QColor(r,g,b));

		k = pSettings->value("AxisStyle",0).toInt();
		SetAxisStyle(k);
		k = pSettings->value("AxisWidth",1).toInt();
		SetAxisWidth(k);

        r = pSettings->value("TitleColorRed",0).toInt();
        g = pSettings->value("TitleColorGreen",0).toInt();
        b = pSettings->value("TitleColorBlue",0).toInt();
		SetTitleColor(QColor(r,g,b));
        r = pSettings->value("LabelColorRed",0).toInt();
        g = pSettings->value("LabelColorGreen",0).toInt();
        b = pSettings->value("LabelColorBlue",0).toInt();
		SetLabelColor(QColor(r,g,b));

		FontName = pSettings->value("TitleFontName","Comic Sans MS").toString();
		lgft.setFamily(FontName);
		lgft.setPointSize(pSettings->value("TitleFontSize",8).toInt());
		SetTitleLogFont(&lgft);

		FontName = pSettings->value("LabelFontName","Comic Sans MS").toString();
		lgft.setFamily(FontName);
		lgft.setPointSize(pSettings->value("LabelFontSize",8).toInt());
		SetLabelLogFont(&lgft);


        r  = pSettings->value("XMajGridColorRed",200).toInt();
        g  = pSettings->value("XMajGridColorGreen",200).toInt();
        b  = pSettings->value("XMajGridColorBlue",200).toInt();
		bs = pSettings->value("XMajGridShow",true).toBool();
		s  = pSettings->value("XMajGridStyle",1).toInt();
		w  = pSettings->value("XMajGridWidth",1).toInt();
		SetXMajGrid(bs,QColor(r,g,b),s,w);

        r  = pSettings->value("YMajGridColorRed",200).toInt();
        g  = pSettings->value("YMajGridColorGreen",200).toInt();
        b  = pSettings->value("YMajGridColorBlue",200).toInt();
		bs = pSettings->value("YMajGridShow",true).toBool();
		s  = pSettings->value("YMajGridStyle",1).toInt();
		w  = pSettings->value("YMajGridWidth",1).toInt();
		SetYMajGrid(bs,QColor(r,g,b),s,w);

        r  = pSettings->value("XMinGridColorRed",200).toInt();
        g  = pSettings->value("XMinGridColorGreen",200).toInt();
        b  = pSettings->value("XMinGridColorBlue",200).toInt();
		ba = pSettings->value("XMinGridAuto",true).toBool();
		bs = pSettings->value("XMinGridShow",false).toBool();
		s  = pSettings->value("XMinGridStyle",2).toInt();
		w  = pSettings->value("XMinGridWidth",1).toInt();
		f  = pSettings->value("XMinGridUnit", 0.01).toDouble();
		SetXMinGrid(bs,ba,QColor(r,g,b),s,w,f);

        r  = pSettings->value("YMinGridColorRed",200).toInt();
        g  = pSettings->value("YMinGridColorGreen",200).toInt();
        b  = pSettings->value("YMinGridColorBlue",200).toInt();
		ba = pSettings->value("YMinGridAuto",true).toBool();
		bs = pSettings->value("YMinGridShow",false).toBool();
		s  = pSettings->value("YMinGridStyle",2).toInt();
		w  = pSettings->value("YMinGridWidth",1).toInt();
		f  = pSettings->value("YMinGridUnit",0.01).toDouble();
		SetYMinGrid(bs,ba,QColor(r,g,b),s,w,f);

        r  = pSettings->value("BorderColorRed",240).toInt();
        g  = pSettings->value("BorderColorGreen",240).toInt();
        b  = pSettings->value("BorderColorBlue",240).toInt();
		s  = pSettings->value("BorderStyle",0).toInt();
		w  = pSettings->value("BorderWidth",2).toInt();
		m_bBorder = pSettings->value("BorderShow", true).toBool();
		SetBorderColor(QColor(r,g,b));
		SetBorderStyle(s);
		SetBorderWidth(w);

		r  = pSettings->value("BackColorRed",0).toInt();
		g  = pSettings->value("BackColorGreen",20).toInt();
		b  = pSettings->value("BackColorBlue",20).toInt();
		SetBkColor(QColor(r,g,b));

		m_bYInverted = pSettings->value("Inverted", false).toBool();

		m_X  = pSettings->value("XVariable",1).toInt();
		m_Y  = pSettings->value("YVariable",0).toInt();

        ///new code DM///
        m_MainVar = pSettings->value("MainVar",0).toInt();
        m_Param = pSettings->value("Param",0).toInt();
        ///end new code DM///

        ///new code JW///
        m_Type = pSettings->value("Type",0).toInt();
        ///end new code JW///

	}
	pSettings->endGroup();
}
