/****************************************************************************

    Reference Foil Class
	Copyright (C) 2003-2009 Andre Deperrois adeperrois@xflr5.com

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

#ifndef CFOIL_H
#define CFOIL_H

#include <QColor>
#include <QPainter>
#include <QPoint>
#include "../Params.h"
#include "CVector.h"
#include <QTextStream>
#include "../StorableObject.h"


class CFoil : public StorableObject
{
public:
	static CFoil* newBySerialize ();
	CFoil(QString name = "< no name >");
	
	void serialize();  // override from StorableObject
	static QStringList prepareMissingObjectMessage();	
	
	void GetFoilName(QString &FoilName);
	int IsPoint(CVector const &Real);
	void DrawFoil(QPainter &painter, double const &alpha, double const &scalex, double const &scaley, QPoint const &Offset);
	void DrawPoints(QPainter &painter, double const &scalex, double const &scaley, QPoint const &Offset);
	void DrawMidLine(QPainter &painter, double const &scalex, double const &scaley, QPoint const &Offset);

	void GetLowerY(double x, double &y, double &normx, double &normy);
	void GetUpperY(double x, double &y, double &normx, double &normy);

	double DeRotate();
	double GetBaseUpperY(double const &x);
	double GetBaseLowerY(double const &x);
	double GetMidY(double const &x);
	double GetLowerY(double x);
	double GetUpperY(double x);
	double GetCamber(double const &x);
	double GetCamberAngle(double const &x);
	double GetCamberSlope(double const &x);
	double GetLength();
	double GetArea();
	double GetTopSlope(double const &x);
	double GetBotSlope(double const &x);
	double NormalizeGeometry();
	bool CompMidLine(bool bParams);

	bool ExportFoil(QTextStream &out);
	bool InitFoil();
	bool SetCamber(double f, double u);

	void CopyFoil(CFoil *pSrcFoil);
	void SetFlap();
	void SetTEFlap();
	void SetLEFlap();
	void SetNaca009();
	void SetLEFlapData(bool bFlap, double xhinge, double yhinge, double angle);
	void SetTEFlapData(bool bFlap, double xhinge, double yhinge, double angle);

	bool Intersect(CVector const &A, CVector const &B, CVector const &C, CVector const &D, CVector *M);
	bool IsBetween(int f, int f1, int f2);
	bool IsBetween(int f, double f1, double f2);

	QString m_FoilDescription;	// a description
	QColor m_FoilColor;	//... and its color

	bool m_bVisible;	//true if the foil is to be displayed
	bool m_bCenterLine;	//true if the foil mid camber line is to be displayed
	bool m_bPoints;		//true if the foil's panels are to be displayed
	int m_iExt, m_iInt;	//the number of points on the upper and lower surfaces of the current foil
	int m_iBaseExt, m_iBaseInt;	//the number of points on the upper and lower surfaces of the base foil

	int m_nFoilStyle;
	int m_nFoilWidth;
	int m_iHighLight;

	int n;				// the number of points of the current foil
	double  x[IBX],  y[IBX];	// the point coordinates of the current foil
	double nx[IBX], ny[IBX];	// the normal vector coordinates of the current foil's points

	//Base geometry;
	int nb;				// the number of points of the base foil
	double  xb[IBX],  yb[IBX];	// the point coordinates of the base foil

	double m_fCamber, m_fXCamber;	//foil camber and max camber position
	double m_fThickness, m_fXThickness;// foil thickness and thickness position

	double m_Gap;			// trailing edge gap
	CVector m_LE, m_TE;		// leading edge and trailing edge points

	CVector m_rpBaseMid[1001];	//base mid camber line points
	CVector m_BaseExtrados[IQX];	//base upper surface points
	CVector m_BaseIntrados[IQX];	//base lower surface points

	CVector m_rpMid[1001];		//mid camber line points
	CVector m_rpExtrados[IQX];	//upper surface points
	CVector m_rpIntrados[IQX];	//lower surface points

	// Trailing edge flap  data
	bool m_bTEFlap;
	double m_TEFlapAngle;
	double m_TEXHinge, m_TEYHinge;
	
	// Leading edge flap  data
	bool m_bLEFlap;
	double m_LEFlapAngle;
	double m_LEXHinge, m_LEYHinge;

    ///new code DM///
    bool m_bisCircular;
    ///end new code DM///
};

extern CFoil *g_pCurFoil;

#endif
