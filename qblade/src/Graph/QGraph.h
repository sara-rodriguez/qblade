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

#ifndef QGraph_H
#define QGraph_H

class QSettings;

#include "Graph.h"


class QGraph : public Graph
{
public: 
	QGraph();
	virtual ~QGraph();

	void DrawGraph(QRect const & rect, QPainter &painter);
	void DrawGraph(QPainter &painter);
	void DrawAxes(QPainter &painter);
	void DrawCurve(int nIndex, QPainter &painter);
	void DrawLegend(QPainter &painter, QPoint &Place, QFont &LegendFont, QColor &LegendColor);
	void DrawTitles(QPainter &painter);
	void DrawXMinGrid(QPainter &painter);
	void DrawYMinGrid(QPainter &painter);
	void DrawXMajGrid(QPainter &painter);
	void DrawYMajGrid(QPainter &painter);
	void DrawXTicks(QPainter &painter);
	void DrawYTicks(QPainter &painter);
	void ExpFormat(double &f, int &exp);
	void ExportToFile(QFile &XFile, int FileType);
	void Highlight(QPainter &painter, CCurve *pCurve, int ref);

	void LoadSettings(QSettings *pSettings);
	void SaveSettings(QSettings *pSettings);
	QPoint GetOffset();


public:
	void *m_pParent;
	bool m_bHighlightPoint;
};

#endif
