/****************************************************************************

	LineCbBox Class
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


#ifndef LINECBBOX_H
#define LINECBBOX_H

#include <QComboBox>

class LineCbBox : public QComboBox
{
	Q_OBJECT


public:
	LineCbBox(void *pParent = NULL);
	LineCbBox(int const &style, int const &width, QColor const &color);
	void paintEvent (QPaintEvent *event);
	void SetLine(int const &style, int const &width, QColor const &color);

private:
	int m_Style, m_Width;
	QColor m_Color;
};

#endif // LINECBBOX_H
