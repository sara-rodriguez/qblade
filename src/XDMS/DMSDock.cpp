/****************************************************************************

	DMSDock Class
		Copyright (C) 2013 Juliane Wendler

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

#include "DMSDock.h"
#include "../Globals.h"
#include "DMS.h"


DMSDock::DMSDock (const QString & title, QMainWindow * parent, Qt::WindowFlags flags)
	: ScrolledDock (title, parent, flags)
{
	g_mainFrame->m_pDMS = new QDMS(g_mainFrame);
	QDMS *pDMS = (QDMS*) g_mainFrame->m_pDMS;
	pDMS->SetupLayout();
	pDMS->Connect();
	pDMS->setAttribute(Qt::WA_DeleteOnClose, false);

	m_contentVBox->addWidget(pDMS->mainWidget);
	setVisible(false);

	addScrolledDock(Qt::LeftDockWidgetArea , parent);

    connect(this,SIGNAL(resized()),pDMS,SLOT(OnResize()));
}

