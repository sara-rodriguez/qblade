/****************************************************************************

	FoilGeomDlg Class
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


#ifndef FOILGEOMDLG_H
#define FOILGEOMDLG_H

#include <QDialog>
#include <QSlider>
#include "../Misc/NumberEdit.h"
#include "../Objects/Foil.h"

class FoilGeomDlg : public QDialog
{
	Q_OBJECT
	friend class QAFoil;
	friend class QXDirect;
	friend class MainFrame;

public:
	FoilGeomDlg();
	void InitDialog();

private slots:
	void OnRestore();
	void OnOK();

	void OnCamberSlide(int pos);
	void OnXCamberSlide(int pos);
	void OnThickSlide(int pos);
	void OnXThickSlide(int pos);
	void OnCamber();
	void OnXCamber();
	void OnThickness();
	void OnXThickness();


private:
	void keyPressEvent(QKeyEvent *event);
	void SetupLayout();
	void Apply();

private:
	QSlider	*m_pctrlCamberSlide, *m_pctrlThickSlide, *m_pctrlXThickSlide, *m_pctrlXCamberSlide;
	NumberEdit *m_pctrlXCamber;
	NumberEdit	*m_pctrlXThickness;
	NumberEdit	*m_pctrlThickness;
	NumberEdit	*m_pctrlCamber;

	QPushButton *OKButton, *CancelButton, *RestoreButton;


private:
	static void* s_pXFoil;

	double m_fCamber;
	double m_fThickness;
	double m_fXCamber;
	double m_fXThickness;
	CFoil* m_pBufferFoil;
	CFoil* m_pMemFoil;

	void *m_pXDirect;
	void *m_pAFoil;

	bool  m_bApplied,m_bAppliedX, m_bModified;

};

#endif // FOILGEOMDLG_H
