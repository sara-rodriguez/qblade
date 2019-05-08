/****************************************************************************

    AboutBEM Class
        Copyright (C) 2010 David Marten david.marten@tu-berlin.de

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

#ifndef ABOUTBEM_H
#define ABOUTBEM_H

#include <QDialog>
#include <QLabel>

class AboutBEM : public QDialog
{
    Q_OBJECT

public:
        explicit AboutBEM(QWidget *parent = NULL);
        virtual ~AboutBEM();


private:
        void SetupLayout();

        QPushButton *OKButton;
        QLabel *m_IconQ5;
        QLabel *m_Rotor;
        QLabel *m_HFI;


};

#endif // AboutQ5_H

