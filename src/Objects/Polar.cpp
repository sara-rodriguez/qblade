/****************************************************************************

    Polar Class
	Copyright (C) 2003 Andre Deperrois adeperrois@xflr5.com

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

#include "Polar.h"

#include <QTextStream>

#include "../Serializer.h"
#include "../MainFrame.h"
#include "OpPoint.h"
#include "../XDirect/XFoil.h"
#include "../Globals.h"
#include <QDate>
#include <QTime>
#include "../Store.h"


CPolar::CPolar(QString name, StorableObject *parent)
    : StorableObject (name, parent)
{
	m_bIsVisible  = true;
	m_bShowPoints = false;
	m_Style = 0;// = PS_SOLID
	m_Width = 1;
	m_ASpec = 0.0;
	m_PolarType = FIXEDSPEEDPOLAR;
	m_ReType = 1;
	m_MaType = 1;
	m_Reynolds = 100000.0;
	m_Mach     = 0.0;
	m_ACrit    = 9.0;
	m_XTop     = 1.0;
	m_XBot     = 1.0;
}

void CPolar::serialize() {
	StorableObject::serialize();
	
	g_serializer.readOrWriteInt (&m_ReType);
	g_serializer.readOrWriteInt (&m_MaType);
	g_serializer.readOrWriteDouble (&m_ASpec);
	g_serializer.readOrWriteDouble (&m_Mach);
	g_serializer.readOrWriteDouble (&m_ACrit);
	g_serializer.readOrWriteDouble (&m_XTop);
	g_serializer.readOrWriteDouble (&m_XBot);
	
	g_serializer.readOrWriteInt (&m_Style);
	g_serializer.readOrWriteInt (&m_Width);
	g_serializer.readOrWriteColor (&m_Color);

	g_serializer.readOrWriteBool (&m_bIsVisible);
	g_serializer.readOrWriteBool (&m_bShowPoints);
		
	g_serializer.readOrWriteDouble (&m_Reynolds);
	int temp = (int) m_PolarType;
	g_serializer.readOrWriteInt (&temp);
	m_PolarType = (enumPolarType) temp;

	g_serializer.readOrWriteDoubleList1D (&m_Alpha);
	g_serializer.readOrWriteDoubleList1D (&m_Cl);
	g_serializer.readOrWriteDoubleList1D (&m_XCp);
	g_serializer.readOrWriteDoubleList1D (&m_Cd);
	g_serializer.readOrWriteDoubleList1D (&m_Cdp);
	g_serializer.readOrWriteDoubleList1D (&m_Cm);
	g_serializer.readOrWriteDoubleList1D (&m_XTr1);
	g_serializer.readOrWriteDoubleList1D (&m_XTr2);
	g_serializer.readOrWriteDoubleList1D (&m_HMom);
	g_serializer.readOrWriteDoubleList1D (&m_Cpmn);
	g_serializer.readOrWriteDoubleList1D (&m_ClCd);
	g_serializer.readOrWriteDoubleList1D (&m_Cl32Cd);
	g_serializer.readOrWriteDoubleList1D (&m_RtCl);
	g_serializer.readOrWriteDoubleList1D (&m_Re);
}




void CPolar::ExportPolar(QTextStream &out, int FileType, bool bDataOnly)
{
	QString Header, strong;
	int j;

	if(!bDataOnly)
	{
		strong = g_mainFrame->m_VersionName + "\n\n";
		out << strong;
		strong =(" Calculated polar for: ");
		strong += getParent()->getName() + "\n\n";
		out << strong;
		strong = QString(" %1 %2").arg(m_ReType).arg(m_MaType);
		if(m_ReType==1) strong += (" Reynolds number fixed       ");
		else if(m_ReType==2) strong += (" Reynolds number ~ 1/sqrt(CL)");
		else if(m_ReType==3) strong += (" Reynolds number ~ 1/CL      ");
		if(m_MaType==1) strong += ("   Mach number fixed         ");
		else if(m_MaType==2) strong += ("   Mach number ~ 1/sqrt(CL)  ");
		else if(m_MaType==3) strong += ("   Mach number ~ 1/CL        ");
		strong +="\n\n";
		out << strong;
		strong=QString((" xtrf =   %1 (top)        %2 (bottom)\n"))
						.arg(m_XTop,0,'f',3).arg(m_XBot,0,'f',3);
		out << strong;

		strong = QString(" Mach = %1     Re = %2 e 6     Ncrit = %3\n\n")
				 .arg(m_Mach,7,'f',3).arg(m_Reynolds/1.e6,9,'f',3).arg(m_ACrit,7,'f',3);
		out << strong;
	}

	if(m_PolarType!=FIXEDAOAPOLAR)
	{
		if(FileType==1) Header = ("  alpha     CL        CD       CDp       Cm    Top Xtr Bot Xtr   Cpmin    Chinge    XCp    \n");
		else            Header = ("alpha,CL,CD,CDp,Cm,Top Xtr,Bot Xtr,Cpmin,Chinge,XCp\n");
		out << Header;
		if(FileType==1)
		{
			Header=QString(" ------- -------- --------- --------- -------- ------- ------- -------- --------- ---------\n");
			out << Header;
		}
		for (j=0; j<m_Alpha.size(); j++)
		{
			if(FileType==1) strong = QString(" %1  %2  %3  %4  %5")
											.arg(m_Alpha[j],7,'f',3)
											.arg(m_Cl[j],7,'f',4)
											.arg(m_Cd[j],8,'f',5)
											.arg(m_Cdp[j],8,'f',5)
											.arg(m_Cm[j],7,'f',4);
			else            strong = QString("%1,%2,%3,%4,%5")
											.arg(m_Alpha[j],7,'f',3)
											.arg(m_Cl[j],7,'f',4)
											.arg(m_Cd[j],8,'f',5)
											.arg(m_Cdp[j],8,'f',5)
											.arg(m_Cm[j],7,'f',4);

			out << strong;
			if(m_XTr1[j]<990.0)
			{
				if(FileType==1) strong=QString("  %1  %2").arg(m_XTr1[j],6,'f',4).arg( m_XTr2[j],6,'f',4);
				else            strong=QString(",%1,%2").arg(m_XTr1[j],6,'f',4).arg( m_XTr2[j],6,'f',4);
				out << strong;
			}
			if(FileType==1) strong=QString("  %1  %2  %3\n").arg(m_Cpmn[j],7,'f',4).arg(m_HMom[j],7,'f',4).arg(m_XCp[j],7,'f',4);
			else            strong=QString(",%1,%2,%3\n").arg(m_Cpmn[j],7,'f',4).arg(m_HMom[j],7,'f',4).arg(m_XCp[j],7,'f',4);
			out << strong;
			}
	}
	else 
	{
		if(FileType==1) Header=QString(("  alpha     Re      CL        CD       CDp       Cm    Top Xtr Bot Xtr   Cpmin    Chinge     XCp    \n"));
		else            Header=QString(("alpha,Re,CL,CD,CDp,Cm,Top Xtr,Bot Xtr,Cpmin,Chinge,XCp\n"));
		out << Header;
		if(FileType==1)
		{
			Header=QString(" ------- -------- -------- --------- --------- -------- ------- ------- -------- --------- ---------\n");
			out << Header;
		}
		for (j=0; j<m_Alpha.size(); j++)
		{
			if(FileType==1) strong=QString(" %1 %2  %3  %4  %5  %6")
											.arg(m_Alpha[j],7,'f',3)
											.arg( m_Re[j],8,'f',0)
											.arg( m_Cl[j],7,'f',4)
											.arg( m_Cd[j],8,'f',5)
											.arg(m_Cdp[j],8,'f',5)
											.arg(m_Cm[j],7,'f',4);
			else            strong=QString(" %1,%2,%3,%4,%5,%6")
											.arg(m_Alpha[j],7,'f',3)
											.arg( m_Re[j],8,'f',0)
											.arg( m_Cl[j],7,'f',4)
											.arg( m_Cd[j],8,'f',5)
											.arg(m_Cdp[j],8,'f',5)
											.arg(m_Cm[j],7,'f',4);
			out << strong;
			if(m_XTr1[j]<990.0)
			{
				if(FileType==1) strong=QString("  %1  %2").arg(m_XTr1[j],6,'f',4).arg(m_XTr2[j],6,'f',4);
				else            strong=QString(",%1,%2").arg(m_XTr1[j],6,'f',4).arg(m_XTr2[j],6,'f',4);
				out << strong;
			}
			if(FileType==1) strong=QString("  %1  %2  %3\n").arg(m_Cpmn[j],7,'f',4).arg(m_HMom[j],7,'f',4).arg(m_XCp[j],7,'f',4);
			else            strong=QString(",%1,%2,%3\n").arg(m_Cpmn[j],7,'f',4).arg(m_HMom[j],7,'f',4).arg(m_XCp[j],7,'f',4);
			out << strong;
        }
	}
	out << "\n\n";
}


void CPolar::ResetPolar()
{
//	Removes all existing OpPoints results from polar
    m_Alpha.clear();
    m_Cl.clear();
    m_Cd.clear();
    m_Cdp.clear();
    m_Cm.clear();
    m_XTr1.clear();
    m_XTr2.clear();
    m_HMom.clear();
    m_Cpmn.clear();
    m_ClCd.clear();
    m_RtCl.clear();
    m_Cl32Cd.clear();
    m_Re.clear();
	m_XCp.clear();
}

CPolar *CPolar::newBySerialize() {
	CPolar* polar = new CPolar ();
	polar->serialize();
	return polar;
}

void CPolar::AddData(OpPoint *pOpPoint)
{
	//Adds the OpPoint data to the Polar object for further use
	if(!pOpPoint->m_bVisc) return;
	m_ACrit = pOpPoint->ACrit;
	AddPoint(pOpPoint->Alpha, pOpPoint->Cd, pOpPoint->Cdp, pOpPoint->Cl, pOpPoint->Cm,
			 pOpPoint->Xtr1, pOpPoint->Xtr2, pOpPoint->m_TEHMom, pOpPoint->Cpmn, pOpPoint->Reynolds,
			 pOpPoint->m_XCP);
}


void CPolar::AddData(void* ptrXFoil)
{
	XFoil *pXFoil = (XFoil*)ptrXFoil;

	if(!pXFoil->lvisc) return;
	double alpha = pXFoil->alfa*180.0/PI;
	m_ACrit = pXFoil->acrit;

	AddPoint(alpha, pXFoil->cd, pXFoil->cdp, pXFoil->cl, pXFoil->cm, pXFoil->xoctr[1],
			 pXFoil->xoctr[2], pXFoil->hmom, pXFoil->cpmn, pXFoil->reinf1, pXFoil->xcp);
	
}

void CPolar::AddPoint(double Alpha, double Cd, double Cdp, double Cl, double Cm, double Xtr1,
					  double Xtr2, double HMom, double Cpmn, double Reynolds, double XCp)
{
	int i;
	bool bInserted = false;
	int size = (int)m_Alpha.size();
	if(size)
	{
		for ( i=0; i<size; i++)
		{
			if(m_PolarType!=FIXEDAOAPOLAR)
			{
                if (fabs(Alpha - m_Alpha[i]) < 0.001)
				{
					// then erase former result
					m_Alpha[i] =  Alpha;
					m_Cd[i]    =  Cd;
					m_Cdp[i]   =  Cdp;
					m_Cl[i]    =  Cl;
					m_Cm[i]    =  Cm;
					m_XTr1[i]  =  Xtr1;
					m_XTr2[i]  =  Xtr2;
					m_HMom[i]  =  HMom;
					m_Cpmn[i]  =  Cpmn;
					m_ClCd[i]  =  Cl/Cd;
					m_XCp[i]   =  XCp;

					if(Cl>0.0)	 m_RtCl[i] = 1.0/sqrt(Cl);
					else         m_RtCl[i] = 0.0;
                    if (Cl>=0.0) m_Cl32Cd[i] = pow(Cl, 1.5)/ Cd;
                    else         m_Cl32Cd[i] = -pow(-Cl, 1.5)/ Cd;

					if(m_PolarType==FIXEDSPEEDPOLAR)	       m_Re[i] =  Reynolds;
					else if (m_PolarType==FIXEDLIFTPOLAR)
                    {
						if(Cl>0.0) m_Re[i] =  Reynolds/ sqrt(Cl);
						else m_Re[i] = 0.0;
					}
					else if (m_PolarType==RUBBERCHORDPOLAR)
                    {
						if(Cl>0.0) m_Re[i] =  Reynolds/(Cl);
						else m_Re[i] = 0.0;
					}


					bInserted = true;
					break;
				}
                else if (Alpha < m_Alpha[i])
                {
                    // sort by crescending alphas
                    m_Alpha.insert(i, Alpha);
                    m_Cd.insert(i, Cd);
                    m_Cdp.insert(i, Cdp);
                    m_Cl.insert(i, Cl);
                    m_Cm.insert(i, Cm);
                    m_XTr1.insert(i, Xtr1);
                    m_XTr2.insert(i, Xtr2);
                    m_HMom.insert(i, HMom);
                    m_Cpmn.insert(i, Cpmn);
                    m_ClCd.insert(i, Cl/Cd);
                    m_XCp.insert(i, XCp);

                    if(Cl>0.0)	 m_RtCl.insert(i, 1.0/sqrt(Cl));
                    else         m_RtCl.insert(i, 0.0);

                    if (Cl>=0.0) m_Cl32Cd.insert(i,pow(Cl, 1.5)/ Cd);
                    else         m_Cl32Cd.insert(i,-pow(-Cl, 1.5)/ Cd);

					if(m_PolarType==FIXEDSPEEDPOLAR)	 m_Re.insert(i, Reynolds);
					else if (m_PolarType==FIXEDLIFTPOLAR)
                    {
                        if(Cl>0) m_Re.insert(i, Reynolds/sqrt(Cl));
                        else m_Re[i] = 0.0;
                    }
					else if (m_PolarType==RUBBERCHORDPOLAR)
                    {
                        if(Cl>0.0) m_Re.insert(i, Reynolds/Cl);
                        else       m_Re.insert(i, 0.0);
                    }

                    bInserted = true;
                    break;
                }
            }
			else
			{
				//m_PolarType 4 polar, sort by Reynolds numbers
                if (fabs(Reynolds - m_Re[i]) < 0.1)
				{
					// then erase former result
					m_Alpha[i] =  Alpha;
					m_Cd[i]    =  Cd;
					m_Cdp[i]   =  Cdp;
					m_Cl[i]    =  Cl;
					m_Cm[i]    =  Cm;
					m_XTr1[i]  =  Xtr1;
					m_XTr2[i]  =  Xtr2;
					m_HMom[i]  =  HMom;
					m_Cpmn[i]  =  Cpmn;
					m_ClCd[i]  =  Cl/Cd;
					m_XCp[i]   =  XCp;

					if(Cl>0.0)	 m_RtCl[i] = 1.0/sqrt(Cl);
					else 	     m_RtCl[i] = 0.0;
                    if (Cl>=0.0) m_Cl32Cd[i] = pow(Cl, 1.5)/ Cd;
                    else         m_Cl32Cd[i] = -pow(-Cl, 1.5)/ Cd;
					m_Re[i] =  Reynolds;

					bInserted = true;
					break;
				}
                else if (Reynolds < m_Re[i]){// sort by crescending Reynolds numbers
                    m_Alpha.insert(i, Alpha);
                    m_Cd.insert(i, Cd);
                    m_Cdp.insert(i, Cdp);
                    m_Cl.insert(i, Cl);
                    m_Cm.insert(i, Cm);
                    m_XTr1.insert(i, Xtr1);
                    m_XTr2.insert(i, Xtr2);
                    m_HMom.insert(i, HMom);
                    m_Cpmn.insert(i, Cpmn);
                    m_ClCd.insert(i, Cl/Cd);
                    m_XCp.insert(i, XCp);

                    if(Cl>0.0)   m_RtCl.insert(i, 1.0/(double)sqrt(Cl));
                    else         m_RtCl.insert(i, 0.0);
                    if (Cl>=0.0) m_Cl32Cd.insert(i, pow(Cl, 1.5)/ Cd);
                    else         m_Cl32Cd.insert(i,-pow(-Cl, 1.5)/ Cd);

                    m_Re.insert(i, Reynolds);

                    bInserted = true;
                    break;
                }
            }
		}
	}
    if(!bInserted)
    {
        // data is appended at the end
        m_Alpha.insert(size, Alpha);
        m_Cd.insert(size, Cd);
        m_Cdp.insert(size, Cdp);
        m_Cl.insert(size, Cl);
        m_Cm.insert(size, Cm);
        m_XTr1.insert(size, Xtr1);
        m_XTr2.insert(size, Xtr2);
        m_HMom.insert(size, HMom);
        m_Cpmn.insert(size, Cpmn);
        m_ClCd.insert(size, Cl/Cd);
        m_XCp.insert(size, XCp);

        if(Cl>0.0)   m_RtCl.insert(size, 1.0/(double)sqrt(Cl));
        else         m_RtCl.insert(size, 0.0);
        if (Cl>=0.0) m_Cl32Cd.insert(size,(double)pow(Cl, 1.5)/ Cd);
        else         m_Cl32Cd.insert(size,-(double)pow(-Cl, 1.5)/ Cd);

		if(m_PolarType==FIXEDSPEEDPOLAR || m_PolarType==FIXEDAOAPOLAR) m_Re.insert(size, Reynolds);
		else if (m_PolarType==FIXEDLIFTPOLAR)
        {
            if(Cl>0) m_Re.insert(size, Reynolds/(double) sqrt(Cl));
            else     m_Re.insert(size, 0.0);
        }
		else if (m_PolarType==RUBBERCHORDPOLAR)
        {
            if(Cl>0.0) m_Re.insert(size, Reynolds/Cl);
            else     m_Re.insert(size, 0.0);
        }
    }
}


void CPolar::Copy(CPolar *pPolar)
{
	int i;
	int size  = (int)m_Alpha.size();
	for(i=size-1; i>=0; i--)
		Remove(i);

	size  = (int)pPolar->m_Alpha.size();
	for(i=0; i<size; i++)
	{
		m_Alpha.insert(i,  pPolar->m_Alpha[i]);
		m_Cd.insert(i,     pPolar->m_Cd[i]);
		m_Cdp.insert(i,    pPolar->m_Cdp[i]);
		m_Cl.insert(i,     pPolar-> m_Cl[i]);
		m_Cm.insert(i,     pPolar->m_Cm[i]);
		m_XTr1.insert(i,   pPolar->m_XTr1[i]);
		m_XTr2.insert(i,   pPolar->m_XTr2[i]);
		m_HMom.insert(i,   pPolar->m_HMom[i]);
		m_Cpmn.insert(i,   pPolar->m_Cpmn[i]);
		m_ClCd.insert(i,   pPolar->m_ClCd[i]);
		m_RtCl.insert(i,   pPolar->m_RtCl[i]);
		m_Cl32Cd.insert(i, pPolar->m_Cl32Cd[i]);
		m_Re.insert(i,     pPolar->m_Re[i]);
		m_XCp.insert(i,    pPolar->m_XCp[i]);
	}
}

void CPolar::ExportPolarNREL(QTextStream &out) {
    /* export the polar in NREL format compatibel with AeroDyn 13 and WT_perf 3 */
    QDate date = QDate::currentDate();
    QTime time = QTime::currentTime();
    double clMaxAngle, cdMin, cdMinAngle;
    double cnSlopeZeroLift, cnStallPositiveAngle, cnStallNegativeAngle;
    double dummy; // used for not needed values
    getClMaximum (dummy, clMaxAngle);
    getCdMinimum (cdMin, cdMinAngle);
    GetLinearizedCn (dummy, cnSlopeZeroLift);
    GetCnAtStallAngles (cnStallPositiveAngle, cnStallNegativeAngle);

    /* Version that really is compatibel with AeroDyn 13 */
    out <<  "AeroDyn airfoil file Created with " << g_mainFrame->m_VersionName << " on "<<date.toString("dd.MM.yyyy") << " at " << time.toString("hh:mm:ss") << " Compatible with AeroDyn v13.0." << endl <<
            "Polar \"" << getName() << "\" on Foil \"" << getParent()->getName() << "\" :: generated by QBlade"<< endl <<
            QString("%1").arg(1, -15) <<
            "Number of airfoil tables in this file" << endl <<
            QString("%1").arg(0, -15) <<
            "Table ID parameter" << endl <<
            QString("%1").arg(clMaxAngle, -15, 'f', 2) <<
            "Stall angle (deg)" << endl <<
            "0              No longer used, enter zero" << endl <<
            "0              No longer used, enter zero" << endl <<
            "0              No longer used, enter zero" << endl <<
            QString("%1").arg(GetZeroLiftAngle(), -15, 'f', 2) <<
            "Angle of attack for zero Cn for linear Cn curve (deg)" << endl <<
            QString("%1").arg(cnSlopeZeroLift, -15, 'f', 5) <<
            "Cn slope for zero lift for linear Cn curve (1/rad)" << endl <<
            QString("%1").arg(cnStallPositiveAngle, -15, 'f', 4) <<
            "Cn at stall value for positive angle of attack for linear Cn curve" << endl <<
            QString("%1").arg(cnStallNegativeAngle, -15, 'f', 4) <<
            "Cn at stall value for negative angle of attack for linear Cn curve" << endl <<
            QString("%1").arg(cdMinAngle, -15, 'f', 2) <<
            "Angle of attack for minimum CD (deg)" << endl <<
            QString("%1").arg(cdMin, -14, 'f', 4) <<
            " Minimum CD value" << endl;

    for (int i = 0; i < m_Alpha.size(); ++i) {
        out << QString("%1").arg(m_Alpha[i], 10, 'f', 2) <<
               QString("%1").arg(m_Cl[i], 10, 'f', 4) <<
               QString("%1").arg(m_Cd[i], 10, 'f', 4) <</*
               QString("%1").arg(m_Cm[i], 10, 'f', 4) << */endl;  // m_Cm is not filled with values
    }

    if (!m_Alpha.size()){
    out << QString(tr("- exported polar did not contain any data -"));
    }
}

void CPolar::GetCnAtStallAngles(double &cnPosStallAlpha, double &cnNegStallAlpha)
{
    //the stall angles are best seen in the peaks of the Cl/Cd over Alpha curve.
    //Only the part between -50 and +50 deg is considered
    QList <double> ClCd;
    QList <double> Alpha;
    QList <double> Cn;
    cnNegStallAlpha = 0;
    cnPosStallAlpha = 0;

    if (m_Cl.size() != m_Cd.size()) return;

    for (int i = 0; i < m_Cl.size(); ++i)
    {
        if ((m_Alpha[i] > -50) && (m_Alpha[i] < 50) && (m_Cd[i] != 0.0) )
        {
            Alpha.append(m_Alpha[i]);
            ClCd.append(m_Cl[i] / m_Cd[i]);
            Cn.append( m_Cl[i]*cos(m_Alpha[i]*PI/180) + m_Cd[i]*sin(m_Alpha[i]*PI/180));
        }

    }

    bool bNegStallFound = 0;
    bool bPosStallFound = 0;

    //get 2 Inflection points
    for (int i = 0; i < ClCd.size()-1; ++i)
    {
       if (!bNegStallFound)
        if (ClCd[i+1] > ClCd[i])
         {
             //double negStallAlpha = Alpha[i];
             bNegStallFound = 1;
             cnNegStallAlpha = Cn[i];
         }


       if (!bPosStallFound && bNegStallFound)
        if (ClCd[i+1] < ClCd[i])
         {
             //double posStallAlpha = Alpha[i];
             bPosStallFound = 1;
             cnPosStallAlpha = Cn[i];
         }

    }
}


void CPolar::GetLinearizedCn(double &Alpha0, double &slope)
{
    //linearize Cn vs. Alpha set of points by least square method
    QList <double> Cn;
    QList <double> Alpha;

    int n;
    double alpha0L = GetZeroLiftAngle();

    //calculate Cn
    //all points in the range -3/+5 deg around zero lift angle is taken into account
    for (int i = 0; i < m_Cl.size(); ++i)
    {
        if ((m_Alpha[i] > (alpha0L - 3)) && (m_Alpha[i] < (alpha0L + 3)))
        {
            Cn.append( m_Cl[i]*cos(m_Alpha[i]*PI/180) + m_Cd[i]*sin(m_Alpha[i]*PI/180));
            Alpha.append(m_Alpha[i]);
        }

    }

    n = (int)Cn.size();

    if(n<=1)
    {
        Alpha0 = 0.0;
        slope = 2.0*PI*PI/180.0;
        return;
    }


    double fn = (double)n;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double b1, b2;

    for (int k=0; k<n; k++)
    {
        sum1 += Cn[k] * Alpha[k];
        sum2 += Alpha[k];
        sum3 += Cn[k];
        sum4 += Alpha[k] * Alpha[k];
    }


    if(fn*sum4 == sum2*sum2 || fn*sum1 == sum2 * sum3) {//very improbable...
        Alpha0 = 0.0;
        slope = 2.0*PI*PI/180.0;
        return;
    }


    b1 = (fn*sum1 - sum2 * sum3)/(fn*sum4 - sum2*sum2);
    b2 = (sum3 - b1 * sum2)/fn;

    slope  = b1*180/PI; //in cn/alpha[rad]
    Alpha0 = -b2/b1;
}

void CPolar::getCdMinimum(double &cdMin, double &cdMinAngle) {
    if (m_Cd.empty()) {
        cdMin = 0.0;
        cdMinAngle = 0.0;
    } else {
        int minIndex = m_Cd.size() / 2;  // guess some value in the middle
        for (int i = 0; i < m_Cd.size(); ++i) {
            // search Cd minimum between -20 and +20 degree
            if (m_Cd[i] < m_Cd[minIndex] && m_Alpha[i] > -20 && m_Alpha[i] < +20) {
                minIndex = i;
            }
        }
        cdMin = m_Cd[minIndex];
        cdMinAngle = m_Alpha[minIndex];
    }
}

void CPolar::getClMaximum(double &clMax, double &clMaxAngle) {
    if (m_Cl.empty()) {
        clMax = 0.0;
        clMaxAngle = 0.0;
    } else {
        int maxIndex = 0;
        for (int i = 0; i < m_Cl.size(); ++i) {
            if (m_Cl[i] > m_Cl[maxIndex]) {
                maxIndex = i;
            }
        }
        clMax = m_Cl[maxIndex];
        clMaxAngle = m_Alpha[maxIndex];
    }
}

QStringList CPolar::prepareMissingObjectMessage() {
	if (g_polarStore.isEmpty()) {
		QStringList message = CFoil::prepareMissingObjectMessage();
		if (message.isEmpty()) {
			message = QStringList(">>> Create a new Polar in the XFOIL Direct Analysis Module");
		}
		message.prepend("- No Polar in Database");
		return message;
	} else {
		return QStringList();
	}
}

void CPolar::Remove(int i)
{
	m_Alpha.removeAt(i);
	m_Cl.removeAt(i);
	m_Cd.removeAt(i);
	m_Cdp.removeAt(i);
	m_Cm.removeAt(i);
	m_XTr1.removeAt(i);
	m_XTr2.removeAt(i);
	m_HMom.removeAt(i);
	m_Cpmn.removeAt(i);
	m_ClCd.removeAt(i);
	m_RtCl.removeAt(i);
	m_Cl32Cd.removeAt(i);
	m_Re.removeAt(i);
	m_XCp.removeAt(i);
}

void CPolar::GetAlphaLimits(double &amin, double &amax)
{
    if(!m_Alpha.size()){
		amin = 0.0;
		amax = 0.0;
	}
    else
    {
        amin = m_Alpha[0];
        amax = m_Alpha[m_Alpha.size()-1];
	}
}

void CPolar::getClLimits(double &Clmin, double &Clmax)
{
    if(!m_Cl.size())
    {
		Clmin = 0.0;
		Clmax = 0.0;
	}
    else
    {
		Clmin = 10000.0;
		Clmax =-10000.0;
		double Cl;
        for (int i=0;i<m_Cl.size(); i++)
        {
            Cl = m_Cl[i];
			if(Clmin>Cl) Clmin = Cl;
			if(Clmax<Cl) Clmax = Cl;
		}
	}
}

double CPolar::GetZeroLiftAngle()
{
	/////////////// new code NM ///////////////
	// only consider reasonable aoa range
	QList<double> clRange;
	QList<double> alphaRange;
	
	for (int i = 0; i < m_Alpha.size(); ++i)
	{
	    if (m_Alpha[i] > -10 && m_Alpha[i] < 10)
	    {
	        clRange.append(m_Cl[i]);
	        alphaRange.append(m_Alpha[i]);
	    }
	
	}
	/////////// end new code NM ///////////////
	
	
	double Clmin =  1000.0;
	double Clmax = -1000.0;
    for (int i=0; i<clRange.size(); i++)
    {
		Clmin = qMin(Clmin, clRange[i]);
		Clmax = qMax(Clmax, clRange[i]);
	}
	if(!(Clmin<0.0) || !(Clmax>0.0))
		return 0.0;
	
	int k=0;
//	double rr  = clRange[k];
//	double rr1 = clRange[k+1];

    while (clRange[k+1]<0.0)
    {
//		rr  = clRange[k];
//		rr1 = clRange[k+1];
		k++;
	}
    if(k+1>=alphaRange.size()) return 0.0;
	double Alpha0 = alphaRange[k] + (alphaRange[k+1]-alphaRange[k])*(0.0-clRange[k])/(clRange[k+1]-clRange[k]);
	return Alpha0;

}

void CPolar::GetLinearizedCl(double &Alpha0, double &slope)
{
	// linearize Cl vs. Alpha set of points by least square method
    int n = (int)m_Cl.size();

    if(n<=1)
    {
		Alpha0 = 0.0;
		slope = 2.0*PI*PI/180.0;
		return;
	}

	double fn = (double)n;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	double sum4 = 0.0;
	double b1, b2;

	for (int k=0; k<n; k++){
		sum1 += m_Cl[k] * m_Alpha[k];
		sum2 += m_Alpha[k];
		sum3 += m_Cl[k];
		sum4 += m_Alpha[k] * m_Alpha[k];
	}
	if(fn*sum4 == sum2*sum2 || fn*sum1 == sum2 * sum3) {//very improbable...
		Alpha0 = 0.0;
		slope = 2.0*PI*PI/180.0;
		return;
	}

	b1 = (fn*sum1 - sum2 * sum3)/(fn*sum4 - sum2*sum2);
	b2 = (sum3 - b1 * sum2)/fn;

	slope  = b1; //in cl/�
	Alpha0 = -b2/b1;

}

void CPolar::GetPolarProperties(QString &PolarProperties, bool bData)
{
	QString strong;
	PolarProperties = getName() +"\n\n";

//	PolarProperties += QObject::tr("Parent foil")+" = "+ m_FoilName+"\n";

//	strong = QString(QObject::tr("Analysis Type")+" = %1\n").arg(m_PolarType);
	PolarProperties.clear();

	strong = QString(QObject::tr("Type")+" = %1").arg(m_PolarType);
	if(m_PolarType==FIXEDSPEEDPOLAR)      strong += " ("+QObject::tr("Fixed speed") +")\n";
	else if(m_PolarType==FIXEDLIFTPOLAR) strong += " ("+QObject::tr("Fixed lift") +")\n";
	else if(m_PolarType==FIXEDAOAPOLAR) strong += " ("+QObject::tr("Fixed angle of attack") +")\n";
	PolarProperties += strong;

	if(m_PolarType==FIXEDSPEEDPOLAR)
	{
		strong = QString(QObject::tr("Reynolds number")+" = %1\n").arg(m_Reynolds,0,'f',0);
		PolarProperties += strong;
		strong = QString(QObject::tr("Mach number") + " = %1\n").arg(m_Mach,5,'f',2);
		PolarProperties += strong;
	}
	else if(m_PolarType==FIXEDLIFTPOLAR)
	{
		strong = QString("Re.sqrt(Cl) = %1\n").arg(m_Reynolds,0,'f',0);
		PolarProperties += strong;
		strong = QString("Ma.sqrt(Cl) = %1\n").arg(m_Mach,5,'f',2);
		PolarProperties += strong;
	}
	else if(m_PolarType==RUBBERCHORDPOLAR)
	{
		strong = QString(QObject::tr("Re.Cl")+" = %1\n").arg(m_Reynolds,0,'f',0);
		PolarProperties += strong;
		strong = QString(QObject::tr("Mach number") + " = %1\n").arg(m_Mach,5,'f',2);
		PolarProperties += strong;
	}
	else if(m_PolarType==FIXEDAOAPOLAR)
	{
		strong = QString(QObject::tr("Alpha")+" = %1"+QString::fromUtf8("deg")+"\n").arg(m_ASpec,7,'f',2);
		PolarProperties += strong;
		strong = QString(QObject::tr("Mach number") + " = %1\n").arg(m_Mach,5,'f',2);
		PolarProperties += strong;
	}


	strong = QString(QObject::tr("NCrit") + " = %1\n").arg(m_ACrit,6,'f',2);
	PolarProperties += strong;

	strong = QString(QObject::tr("Forced top trans.") + " = %1\n").arg(m_XTop,6,'f',2);
	PolarProperties += strong;

	strong = QString(QObject::tr("Forced bottom trans.") + " = %1\n").arg(m_XBot,6,'f',2);
	PolarProperties += strong;

	strong = QString(QObject::tr("Number of data points") +" = %1").arg(m_Alpha	.size());
	PolarProperties += "\n" +strong;

	if(!bData) return;
	QTextStream out;
	strong.clear();
	out.setString(&strong);
	ExportPolar(out, g_mainFrame->m_ExportFileType, true);
	PolarProperties += "\n"+strong;
}
