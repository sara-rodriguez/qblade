/****************************************************************************

    Quaternion Class
	Copyright (C) 2008-2009 Andre Deperrois adeperrois@xflr5.com

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

#ifndef ARCBALL_H
#define ARCBALL_H

#include "CVector.h"
#include "../Params.h"
class Quaternion
{
private:
	double theta;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t15, t19, t20, t24;
	CVector R;	

public:
	double a, qx, qy,qz;
	
	void QuattoMat(double m[][4]);
	void Normalize();


	void operator*=(Quaternion Q);
	void operator ~();
	void operator =(Quaternion Q);
	Quaternion operator *(Quaternion Q);

	//inline constructors
	Quaternion(void)
	{
		a=0.0; qx= 0.0; qy=0.0; qz = 0.0;
		theta = 0.0;
		Settxx();
	};

	Quaternion(double const &t, double const &x, double const &y, double const &z)
	{
		a=t; qx= x; qy=y; qz = z;
		theta = 2.0*acos(t);
		Settxx();
	};

	Quaternion(double const &Angle, CVector const &R)
	{	
		CVector N;
		N = R;
		N.Normalize();
		theta = Angle*PI/180.0;

		a = cos(theta/2.0);
		double sina = sin(theta/2.0);

		qx = N.x*sina;
		qy = N.y*sina;
		qz = N.z*sina;
		Settxx();
	};

	//inline functions
	void Conjugate(CVector const &Vin, CVector &Vout)
	{
		Vout.x = 2.0*( (t8 + t10)*Vin.x + (t6 -  t4)*Vin.y + (t3 + t7)*Vin.z ) + Vin.x;
		Vout.y = 2.0*( (t4 +  t6)*Vin.x + (t5 + t10)*Vin.y + (t9 - t2)*Vin.z ) + Vin.y;
		Vout.z = 2.0*( (t7 -  t3)*Vin.x + (t2 +  t9)*Vin.y + (t5 + t8)*Vin.z ) + Vin.z;
	};
	void Conjugate(CVector &V)
	{
		R.x = V.x;
		R.y = V.y;
		R.z = V.z;
	
		V.x = 2.0*( (t8 + t10)*R.x + (t6 -  t4)*R.y + (t3 + t7)*R.z ) + R.x;
		V.y = 2.0*( (t4 +  t6)*R.x + (t5 + t10)*R.y + (t9 - t2)*R.z ) + R.y;
		V.z = 2.0*( (t7 -  t3)*R.x + (t2 +  t9)*R.y + (t5 + t8)*R.z ) + R.z;
	};
	void Conjugate(double &x, double &y, double &z)
	{
		R.x = x;
		R.y = y;
		R.z = z;
	
		x = 2.0*( (t8 + t10)*R.x + (t6 -  t4)*R.y + (t3 + t7)*R.z ) + R.x;
		y = 2.0*( (t4 +  t6)*R.x + (t5 + t10)*R.y + (t9 - t2)*R.z ) + R.y;
		z = 2.0*( (t7 -  t3)*R.x + (t2 +  t9)*R.y + (t5 + t8)*R.z ) + R.z;
	};
	
	void Settxx()
	{
		t2 =   a*qx;
		t3 =   a*qy;
		t4 =   a*qz;
		t5 =  -qx*qx;
		t6 =   qx*qy;
		t7 =   qx*qz;
		t8 =  -qy*qy;
		t9 =   qy*qz;
		t10 = -qz*qz;	
	};


	void Set(double const &real, double const &x, double const &y, double const &z)
	{	
		a = real;
		qx = x;
		qy = y;
		qz = z;
		Settxx();
	};

	void Set(double const &Angle, CVector const &R)
	{	
		CVector N;
		N = R;
		N.Normalize();
		theta = Angle*PI/180.0;

		a = cos(theta/2.0);
		double sina = sin(theta/2.0);

		qx = N.x*sina;
		qy = N.y*sina;
		qz = N.z*sina;
		Settxx();
	};
};
#endif
