#ifndef NOISEOPPOINT_H
#define NOISEOPPOINT_H
#include <QVector>//Sara

class OpPoint;


/* This class is a conveniance wrapper for OpPoint. Its purpose is to provide "fake" OpPoints with values only for
 * alpha and reynolds and to reduce the interface to some getter functions. If the m_opPoint pointer is set, this class
 * behaves like a wrapper. If it is not set, this class returns its own values for alpha and reynolds.
 * */

class NoiseOpPoint
{
public:
	NoiseOpPoint(OpPoint *opPoint);
    NoiseOpPoint(double reynolds, double alpha);
	
    double getMach();//Sara
	double getReynolds();
	double getAlphaDegree();
	double getAlphaDegreeAbsolute();
	
	// the following are only available for true OpPoints and crash for fake OpPoints
	int getNSide1();
	int getNSide2();    
	double getXValue(int index, int topOrBot);
	double getDstrAt(int x, int y);

    //Sara
    double getAlphaAt(int x, int y);
    double getReynoldsAt(int x, int y);
    double getMachAt(int x, int y);
    //Sara
	
private:
    double m_reynolds, m_alpha, m_sects, x, m_rot_speed, m_u_wind_speed, m_TSRtd, m_dstar_user, m_rot_speed_calc, m_u_wind_speed_calc, m_TSR_calc, m_obs_x_pos, m_obs_y_pos, m_obs_z_pos, m_mach, m_Tt, m_Ts, m_As; //Sara
	OpPoint *m_opPoint;
//Sara
bool m_rot_speed_check, m_u_wind_speed_check, m_TSR_check;
int m_Lowson_type;
//Sara
};

#endif // NOISEOPPOINT_H
