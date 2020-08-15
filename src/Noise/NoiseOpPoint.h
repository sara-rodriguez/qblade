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
    NoiseOpPoint(double reynolds, double mach, double alpha);
	
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
    double m_reynolds, m_mach, m_alpha, x, m_rot_speed, m_u_wind_speed, m_TSRtd, m_obs_x_pos, m_obs_y_pos, m_obs_z_pos,  m_obs_x_pos_rotor, m_obs_y_pos_rotor, m_obs_z_pos_rotor, m_tower_to_hub_distance, m_initial_azimuth, m_time, m_shear_roughness, m_shear_height, m_shear_speed, m_yaw_angle, m_valRel_TE, m_valReu_TE, m_valMal_TE, m_valMau_TE, m_valAOAl_TE, m_valAOAu_TE, m_valRel_LE, m_valReu_LE, m_valMal_LE, m_valMau_LE; //Sara
	OpPoint *m_opPoint;
//Sara
bool m_rot_speed_check, m_u_wind_speed_check, m_TSR_check, m_shear_check, m_valRel_TE_check, m_valReu_TE_check, m_valMal_TE_check, m_valMau_TE_check, m_valAOAl_TE_check, m_valAOAu_TE_check, m_valRel_LE_check, m_valReu_LE_check, m_valMal_LE_check, m_valMau_LE_check;//Sara
int m_Lowson_type, m_number_loops, m_timesteps, m_qs3DSim;
//Sara
};

#endif // NOISEOPPOINT_H
