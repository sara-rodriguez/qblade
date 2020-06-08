#ifndef NOISEPARAMETER_H
#define NOISEPARAMETER_H

#include <QList>
#include <QVector>

class OpPoint;
class NoiseOpPoint;


/* This class holds all parameters for a 2D noise simulation. Its only purpose is to reduce the amount of variables
 * in the NoiseSimulation class and keep its namespace clean.
 * */

class NoiseParameter
{
public:
	enum OpPointSource {OnePolar, MultiplePolars, OriginalBpm};
	enum TransitionType {FullyTurbulent, TransitionFlow};

    NoiseParameter();
	void serialize ();
	void restorePointers ();

	QList<NoiseOpPoint*> prepareNoiseOpPointList();  // transfers ownership to caller, don't forget to delete the list

	OpPointSource opPointSource;  // where the analyzed OpPoints originate from
	QVector<OpPoint*> analyzedOpPoints;  // list of the analyzed oppoints
		
	double wettedLength;
	double distanceObsever;
	double directivityGreek;
	double directivityPhi;
	bool highFreq;  // NM the two frequency members are not included in the GUI
	bool lowFreq;
	bool suctionSide; //SPLs
	bool pressureSide; //SPLp
	bool separatedFlow; //SPLa

    //Sara
    //LE
    int Lowson_type;
    //Sara
	
	//XFoil correlation
	double dStarChordStation;
	double dStarScalingFactor;
	double eddyConvectionMach;
	double originalMach;
	double originalChordLength;
	double originalVelocity;
	
	//Original BPM correlations
	double aoa;
	double chordBasedReynolds;
	TransitionType transition;
    
    //3D Sara
    bool rot_speed_check;
    bool u_wind_speed_check;
    bool shear_check;
    bool TSR_check;

    int qs3DSim;
    int phi_type;
    int theta_type;
    int dstar_type;
    int state_ss_us;
    int number_loops;
    int timesteps;
    int rotation_type;
    int anglesteps;

    double obs_x_pos;
    double obs_y_pos;
    double obs_z_pos;
    double obs_x_pos_rotor;
    double obs_y_pos_rotor;
    double obs_z_pos_rotor;
    double tower_height;
    double rot_speed;
    double u_wind_speed;
    double TSRtd;
    double tower_to_hub_distance;
    double initial_azimuth;
    double shear_roughness;
    double shear_height;
    double shear_speed;
    double time;
    double yaw_angle;
    double D_starred_index_user[40];
    double D_starred_S_user[40];
    double D_starred_P_user[40];
        //Sara

    //Turbulence parameters - Alexandre MOD
    double TurbulenceIntensity;
    double IntegralLengthScale;
};

#endif // NOISEPARAMETER_H
