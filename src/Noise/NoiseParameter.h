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
    bool Lowson;
    bool Amiet;
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
    double sects;
    double rot_speed;
    double u_wind_speed;
    double TSRtd;
    double dstar_type;
    double phi_type;
    double rot_speed_calc;
    double u_wind_speed_calc;
    double TSR_calc;
    double obs_x_pos;
    double obs_y_pos;
    double obs_z_pos;
    bool rot_speed_check;
    bool u_wind_speed_check;
    bool TSR_check;

    //Turbulence parameters - Alexandre MOD
    double TurbulenceIntensity;
    double IntegralLengthScale;
};

#endif // NOISEPARAMETER_H
