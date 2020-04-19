#ifndef NOISECALCULATION_H
#define NOISECALCULATION_H

#include <QVector>
#include <QtMath> //Sara

class NoiseParameter;
class NoiseOpPoint;


/* This class performes the calculations required by the noise simulation. Its purpose is to reduce the amount of
 * variables in the NoiseSimulation class and keep its namespace clean.
 * WARNING: When serializing this class only the publicly accsessible members are taken into account to avoid the
 * overhead of storing intermediate results. Therefore, after loading a calculation the calculate() method must be
 * run to restore the intermediate results.
 * */

class NoiseCalculation
{
public:
	enum AirfoilSide {PressureSide, SuctionSide};
	
	typedef QVector< QVector<double> > TwoDVector;
	
	static constexpr double SWITCHING_ANGLE2 = 12.5;
    static constexpr int FREQUENCY_TABLE_SIZE = 34; //Alexandre MOD
	static const double AWeighting[FREQUENCY_TABLE_SIZE];  // 1/3 octave band frequency
	static const double BWeighting[FREQUENCY_TABLE_SIZE];
	static const double CWeighting[FREQUENCY_TABLE_SIZE];
    static const QVector<double> CENTRAL_BAND_FREQUENCY;
	
	NoiseCalculation();
	void serialize ();
	
	void setNoiseParam (NoiseParameter *parameter) { m_parameter = parameter; }
	void calculate();  // can throw NoiseException

    void verifydeltafor3d();//Sara
    void calculateqs3d_graphics();//Sara
    void calculateqs3d_blade();//Sara
    void calculateqs3d_rotor();//Sara
    void unsteady_angles_calc();//Sara
    void positive_graphs();//Sara urgente

    //Sara
    //for quasi 3d rotor
    TwoDVector unsteady_angles() const { return m_unsteady_angles; }
    //Sara
	
	// NM the arrays containing the graph data
	TwoDVector SPLadB() const { return m_SPLadB; }
	TwoDVector SPLsdB() const { return m_SPLsdB; }
	TwoDVector SPLpdB() const { return m_SPLpdB; }
	TwoDVector SPLdB() const { return m_SPLdB; }
	TwoDVector SPLdBAW() const { return m_SPLdBAW; }
	TwoDVector SPLdBBW() const { return m_SPLdBBW; }
	TwoDVector SPLdBCW() const { return m_SPLdBCW; }
//Alexandre MOD
    TwoDVector SPL_LEdB() const { return m_SPL_LEdB; }
    TwoDVector SPL_LEdBAW() const { return m_SPL_LEdBAW; }
    TwoDVector SPL_LEdBBW() const { return m_SPL_LEdBBW; }
    TwoDVector SPL_LEdBCW() const { return m_SPL_LEdBCW; }
    //Alexandre MOD

//Sara
    TwoDVector SPLadB3d() const { return m_SPLadB3d; }
    TwoDVector SPLsdB3d() const { return m_SPLsdB3d; }
    TwoDVector SPLpdB3d() const { return m_SPLpdB3d; }
    TwoDVector SPLdB3d() const { return m_SPLdB3d; }
    TwoDVector SPLdBAW3d() const { return m_SPLdBAW3d; }
    TwoDVector SPLdBBW3d() const { return m_SPLdBBW3d; }
    TwoDVector SPLdBCW3d() const { return m_SPLdBCW3d; }
    TwoDVector SPL_LEdB3d() const { return m_SPL_LEdB3d; }
    TwoDVector SPL_LEdBAW3d() const { return m_SPL_LEdBAW3d; }
    TwoDVector SPL_LEdBBW3d() const { return m_SPL_LEdBBW3d; }
    TwoDVector SPL_LEdBCW3d() const { return m_SPL_LEdBCW3d; }
    TwoDVector SPLadB3d_final() const { return m_SPLadB3d_final; }
    TwoDVector SPLsdB3d_final() const { return m_SPLsdB3d_final; }
    TwoDVector SPLpdB3d_final() const { return m_SPLpdB3d_final; }
    TwoDVector SPLdB3d_final() const { return m_SPLdB3d_final; }
    TwoDVector SPLdBAW3d_final() const { return m_SPLdBAW3d_final; }
    TwoDVector SPLdBBW3d_final() const { return m_SPLdBBW3d_final; }
    TwoDVector SPLdBCW3d_final() const { return m_SPLdBCW3d_final; }
    TwoDVector SPL_LEdB3d_final() const { return m_SPL_LEdB3d_final; }
    TwoDVector SPL_LEdBAW3d_final() const { return m_SPL_LEdBAW3d_final; }
    TwoDVector SPL_LEdBBW3d_final() const { return m_SPL_LEdBBW3d_final; }
    TwoDVector SPL_LEdBCW3d_final() const { return m_SPL_LEdBCW3d_final; }

    TwoDVector SPLadB3d_rotor() const { return m_SPLadB3d_rotor; }
    TwoDVector SPLsdB3d_rotor() const { return m_SPLsdB3d_rotor; }
    TwoDVector SPLpdB3d_rotor() const { return m_SPLpdB3d_rotor; }
    TwoDVector SPLdB3d_rotor() const { return m_SPLdB3d_rotor; }
    TwoDVector SPLdBAW3d_rotor() const { return m_SPLdBAW3d_rotor; }
    TwoDVector SPLdBBW3d_rotor() const { return m_SPLdBBW3d_rotor; }
    TwoDVector SPLdBCW3d_rotor() const { return m_SPLdBCW3d_rotor; }
    TwoDVector SPL_LEdB3d_rotor() const { return m_SPL_LEdB3d_rotor; }
    TwoDVector SPL_LEdBAW3d_rotor() const { return m_SPL_LEdBAW3d_rotor; }
    TwoDVector SPL_LEdBBW3d_rotor() const { return m_SPL_LEdBBW3d_rotor; }
    TwoDVector SPL_LEdBCW3d_rotor() const { return m_SPL_LEdBCW3d_rotor; }
    TwoDVector SPLadB3d_final_rotor() const { return m_SPLadB3d_final_rotor; }
    TwoDVector SPLsdB3d_final_rotor() const { return m_SPLsdB3d_final_rotor; }
    TwoDVector SPLpdB3d_final_rotor() const { return m_SPLpdB3d_final_rotor; }
    TwoDVector SPLdB3d_final_rotor() const { return m_SPLdB3d_final_rotor; }
    TwoDVector SPLdBAW3d_final_rotor() const { return m_SPLdBAW3d_final_rotor; }
    TwoDVector SPLdBBW3d_final_rotor() const { return m_SPLdBBW3d_final_rotor; }
    TwoDVector SPLdBCW3d_final_rotor() const { return m_SPLdBCW3d_final_rotor; }
    TwoDVector SPL_LEdB3d_final_rotor() const { return m_SPL_LEdB3d_final_rotor; }
    TwoDVector SPL_LEdBAW3d_final_rotor() const { return m_SPL_LEdBAW3d_final_rotor; }
    TwoDVector SPL_LEdBBW3d_final_rotor() const { return m_SPL_LEdBBW3d_final_rotor; }
    TwoDVector SPL_LEdBCW3d_final_rotor() const { return m_SPL_LEdBCW3d_final_rotor; }

    double Final_qs3d_alpha;
    double Final_qs3d_S;
    double Final_qs3d_P;
    double Final_qs3d_LE;
    double Final_qs3d;
    double Final_qs3d_alpha_rotor;
    double Final_qs3d_S_rotor;
    double Final_qs3d_P_rotor;
    double Final_qs3d_LE_rotor;
    double Final_qs3d_rotor;

//Sara
	
	// NM apparently needed for export as .txt only
	QVector<double> OASPL() const { return m_OASPL; }
	QVector<double> OASPLA() const { return m_OASPLA; }
	QVector<double> OASPLB() const { return m_OASPLB; }
	QVector<double> OASPLC() const { return m_OASPLC; }
	QVector<double> SPLALOG() const { return m_SPLALOG; }
	QVector<double> SPLSLOG() const { return m_SPLSLOG; }
	QVector<double> SPLPLOG() const { return m_SPLPLOG; }
//Alexandre MOD
        QVector<double> SPLLEdBAW() const { return m_SPLLEdBAW; }
        QVector<double> SPLLEdBBW() const { return m_SPLLEdBBW; }
        QVector<double> SPLLEdBCW() const { return m_SPLLEdBCW; }
        QVector<double> SPLlogLE() const { return m_SPLlogLE; }
//Alexandre MOD

        //Sara
        QVector<double> OASPL3d() const { return m_OASPL3d; }
        QVector<double> OASPLA3d() const { return m_OASPLA3d; }
        QVector<double> OASPLB3d() const { return m_OASPLB3d; }
        QVector<double> OASPLC3d() const { return m_OASPLC3d; }
        QVector<double> SPLALOG3d() const { return m_SPLALOG3d; }
        QVector<double> SPLSLOG3d() const { return m_SPLSLOG3d; }
        QVector<double> SPLPLOG3d() const { return m_SPLPLOG3d; }
        QVector<double> SPLLEdBAW3d() const { return m_SPLLEdBAW3d;}
        QVector<double> SPLLEdBBW3d() const { return m_SPLLEdBBW3d;}
        QVector<double> SPLLEdBCW3d() const { return m_SPLLEdBCW3d;}
        QVector<double> SPLlogLE3d() const { return m_SPLlogLE3d; }

        QVector<double> OASPL3d_rotor() const { return m_OASPL3d_rotor; }
        QVector<double> OASPLA3d_rotor() const { return m_OASPLA3d_rotor; }
        QVector<double> OASPLB3d_rotor() const { return m_OASPLB3d_rotor; }
        QVector<double> OASPLC3d_rotor() const { return m_OASPLC3d_rotor; }
        QVector<double> SPLALOG3d_rotor() const { return m_SPLALOG3d_rotor; }
        QVector<double> SPLSLOG3d_rotor() const { return m_SPLSLOG3d_rotor; }
        QVector<double> SPLPLOG3d_rotor() const { return m_SPLPLOG3d_rotor; }
        QVector<double> SPLLEdBAW3d_rotor() const { return m_SPLLEdBAW3d;}
        QVector<double> SPLLEdBBW3d_rotor() const { return m_SPLLEdBBW3d;}
        QVector<double> SPLLEdBCW3d_rotor() const { return m_SPLLEdBCW3d;}
        QVector<double> SPLlogLE3d_rotor() const { return m_SPLlogLE3d_rotor; }

        QVector<double> m_DStarInterpolatedS3d;
        QVector<double> m_DStarInterpolatedP3d;
        QVector<double> m_AlphaInterpolated3d;
        QVector<double> m_ReynoldsInterpolated3d;
        QVector<double> m_MachInterpolated3d;
        double c_const;
        double d_const;
        int m_Lowson_type;
    //Sara
	
private:
    void setupVectors();
    void setupVectorsqs3d();//Sara

	// calculation sub-functions
	double getK1(NoiseOpPoint* nop);
    double getDStarInterpolated(bool top, NoiseOpPoint *nop);  // can throw NoiseException
    double getDStarInterpolated3d(bool top, double chord,NoiseOpPoint *nop);  // Sara
    double getDH();
    double getDL();
    double getSt1();
    double getSt2(NoiseOpPoint *nop);
    double getBPMThickness(NoiseOpPoint *nop, AirfoilSide as);
    void preCalcA1(NoiseOpPoint* nop);
    void preCalcSPLa(NoiseOpPoint* nop);
    void preCalcSPLs(NoiseOpPoint* nop);
    void preCalcSPLp(NoiseOpPoint* nop);
    void calcSPLa(double alpha,int posOpPoint,int posFreq);
    void calcSPLs(int posOpPoint,int posFreq);
    void calcSPLp(int posOpPoint,int posFreq);
    void LECalc(int posOpPoint, int posFreq); //Alexandre MOD
    //Sara
    double calcXRS(double a, double XB, double YB);
    double calcYRS(double a, double XB, double YB);
    double calcZRS(double ZB, double r_0, double r_1);
    double calcInt_a(double YRS, double c_0, double c_1);
    double calcXRT(double XRS);
    double calcYRT(double b, double calc_int_a, double ZRS);
    double calcZRT(double b, double calc_int_a, double ZRS);
    double calcR_e(double XRT, double YRT, double ZRT);
    double calcTheta_e(double XRT, double YRT, double ZRT);
    double calcPhi_e(double XRT, double ZRT);
    double calcFirstTerm(double Mach, double L, double D, double D_starred, double dist_obs);
    double calcDh(double Mach, double theta_e, double phi_e,double EddyMach);
    double calcDl(double Mach, double theta_e, double phi_e);
    //Sara

    NoiseParameter *m_parameter;
	
    //For general
    double m_DStarInterpolatedS;
    double m_DStarInterpolatedP;

    double m_DStarFinalS;
    double m_DStarFinalP;
    double m_EddyMachNumber;
    double m_SwAlpha1;
    double m_SwAlpha;
    bool m_AlphaBigSw;
    //Turbulent Inflow
    double m_IntegralLengthScale; //Alexandre MOD
    double m_TurbulenceIntensity; //Alexandre MOD

    bool m_CalcSeparatedFlow;
    bool m_CalcSuctionSide;
    bool m_CalcPressureSide;

    //Sara
    bool m_CalcLowson;
    //Sara

    double m_A1Ar;

    QVector<double> m_OASPL;
    QVector<double> m_OASPLA;
    QVector<double> m_OASPLB;
    QVector<double> m_OASPLC;
    QVector<double> m_SPLALOG;
    QVector<double> m_SPLSLOG;
    QVector<double> m_SPLPLOG;

//Sara
    QVector<double> m_SPLLEdB;
    QVector<double> m_SPLLEdBAW;
    QVector<double> m_SPLLEdBBW;
    QVector<double> m_SPLLEdBCW;
    QVector<double> m_SPLlogLE;

    QVector<double> m_SPLLEdB_rotor;
    QVector<double> m_SPLLEdBAW_rotor;
    QVector<double> m_SPLLEdBBW_rotor;
    QVector<double> m_SPLLEdBCW_rotor;
    QVector<double> m_SPLlogLE_rotor;

    QVector<double> m_OASPL3d;
    QVector<double> m_OASPLA3d;
    QVector<double> m_OASPLB3d;
    QVector<double> m_OASPLC3d;
    QVector<double> m_SPLALOG3d;
    QVector<double> m_SPLSLOG3d;
    QVector<double> m_SPLPLOG3d;
    QVector<double> m_SPLLEdB3d;
    QVector<double> m_SPLLEdBAW3d;
    QVector<double> m_SPLLEdBBW3d;
    QVector<double> m_SPLLEdBCW3d;
    QVector<double> m_SPLlogLE3d;

    QVector<double> m_OASPL3d_rotor;
    QVector<double> m_OASPLA3d_rotor;
    QVector<double> m_OASPLB3d_rotor;
    QVector<double> m_OASPLC3d_rotor;
    QVector<double> m_SPLALOG3d_rotor;
    QVector<double> m_SPLSLOG3d_rotor;
    QVector<double> m_SPLPLOG3d_rotor;
    QVector<double> m_SPLLEdB3d_rotor;
    QVector<double> m_SPLLEdBAW3d_rotor;
    QVector<double> m_SPLLEdBBW3d_rotor;
    QVector<double> m_SPLLEdBCW3d_rotor;
    QVector<double> m_SPLlogLE3d_rotor;
    //Sara

    //For SPLa
    TwoDVector m_SPLadB; //Store db of SPL alpha
    TwoDVector m_SPLadBAW; //Store db of SPL alpha + A-Weighting
    TwoDVector m_SPLadBBW; //Store db of SPL alpha + A-Weighting
    TwoDVector m_SPLadBCW; //Store db of SPL alpha + A-Weighting
    double m_SplaFirstTerm;
    double m_SplaSt1;
    double m_SplaSt2;
    double m_SplaGamma;
    double m_SplaBeta;
    double m_SplaBetaZero;
    double m_SplaK1;
    double m_SplaK2;
    double m_SplaBr;
    double m_SplaBMax;
    double m_SplaBMin;
    double m_SplaBo;
    double m_SplaAr;
    double m_SplaAMax;
    double m_SplaAMin;
    double m_SplaAo;
    double m_ChordBasedReynolds;

    //Sara
    bool m_rot_speed_check;
    bool m_u_wind_speed_check;
    bool m_TSRtd_check;
    double m_sects;
    double m_rot_speed;
    double m_u_wind_speed;
    double m_TSRtd;
    double m_dstar_user;
    double x;
    double m_rot_speed_calc;
    double m_u_wind_speed_calc;
    double m_TSR_calc;
    double m_obs_x_pos;
    double m_obs_y_pos;
    double m_obs_z_pos;
    double m_obs_x_pos_rotor;
    double m_obs_y_pos_rotor;
    double m_obs_z_pos_rotor;
    double m_tower_height;
    double m_tower_to_hub_distance;
    int m_initial_azimuth;
    int m_yaw_angle;
    int m_dstar_type;
    int m_state_ss_us;
    int m_step_type;
    int m_anglesteps;
    int m_phi_type;
    int m_theta_type;
    //Sara

    //For SPLs
    TwoDVector m_SPLsdB; //Store db of SPLs
    TwoDVector m_SPLsdBAW; //Store db of SPLs + A-Weighting
    TwoDVector m_SPLsdBBW; //Store db of SPLs + B-Weighting
    TwoDVector m_SPLsdBCW; //Store db of SPLs + C-Weighting
    double m_SplsFirstTerm;
    double m_SplsSt1;
    double m_SplsSt2;
    double m_SplsK1;
    double m_SplsSt1Bar;
    double m_SplsK13;

    //For SPLp
    TwoDVector m_SPLpdB; //Store db of SPLp
    TwoDVector m_SPLpdBAW; //Store db of SPLp + A-Weighting
    TwoDVector m_SPLpdBBW; //Store db of SPLp + B-Weighting
    TwoDVector m_SPLpdBCW; //Store db of SPLp + C-Weighting
    double m_SplpFirstTerm;
    double m_SplpSt1;
    double m_SplpK1;
    double m_SplpK13;
    double m_SplpDeltaK1;
    double m_ReynoldsBasedDisplacement;

    //For SPL
    TwoDVector m_SPLdB; //Store db of SPL
    TwoDVector m_SPLdBAW; //Store db of SPL + A-Weighting
    TwoDVector m_SPLdBBW; //Store db of SPL + B-Weighting
    TwoDVector m_SPLdBCW; //Store db of SPL + C-Weighting

    //For LE - Alexandre MOD
    TwoDVector m_SPL_LEdB; //Store db of SPL_LE
    TwoDVector m_SPL_LEdBAW; //Store db of SPL_LE + A-Weighting
    TwoDVector m_SPL_LEdBBW; //Store db of SPL_LE + B-Weighting
    TwoDVector m_SPL_LEdBCW; //Store db of SPL_LE + C-Weighting
    double m_originalVelocity;
    double m_originalChordLength;
    double m_distanceObserver;
    double m_originalMach;
    double m_wettedLength;

    //Sara
    TwoDVector m_SPLadB3d;
    TwoDVector m_SPLadBAW3d;
    TwoDVector m_SPLadBBW3d;
    TwoDVector m_SPLadBCW3d;
    TwoDVector m_SPLpdB3d;
    TwoDVector m_SPLpdBAW3d;
    TwoDVector m_SPLpdBBW3d;
    TwoDVector m_SPLpdBCW3d;
    TwoDVector m_SPLdB3d;
    TwoDVector m_SPLdBAW3d;
    TwoDVector m_SPLdBBW3d;
    TwoDVector m_SPLdBCW3d;
    TwoDVector m_SPL_LEdB3d;
    TwoDVector m_SPL_LEdBAW3d;
    TwoDVector m_SPL_LEdBBW3d;
    TwoDVector m_SPL_LEdBCW3d;
    TwoDVector m_SPLsdB3d;
    TwoDVector m_SPLsdBAW3d;
    TwoDVector m_SPLsdBBW3d;
    TwoDVector m_SPLsdBCW3d;

    TwoDVector m_SPLadB3d_rotor;
    TwoDVector m_SPLadBAW3d_rotor;
    TwoDVector m_SPLadBBW3d_rotor;
    TwoDVector m_SPLadBCW3d_rotor;
    TwoDVector m_SPLpdB3d_rotor;
    TwoDVector m_SPLpdBAW3d_rotor;
    TwoDVector m_SPLpdBBW3d_rotor;
    TwoDVector m_SPLpdBCW3d_rotor;
    TwoDVector m_SPLdB3d_rotor;
    TwoDVector m_SPLdBAW3d_rotor;
    TwoDVector m_SPLdBBW3d_rotor;
    TwoDVector m_SPLdBCW3d_rotor;
    TwoDVector m_SPL_LEdB3d_rotor;
    TwoDVector m_SPL_LEdBAW3d_rotor;
    TwoDVector m_SPL_LEdBBW3d_rotor;
    TwoDVector m_SPL_LEdBCW3d_rotor;
    TwoDVector m_SPLsdB3d_rotor;
    TwoDVector m_SPLsdBAW3d_rotor;
    TwoDVector m_SPLsdBBW3d_rotor;
    TwoDVector m_SPLsdBCW3d_rotor;

    TwoDVector m_SPLadB3d_final;
    TwoDVector m_SPLpdB3d_final;
    TwoDVector m_SPLdB3d_final;
    TwoDVector m_SPLdBAW3d_final;
    TwoDVector m_SPLdBBW3d_final;
    TwoDVector m_SPLdBCW3d_final;
    TwoDVector m_SPL_LEdB3d_final;
    TwoDVector m_SPL_LEdBAW3d_final;
    TwoDVector m_SPL_LEdBBW3d_final;
    TwoDVector m_SPL_LEdBCW3d_final;
    TwoDVector m_SPLsdB3d_final;

    TwoDVector m_SPLadB3d_final_rotor;
    TwoDVector m_SPLpdB3d_final_rotor;
    TwoDVector m_SPLdB3d_final_rotor;
    TwoDVector m_SPLdBAW3d_final_rotor;
    TwoDVector m_SPLdBBW3d_final_rotor;
    TwoDVector m_SPLdBCW3d_final_rotor;
    TwoDVector m_SPL_LEdB3d_final_rotor;
    TwoDVector m_SPL_LEdBAW3d_final_rotor;
    TwoDVector m_SPL_LEdBBW3d_final_rotor;
    TwoDVector m_SPL_LEdBCW3d_final_rotor;
    TwoDVector m_SPLsdB3d_final_rotor;

    //for quasi 3d rotor
    TwoDVector m_unsteady_angles;
    //Sara
};

#endif // NOISECALCULATION_H
