#include "NoiseSimulation.h"
#include "NoiseCalculation.h" //Sara
#include "NoiseParameter.h" //Sara

#include "../ParameterViewer.h"
#include "../Store.h"
#include "../Objects/Polar.h"
#include "../Objects/OpPoint.h"
#include "../Graph/NewCurve.h"
#include "NoiseModule.h"
#include "../ColorManager.h"
#include "NoiseOpPoint.h"
#include "../XBEM/BEM.h"
//Sara
#include "../XBEM/TData.h"
#include "../XBEM/TBEMData.h"
#include "../XDMS/DData.h"
#include "../XBEM/Blade.h"
#include "NoiseException.h"
#include "NoiseOpPoint.h"
#include "NoiseCreatorDialog.h"
#include <cmath>
#include <sstream>
#include <string>
//Sara

NoiseSimulation *NoiseSimulation::newBySerialize() {
    NoiseSimulation *simulation = new NoiseSimulation;
    simulation->serialize();
    return simulation;
}

NoiseSimulation::NoiseSimulation(ParameterViewer<Parameter::NoiseSimulation> *viewer) {
    viewer->storeObject(this);
    pen()->setColor(g_colorManager.getLeastUsedColor(&g_noiseSimulationStore));
}

void NoiseSimulation::serialize() {
    StorableObject::serialize();
    ShowAsGraphInterface::serialize();

    m_parameter.serialize();
    m_calculation.serialize();
}

void NoiseSimulation::restorePointers() {
    StorableObject::restorePointers();

    m_parameter.restorePointers();
}

QPen NoiseSimulation::doGetPen(int curveIndex, int highlightedIndex, bool) {
    if (curveIndex == -1 || highlightedIndex == -1) {
        return m_pen;
    } else {
        QPen pen (m_pen);

        if (g_noiseModule->isColorByOpPoint()) {
            pen.setColor (g_colorManager.getColor(curveIndex));
        }

        if (curveIndex == highlightedIndex) {
            pen.setWidth(pen.width()+2);
        }

        return pen;
    }
}

NewCurve *NoiseSimulation::newCurve(QString xAxis, QString yAxis, NewGraph::GraphType /*graphType*/, int opPointIndex) {
    if (xAxis == "" || yAxis == "")
        return NULL;

    bool zeroY = false;
    QVector<double> xVector, yVector;  // because QVector is internally shared there should be no copying
    for (int i = 0; i < 2; ++i) {
        const int index = getAvailableVariables().indexOf(i == 0 ? xAxis : yAxis);
        QVector<double> *vector = (i == 0 ? &xVector : &yVector);
        switch (index) {
        case 0: *vector = NoiseCalculation::CENTRAL_BAND_FREQUENCY; break;
        case 1: *vector = m_calculation.SPLadB()[opPointIndex]; zeroY = true; break;
        case 2: *vector = m_calculation.SPLsdB()[opPointIndex]; zeroY = true; break;
        case 3: *vector = m_calculation.SPLpdB()[opPointIndex]; zeroY = true; break;
        case 4: *vector = m_calculation.SPLdB()[opPointIndex]; break;
        case 5: *vector = m_calculation.SPLdBAW()[opPointIndex]; break;
        case 6: *vector = m_calculation.SPLdBBW()[opPointIndex]; break;
        case 7: *vector = m_calculation.SPLdBCW()[opPointIndex]; break;
        case 8: *vector = m_calculation.SPL_LEdB()[opPointIndex]; break; //Alexandre MOD
        default: return nullptr;
        }
    }

    NewCurve *curve = new NewCurve (this);
//    curve->setAllPoints(xVector.data(), yVector.data(), xVector.size());
    for (int i = 0; i < xVector.size(); ++i) {  // zero the y values lower 0 for certain outputs
        curve->addPoint(xVector[i], (zeroY && yVector[i] < 0 ? 0.0 : yVector[i]));
    }
    return curve;
}

QStringList NoiseSimulation::getAvailableVariables(NewGraph::GraphType /*graphType*/) {
    QStringList variables;

    // WARNING: when changing any variables list, change newCurve as well!
    variables << "Freq [Hz]" << "SPL_alpha" << "SPL_S" << "SPL_P" << "SPL (dB)" << "SPL (dB(A))" << "SPL (dB(B))"
              << "SPL (dB(C))" << "SPL_LE (dB)"; //Alexandre MOD

    return variables;
}

QStringList NoiseSimulation::prepareMissingObjectMessage() {
    if (g_noiseSimulationStore.isEmpty()) {
        QStringList message = CPolar::prepareMissingObjectMessage();
        if (message.isEmpty()) {
            if (g_mainFrame->m_iApp == NOISEMODULE) {
                message = QStringList(">>> Click 'New' to create a new Noise Simulation");
            } else {
                message = QStringList(">>> unknown hint");
            }
        }
        message.prepend("- No Noise Simulation in Database");
        return message;
    } else {
        return QStringList();
    }
}

void NoiseSimulation::simulate() {
    m_calculation.setNoiseParam(&m_parameter);
    m_calculation.calculate();
}

void NoiseSimulation::exportCalculation(QTextStream &stream) {
    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Noise prediction file export" << endl;
    stream << endl;
    if(m_parameter.Lowson_type==1){stream << "leading edge model: Von Kármán" <<endl;}
    if(m_parameter.Lowson_type==2){stream << "leading edge model: Rapid Distortion" <<endl;}//era m_calculation.m_
//    stream << "constante D: " << m_calculation.d_const << endl;
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();
    for (int i = 0; i < noiseOpPoints.size(); ++i) {
        stream << qSetFieldWidth(0);
        stream << "Alpha: " << noiseOpPoints[i]->getAlphaDegree() <<
                  ", Re = " << noiseOpPoints[i]->getReynolds() << endl;
        stream << "OASPL: " << m_calculation.OASPL()[i] << " dB" << endl;
        stream << "OASPL (A): " << m_calculation.OASPLA()[i] << " dB(A)" << endl;
        stream << "OASPL (B): " << m_calculation.OASPLB()[i] << " dB(B)" << endl;
        stream << "OASPL (C): " << m_calculation.OASPLC()[i] << " dB(C)" << endl;
        stream << "SPL_a: " << m_calculation.SPLALOG()[i] << "" << endl;
        stream << "SPL_s: " << m_calculation.SPLSLOG()[i] << "" << endl;
        stream << "SPL_p: " << m_calculation.SPLPLOG()[i] << "" << endl;
        if(m_parameter.Lowson_type!=0){
        stream << "SPL_LE: " << m_calculation.SPLlogLE()[i] << " dB" << endl;
        stream << "SPL_LE (A): " << m_calculation.SPLLEdBAW()[i] << " dB(A)" << endl;
        stream << "SPL_LE (B): " << m_calculation.SPLLEdBBW()[i] << " dB(B)" << endl;
        stream << "SPL_LE (C): " << m_calculation.SPLLEdBCW()[i] << " dB(C)" << endl;}
        stream << endl;
               if(m_parameter.Lowson_type!=0){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" <<
                  "SPL (dB)" <<
                  "SPLa" <<
                  "SPLs" <<
                  "SPLp" <<
                  "SPL (dB(A))" <<
                  "SPL (dB(B))" <<
                  "SPL (dB(C))" <<
                  "SPL_LE (dB)" <<
                  "SPL_LE (dB(A))" <<
                  "SPL_LE (dB(B))" <<
                  "SPL_LE (dB(C))" <<endl; //Alexandre MOD
               }
               else{
                   stream << qSetFieldWidth(14) <<
                             "Freq [Hz]" <<
                             "SPL (dB)" <<
                             "SPLa" <<
                             "SPLs" <<
                             "SPLp" <<
                             "SPL (dB(A))" <<
                             "SPL (dB(B))" <<
                             "SPL (dB(C))" <<endl; //Alexandre MOD
               }

        for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

 if(m_parameter.Lowson_type!=0){
            stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] <<
                      m_calculation.SPLdB()[i][j] <<
                      m_calculation.SPLadB()[i][j] <<
                      m_calculation.SPLsdB()[i][j] <<
                      m_calculation.SPLpdB()[i][j] <<
                      m_calculation.SPLdBAW()[i][j] <<
                      m_calculation.SPLdBBW()[i][j] <<
                      m_calculation.SPLdBCW()[i][j] <<
                      m_calculation.SPL_LEdB()[i][j] <<
                      m_calculation.SPL_LEdBAW()[i][j] << //Sara
                      m_calculation.SPL_LEdBBW()[i][j] << //Sara
                      m_calculation.SPL_LEdBCW()[i][j] << //Sara
                      endl; //Alexandre MOD
        }
else{
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] <<
               m_calculation.SPLdB()[i][j] <<
               m_calculation.SPLadB()[i][j] <<
               m_calculation.SPLsdB()[i][j] <<
               m_calculation.SPLpdB()[i][j] <<
               m_calculation.SPLdBAW()[i][j] <<
               m_calculation.SPLdBBW()[i][j] <<
               m_calculation.SPLdBCW()[i][j] << endl;
 }
        }

        stream << endl;
        stream << endl;
    }
    qDeleteAll(noiseOpPoints);
}

//Sara
void NoiseSimulation::exportqs3DLog(QTextStream &stream) {

    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    const double AWeighting[] = {-44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1,
                                                                       -13.4, -10.9,  -8.6,  -6.6,  -4.8,  -3.2,  -1.9,  -0.8,
                                                                         0.0,   0.6,   1.0,   1.2,   1.3,   1.2,   1.0,   0.5,
                                                                        -0.1,  -1.1,  -2.5,  -4.3,  -6.6,-  9.3};

    const double BWeighting[] = {-20.4, -17.1, -14.2, -11.6,  -9.3,  -7.4,  -5.6,  -4.2,
                                                                        -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.1,   0.0,
                                                                         0.0,   0.0,   0.0,  -0.1,  -0.2,  -0.4,  -0.7,  -1.2,
                                                                        -1.9,  -2.9,  -4.3,  -6.1,  -8.4, -11.1};

    const double CWeighting[] = {-4.4,  -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.2,
                                                                        -0.1,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                                                                         0.0,   0.0,  -0.1,  -0.2,  -0.3,  -0.5,  -0.8,  -1.3,
                                                                        -2.0,  -3.0,  -4.4,  -6.2,  -8.5, -11.2};

    const double Frequency[]= {25, 31.5, 40, 50, 63, 80, 100, 125, 160,200, 250, 315, 400,500, 630, 800, 1000,1250, 1600, 2000, 2500, 3150, 4000, 5000,6300, 8000, 10000, 12500, 16000, 20000};

    stream << "Quasi 3D Noise Log" << endl;
    stream << endl;

QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
double lend  =   pSimuWidget->m_pctrlLELineEdit->getValue();
double z=lstart;
double approaxing_wind_speed = m_parameter.originalVelocity;

QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
double blade_pitch=pbem->m_pctrlFixedPitch->getValue();

    foreach(BData * bdata, pbem->m_pBEMData->GetBData()){
    int number_of_segments = bdata->m_pos.size();
    double rho = pbem->dlg_rho;
    double dynamic_visc = pbem->dlg_visc;
    double cin_visc = dynamic_visc/rho;
    double K_air = 1.4;
    double R_air = 286.9;
    double T_std_cond = pbem->dlg_temp;
    double P_std_cond = 101300;
    double lambda = pbem->dlg_lambda;
    int mpos_size = bdata->m_pos.size(); //total number of segments
    double finalradius = bdata->m_pos.value(mpos_size-1);
    double nom_tg_speed = bdata->windspeed*lambda;
    double omega = nom_tg_speed/finalradius;
    double rotation = 60/(M_PI*100/nom_tg_speed);

    //definitions
    double axial_ind_fact[number_of_segments];
    double axial_ind_fact_n[number_of_segments];
    double axial_velocity[number_of_segments];
    double tangential_speed[number_of_segments];
    double resultant_local_speed[number_of_segments];
    double chord[number_of_segments];
    double Reynolds[number_of_segments];
    double Reynolds_BEM[number_of_segments];
    double Reynolds_polar[number_of_segments];
    double Reynolds_error[number_of_segments];
    double Mach_BEM[number_of_segments];
    double Mach_polar[number_of_segments];
    double Mach_error[number_of_segments];
    double alpha[number_of_segments];
    double alpha_BEM[number_of_segments];
    double alpha_polar[number_of_segments];
    double alpha_error[number_of_segments];
    double phi_BEM[number_of_segments];
    double theta_BEM[number_of_segments];
    double cl_cd[number_of_segments];
    double r_R[number_of_segments];
    double c_Rx[number_of_segments];
    double D_starred_C_HT[number_of_segments];
    double D_starred_HT[number_of_segments];
    double D_starred_C_N[number_of_segments];
    double D_starred_N[number_of_segments];
    double Dh[number_of_segments];
    double Dl[number_of_segments];
    double Mach[number_of_segments];
    double corr_fact[number_of_segments];
    double D_starred_HT_S[number_of_segments];
    double D_starred_HT_P[number_of_segments];
    double D_starred_N_S[number_of_segments];
    double D_starred_N_P[number_of_segments];
    double L[number_of_segments];
    double SwAlpha[number_of_segments];
    double SwAlpha_1[number_of_segments];
    double SwAlpha_2[number_of_segments];
    double observer_position = 10;
    double gamma[number_of_segments];
    double gamma0[number_of_segments];
    double beta[number_of_segments];
    double beta0[number_of_segments];
    double gamma0_gamma_min[number_of_segments];
    double gamma0_gamma_plus[number_of_segments];
    double K1[number_of_segments];
    double K2[number_of_segments];
    double EddyMach_calc[number_of_segments];
    double dist_z[number_of_segments];
    double dist_y[number_of_segments];
    double dist_obs[number_of_segments];
    double D_starred_S[number_of_segments];
    double D_starred_P[number_of_segments];
    double first_term_Dh_S[number_of_segments];
    double first_term_Dl_S[number_of_segments];
    double first_term_Dh_P[number_of_segments];
    double St1[number_of_segments];
    double St2[number_of_segments];
    double b0[number_of_segments];
    double B_min_b0[number_of_segments];
    double B_max_b0[number_of_segments];
    double BR_b0[number_of_segments];
    double RCmod[number_of_segments];
    double ao_Rc[number_of_segments];
    double A_min_ao[number_of_segments];
    double A_max_ao[number_of_segments];
    double K1_3[number_of_segments];
    double AR_ao[number_of_segments];
    double St1_bar[number_of_segments];
    double Re_disp_thick[number_of_segments];
    double delta_K1[number_of_segments];
    double phi_rad[number_of_segments];
    double theta_rad[number_of_segments];
    double b[number_of_segments];
    double a[number_of_segments];
    double XRS[number_of_segments];
    double YRS[number_of_segments];
    double ZRS[number_of_segments];
    double XRT[number_of_segments];
    double YRT[number_of_segments];
    double ZRT[number_of_segments];
    double re[number_of_segments];
    double theta[number_of_segments];
    double phi[number_of_segments];
    double calc_int_a[number_of_segments];
    double psi_e[number_of_segments];
    double theta_e[number_of_segments];
    double r_rt[number_of_segments];
    double r_e[number_of_segments];
    double r_1[number_of_segments];
    double c_1[number_of_segments];
    double r_0[number_of_segments];
    double c_0[number_of_segments];
    double twist[number_of_segments];
    double local_twist[number_of_segments];

    double ri[number_of_segments];
    double ri_1[number_of_segments];
    double ri_2[number_of_segments];
    double A[number_of_segments];
    double phi_lin[number_of_segments];

    double DStarXFoilS[number_of_segments];
    double DStarXFoilP[number_of_segments];

    double sp_OASPL_alpha=0;
    double splog_OASPL_alpha=0;
    double st_OASPL_alpha=0;
    double stlog_OASPL_alpha=0;
    double sp_OASPL_S=0;
    double splog_OASPL_S=0;
    double st_OASPL_S=0;
    double stlog_OASPL_S=0;
    double sp_OASPL_P=0;
    double splog_OASPL_P=0;
    double st_OASPL_P=0;
    double stlog_OASPL_P=0;
    double sp_OASPL=0;
    double splog_OASPL=0;
    double st_OASPL=0;
    double stlog_OASPL=0;
    double sp_dBA=0;
    double splog_dBA=0;
    double st_dBA=0;
    double stlog_dBA=0;
    double sp_dBB=0;
    double splog_dBB=0;
    double st_dBB=0;
    double stlog_dBB=0;
    double sp_dBC=0;
    double splog_dBC=0;
    double st_dBC=0;
    double stlog_dBC=0;
    double SPL_blade_total=0;//Sara new
    double SPL_blade_total_aux=0;//Sara new

    double r_R0  =  0.05; double c_R0 = 0.05500;
    double r_R1  =  0.25; double c_R1 = 0.07500;
    double r_R2  =  1.00; double c_R2 = 0.02000;

        QString str= QString::number(z, 'f', 1);

        if((z<=m_parameter.TSRtd) & (z>=m_parameter.TSRtd)){

        stream << "Tip Speed Ratio: " << str << endl;
        stream << endl;

        stream << qSetFieldWidth(14)  <<
                  "Section"  << ";" <<
                  "Radius [m]"  << ";" <<
                  "r/R"  << ";" <<
                  "Chord [m]" << ";" <<
                  "Re polar"  << ";" <<
                  "Re calc"  << ";" <<
                  "Re error [%]"  << ";" <<
                  "Mach polar"  << ";" <<
                  "Mach calc"  << ";" <<
                  "Mach error [%]"   << ";" <<
                  "AOA polar [deg]"   <<";"   <<
                  "AOA calc [deg]"   <<";"   <<
                  "AOA error [%]"   <<";"   <<
                  "Error number[1]" <<  ";" <<
                  endl;
}

        for (int i = 0; i < number_of_segments; ++i) {

sp_OASPL_alpha=0;
sp_OASPL_S=0;
sp_OASPL_P=0;
sp_OASPL=0;
sp_dBA=0;
sp_dBB=0;
sp_dBC=0;
splog_OASPL_alpha=0;
splog_OASPL_S=0;
splog_OASPL_P=0;
splog_OASPL=0;
splog_dBA=0;
splog_dBB=0;
splog_dBC=0;

            // definitions
            axial_ind_fact[i] = bdata->m_a_axial.value(i);
            axial_ind_fact_n[i] = bdata->m_a_axial.value(i+1);

            if (i<number_of_segments/2) {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact[i]);}
            else {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact_n[i]);}

            tangential_speed[i] = omega*bdata->m_pos.value(i)*(1.+bdata->m_a_tangential.value(i));
            resultant_local_speed[i] = qSqrt(pow(axial_velocity[i],2)+pow(tangential_speed[i],2));
            chord[i] = bdata->m_c_local.value(i);
            Reynolds[i] = bdata->m_Reynolds.value(i);

            CPolar *pCPolar = (CPolar *) g_mainFrame->m_pctrlPolar;
            Reynolds_BEM[i]=bdata->m_Reynolds.value(i);
            Reynolds_polar[i]=noiseOpPoints[i]->getReynolds();
            Reynolds_error[i]=qFabs(Reynolds_polar[i]-Reynolds_BEM[i])/Reynolds_BEM[i]*100.;

            Mach_polar[i]=noiseOpPoints[i]->getMach();
            Mach[i]=bdata->m_Mach.value(i);
            Mach_BEM[i] = bdata->m_Mach.value(i);

            Mach_error[i]=qFabs(Mach_polar[i]-Mach_BEM[i])/Mach_BEM[i]*100.;
            alpha_BEM[i] = bdata->m_alpha.value(i);

            D_starred_N_S[i]=0;
            D_starred_HT_S[i]=0;

double aux_alpha_polar[noiseOpPoints.size()];
for (int k=0;k<noiseOpPoints.size();++k){
aux_alpha_polar[k]=qFabs(alpha_BEM[i]-noiseOpPoints[k]->getAlphaDegreeAbsolute());
}

double aux_alpha_polar_set=aux_alpha_polar[0];

for (int k=1;k<noiseOpPoints.size();++k){
    if(aux_alpha_polar_set>aux_alpha_polar[k])
   {aux_alpha_polar_set=aux_alpha_polar[k];
alpha_polar[i]=noiseOpPoints[k]->getAlphaDegreeAbsolute();
    }
}

            alpha[i]=alpha_BEM[i];

            phi_BEM[i] = bdata->m_phi.value(i);
            theta_BEM[i] = bdata->m_theta.value(i);
            cl_cd[i] =  bdata->m_LD.value(i);
            r_R[i] = bdata->m_pos.value(i)/finalradius;

            if (r_R[i] <= r_R0) {c_Rx[i] = c_R0;}
            if ((r_R[i] > r_R0) & (r_R[i] < r_R1)) {c_Rx[i] = (r_R[i]-r_R0)*(c_R1-c_R0)/(0.25-r_R0)+c_R0;}
            if ((r_R[i] <= r_R1) & (r_R[i] >= r_R1)) {c_Rx[i] = c_R1;}
            if ((r_R[i] > r_R1) & (r_R[i] < r_R2)) {c_Rx[i] = (r_R[i]-r_R1)*(c_R2-c_R1)/(r_R2-r_R1)+c_R1;}
            if (r_R[i] >= r_R2) {c_Rx[i] = c_R2;}

            QString c_R= QString::number(c_Rx[i], 'f', 5);
            QString Mach_error_x= QString::number(Mach_error[i], 'f', 2);
            QString Reynolds_error_x= QString::number(Reynolds_error[i], 'f', 2);
            QString alpha_error_x= QString::number(alpha_error[i], 'f', 2);

//            Mach[i]=Mach_calc[i];

//heavy tripping
if ((alpha[i]<=0) & (alpha[i]>=0)){
if (Reynolds[i]>300000){
    D_starred_C_HT[i]=pow(10,(3.411-1.5397*log10(Reynolds[i])+0.1059*pow(log10(Reynolds[i]),2)));
}
else {D_starred_C_HT[i]=0.0601*(pow(Reynolds[i],(-0.114)));}

D_starred_HT[i]=chord[i]*D_starred_C_HT[i];

//natural transition
    D_starred_C_N[i]=pow(10,(3.0187-1.5397*log10(Reynolds[i])+0.1059*pow(log10(Reynolds[i]),2)));

D_starred_N[i]=D_starred_C_N[i]*chord[i];

D_starred_HT_S[i]=D_starred_HT[i];
D_starred_HT_P[i]=D_starred_HT[i];
D_starred_N_S[i]=D_starred_N[i];
D_starred_N_P[i]=D_starred_N[i];
}
else{
//alpha !=0 pressure side
if ((alpha[i]<0) & (alpha[i]>0)){
corr_fact[i]=pow(10,(-0.0432*alpha[i]+0.00113*pow(alpha[i],2)));
D_starred_HT_P[i]=D_starred_HT[i]*corr_fact[i];
D_starred_N_P[i]=D_starred_N[i]*corr_fact[i];
}

//alpha !=0 suction side heavy tripping
if ((alpha[i]>0) & (alpha[i]<=5)){
corr_fact[i]=pow(10,(0.0679*alpha[i]));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

if ((alpha[i]>5) & (alpha[i]<=12.5)){
corr_fact[i]=0.381*(pow(10,(0.1516*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

if ((alpha[i]>12.5) & (alpha[i]<=25)){
corr_fact[i]=14.296*(pow(10,(0.0258*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

//alpha !=0 suction side natural transition
if ((alpha[i]>0) & (alpha[i]<=7.5)){
corr_fact[i]=pow(10,(0.0679*alpha[i]));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

if ((alpha[i]>7.5) & (alpha[i]<=12.5)){
corr_fact[i]=0.0162*(pow(10,(0.3066*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

if ((alpha[i]>12.5) & (alpha[i]<=25)){
corr_fact[i]=54.42*(pow(10.,(0.0258*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}
}

//For D* Xfoil

DStarXFoilS[i]=m_calculation.m_DStarInterpolatedS3d[i];
DStarXFoilP[i]=m_calculation.m_DStarInterpolatedP3d[i];

//Length of Wetted Trailing Edge
L[i]=bdata->m_pos.value(i+1)-bdata->m_pos.value(i);

//Calculate the Switching Angle
SwAlpha_1[i]=23.43*Mach[i]+4.651;
SwAlpha_2[i]=12.5;

if (SwAlpha_1[i]<SwAlpha_2[i]){SwAlpha[i]=SwAlpha_1[i];}
else {SwAlpha[i]=SwAlpha_2[i];}

double EddyMach = m_parameter.eddyConvectionMach;

Dh[i]=(2.*pow(sin(qDegreesToRadians(theta[i]/2.)),2)*pow(sin(qDegreesToRadians(phi[i])),2))/pow(1+Mach[i]*cos(qDegreesToRadians(theta[i]))*(1.+(Mach[i]-Mach[i]*EddyMach)*cos(qDegreesToRadians(phi[i]))),2);

Dl[i]=(2.*pow(sin(qDegreesToRadians(theta[i])),2)*pow(sin(qDegreesToRadians(phi[i])),2))/pow((1+(Mach[i]*cos(qDegreesToRadians(theta[i])))),4);

gamma[i]=27.094*Mach[i]+3.32;

gamma0[i]=SwAlpha_1[i];

beta[i]=72.65*Mach[i]+10.74;

beta0[i]=-34.19*Mach[i]-13.82;

gamma0_gamma_min[i]=gamma0[i]-gamma[i];

gamma0_gamma_plus[i]=gamma0[i]+gamma[i];

if (bdata->m_Reynolds.value(i)<247000)
{K1[i]=-4.31*log10(bdata->m_Reynolds.value(i))+156.3;}
else if (bdata->m_Reynolds.value(i)>800000)
{K1[i]=128.5;}
else {K1[i]=-9.*log10(bdata->m_Reynolds.value(i))+181.6;}

if (alpha[i]<gamma0_gamma_min[i])
{K2[i]=K1[i]-1000.;}
else if (alpha[i]>gamma0_gamma_plus[i])
{K2[i]=K1[i]-12.;}
else
{K2[i]=K1[i]+(sqrt(pow(beta[i],2)-pow((beta[i]/gamma[i]),2)*pow((alpha[i]-gamma0[i]),2)))+beta0[i];}

double EddyMach_perc=EddyMach;

EddyMach_calc[i]=Mach[i]*EddyMach_perc;

//delta starred type, if natural transition or heavy-tripping
if (m_parameter.dstar_type==0){
    D_starred_S[i]=D_starred_N_S[i];
    D_starred_P[i]=D_starred_N_P[i];
}

else if (m_parameter.dstar_type==1){
FoilPolarDlg *pFoilPolarDlg = (FoilPolarDlg *) g_mainFrame->m_pctrlXDirectWidget;

    double TopTrip=pFoilPolarDlg->m_XTopTr;
    double BotTrip=pFoilPolarDlg->m_XBotTr;

    if((TopTrip<=0 & TopTrip>=0) & (BotTrip<=0 & BotTrip>=0)) {
//        natural transition
    D_starred_S[i]=D_starred_HT_S[i];
    D_starred_P[i]=D_starred_HT_P[i];
}
else {
        //heavy tripping
            D_starred_S[i]=DStarXFoilS[i];
            D_starred_P[i]=DStarXFoilP[i];
        }}
        else if (m_parameter.dstar_type==2){
    //user
    NoiseParameter *pNoiseParameter = (NoiseParameter *) g_mainFrame->m_pSimuWidget;
    D_starred_S[i]=pNoiseParameter->D_starred_S_user[i];
    D_starred_P[i]=pNoiseParameter->D_starred_P_user[i];
        }

double B=0;
double XB=0;
double YB=0;
double ZB=0;

//phi type, fixed 90º or free
if (m_parameter.phi_type==0){
    //by the quasi-3D spreadsheet
    ri[i]=bdata->m_pos.value(i);
    ri_1[i]=bdata->m_pos.value(i+1);
    ri_2[i]=bdata->m_pos.value(i+2);
    B=bdata->m_pos.value(number_of_segments-1)/2.;
    A[i]=ri_1[i]+(ri_2[i]-ri_1[i])/2.;
    re[i]=sqrt(pow((A[i]-B),2)+pow(m_parameter.distanceObsever,2));
    phi_lin[i]=qRadiansToDegrees(qAsin(m_parameter.distanceObsever/re[i]));
    phi[i]=180.-phi_lin[i];
    theta[i]=90;
    dist_obs[i]=ri_1[i];
}
else if (m_parameter.phi_type==1){
    //    re phi and theta calculation p 77 C_Project_Log_Text_Jan_16.pdf
    //    Input X e , Y e , Z e
    //    Attribute their respective values to X B , Y B , Z B
        XB=m_parameter.obs_x_pos;
        YB=m_parameter.obs_y_pos;
        ZB=m_parameter.obs_z_pos;

    //For the blade, get the Pitch angle, θp - blade pitch

    //For each blade segment, get: C r i+1 , C r i , r i+1 , r i , β
    //Cri+1 = c_1
    //r+1 = r_1
    //Cri = c_0
    //ri = r_0
    //beta = local_twist

    if((i<=number_of_segments) & (i>=number_of_segments)){
c_1[i]=0;
r_1[i]=0;
    }
    else
    {
c_1[i]=bdata->m_c_local.value(i+1);
r_1[i]=bdata->m_pos.value(i+1);
    }}

c_0[i]=bdata->m_c_local.value(i);
r_0[i]=bdata->m_pos.value(i);

local_twist[i]=theta_BEM[i];

    b[i]=qRadiansToDegrees(qAtan((c_1[i]-c_0[i])/(r_1[i]-r_0[i])));

//    the angle a is the total angle between the Y B Z B blade reference system plane and the local midsection chord line
    a[i]=local_twist[i]+blade_pitch;

    XRS[i]=XB*cos(qDegreesToRadians(a[i]))+YB*sin(qDegreesToRadians(a[i]));
    YRS[i]=-XB*sin(qDegreesToRadians(a[i]))+YB*cos(qDegreesToRadians(a[i]));
    ZRS[i]=ZB-(bdata->m_pos.value(i)-bdata->m_pos.value(i-1))/2.;

//    Calculate Y RS − 0.75 ∗ (C r i+1 − C r i )/2
    calc_int_a[i]=(YRS[i]-0.75*(c_1[i]-c_0[i])/2.);

    XRT[i]=XRS[i];
    YRT[i]=cos(qDegreesToRadians(b[i]))*calc_int_a[i]+sin(qDegreesToRadians(b[i]))*ZRS[i];
    ZRT[i]=-sin(qDegreesToRadians(b[i]))*calc_int_a[i]+cos(qDegreesToRadians(b[i]))*ZRS[i];

    r_e[i]=sqrt(pow(XRT[i],2)+pow(YRT[i],2)+pow(ZRT[i],2));
    r_rt[i]=r_e[i];
    theta_e[i]=qRadiansToDegrees(qAtan(ZRT[i]/YRT[i]));

    psi_e[i]=qRadiansToDegrees(qAtan(XRT[i]/ZRT[i]));

    dist_obs[i]=re[i];

   phi_rad[i]=qDegreesToRadians(phi[i]);
   theta_rad[i]=qDegreesToRadians(theta[i]);

   alpha_error[i]=qFabs(alpha_polar[i]-alpha_BEM[i])/alpha_BEM[i]*100.;

QString observations_x("");
if (!((Reynolds[i] >9.5*pow(10,5)) & (Reynolds[i]<2.5*pow(10,6))) || Mach[i]>=0.19 || !((alpha[i]<=0) & (alpha[i]>=0)) || !((alpha[i]<=5) & (alpha[i]>=5)) || !(alpha[i]<=10 & alpha[i]>=10)){
observations_x.append("1 ");}
if((Reynolds[i] <=4.8*pow(10,4)) || (Reynolds[i]>=2.5*pow(10,6)) || Mach[i]<0.208 || !((alpha[i]<=0) & (alpha[i]>=0))){
observations_x.append("2 ");}
if (Reynolds[i] >=3*pow(10,6) || Mach[i]<0.208 || alpha[i]<=19.8){
observations_x.append("3 ");}
if (Reynolds[i] >1.5*pow(10,6) || Mach[i]<0.208 || alpha[i]<=19.8){
observations_x.append('4');}

//uncomment to input data
 if((z<=m_parameter.TSRtd) & (z>=m_parameter.TSRtd)){
        stream << qSetFieldWidth(14)  <<
                      (i+1) << ";" <<
                      bdata->m_pos.value(i) << ";" <<
                      r_R[i] << ";" <<
                      chord[i] << ";" <<
                      Reynolds_polar[i] << ";" <<
                      Reynolds_BEM[i]  << ";" <<
                      Reynolds_error_x  << ";" <<
                      Mach_polar[i] << ";" <<
                      Mach_BEM[i]  <<  ";" <<
                      Mach_error_x  << ";" <<
                      alpha_polar[i] << ";" <<
                      alpha_BEM[i]   <<  ";" <<
                      alpha_error_x    <<  ";" <<
                  observations_x   <<  ";" <<
                      endl;

}}
         z=z+ldelta;

}
     stream << endl;
     stream << "[1]Out of range in accordance of: 1-Brooks & Hodgson 1981. 2 - Brooks & Marcolini 1985. 3 - Brooks & Marcolini 1986. 4 - Brooks, Pope & Marcolini 1989." << endl;
}

void NoiseSimulation::exportqs3DCalculationComplete(QTextStream &stream)
{
        stream.setRealNumberNotation(QTextStream::FixedNotation);
        stream.setRealNumberPrecision(5);

        const double AWeighting[] = {-44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1,
                                                                           -13.4, -10.9,  -8.6,  -6.6,  -4.8,  -3.2,  -1.9,  -0.8,
                                                                             0.0,   0.6,   1.0,   1.2,   1.3,   1.2,   1.0,   0.5,
                                                                            -0.1,  -1.1,  -2.5,  -4.3,  -6.6,-  9.3};

        const double BWeighting[] = {-20.4, -17.1, -14.2, -11.6,  -9.3,  -7.4,  -5.6,  -4.2,
                                                                            -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.1,   0.0,
                                                                             0.0,   0.0,   0.0,  -0.1,  -0.2,  -0.4,  -0.7,  -1.2,
                                                                            -1.9,  -2.9,  -4.3,  -6.1,  -8.4, -11.1};

        const double CWeighting[] = {-4.4,  -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.2,
                                                                            -0.1,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                                                                             0.0,   0.0,  -0.1,  -0.2,  -0.3,  -0.5,  -0.8,  -1.3,
                                                                            -2.0,  -3.0,  -4.4,  -6.2,  -8.5, -11.2};

        const double Frequency[]= {25, 31.5, 40, 50, 63, 80, 100, 125, 160,200, 250, 315, 400,500, 630, 800, 1000,1250, 1600, 2000, 2500, 3150, 4000, 5000,6300, 8000, 10000, 12500, 16000, 20000};

        stream << "Quasi 3D Noise Log" << endl;
        stream << endl;

NoiseCreatorDialog *pnoisecreatordialog = (NoiseCreatorDialog *) g_mainFrame->m_pSimuWidget;

    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

    SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
    double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
    double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
    double lend  =   pSimuWidget->m_pctrlLELineEdit->getValue();
    double z=lstart;
    double approaxing_wind_speed = m_parameter.originalVelocity;

    QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
    double blade_pitch=pbem->m_pctrlFixedPitch->getValue();

        foreach(BData * bdata, pbem->m_pBEMData->GetBData()){
        int number_of_segments = bdata->m_pos.size();
        double rho = pbem->dlg_rho;
        double dynamic_visc = pbem->dlg_visc;
        double cin_visc = dynamic_visc/rho;
        double K_air = 1.4;
        double R_air = 286.9;
        double T_std_cond = pbem->dlg_temp;
        double P_std_cond = 101300;
        double lambda = pbem->dlg_lambda;
        int mpos_size = bdata->m_pos.size(); //total number of segments
        double finalradius = bdata->m_pos.value(mpos_size-1);
        double nom_tg_speed = bdata->windspeed*lambda;
        double omega = nom_tg_speed/finalradius;
        double rotation = 60/(M_PI*100/nom_tg_speed);
        double c_0_le = 34000;
        double c_const_rd_le = 19./6.;
        double d_const_rd_le = 65.95;
        double c_const_vk_le = 7./3.;
        double d_const_vk_le = 58.4;
        double c_const_le=0;
        double d_const_le=0;
        double aux_le;
        double aux_1_le;
        double aux_4_le;
        double aux_5_le;
        double u_le;
        double c_le;
        double I_le;
        double lambda_le;
        double r_e_le;
        double beta_le;
        double L_le;
        double D_L_le;
        double K_le;
        double S_le;
        double LFC_le;
        double aux0_le;
        double aux1_le;
        double aux4_le;
        double aux5_le;

        QString SPL_alpha_val("");
        QString SPL_dB_S_val("");
        QString SPL_dB_P_val("");
        QString SPL_S_val("");
        QString SPL_P_val("");
        QString SPL_dB_val("");
        QString SPL_A_val("");
        QString SPL_B_val("");
        QString SPL_C_val("");
        QString SPL_LedB_val("");
        QString SPL_LedB_valAW("");
        QString SPL_LedB_valBW("");
        QString SPL_LedB_valCW("");

        //definitions
        double axial_ind_fact[number_of_segments];
        double axial_ind_fact_n[number_of_segments];
        double axial_velocity[number_of_segments];
        double tangential_speed[number_of_segments];
        double resultant_local_speed[number_of_segments];
        double chord[number_of_segments];
        double Reynolds[number_of_segments];
        double Reynolds_BEM[number_of_segments];
        double Reynolds_polar[number_of_segments];
        double Reynolds_error[number_of_segments];
        double Mach_BEM[number_of_segments];
        double Mach_polar[number_of_segments];
        double Mach_error[number_of_segments];
        double alpha[number_of_segments];
        double alpha_BEM[number_of_segments];
        double alpha_polar[number_of_segments];
        double alpha_error[number_of_segments];
        double phi_BEM[number_of_segments];
        double theta_BEM[number_of_segments];
        double cl_cd[number_of_segments];
        double r_R[number_of_segments];
        double c_Rx[number_of_segments];
        double D_starred_C_HT[number_of_segments];
        double D_starred_HT[number_of_segments];
        double D_starred_C_N[number_of_segments];
        double D_starred_N[number_of_segments];
        double Dh[number_of_segments];
        double Dl[number_of_segments];
        double Mach[number_of_segments];
        double corr_fact[number_of_segments];
        double D_starred_HT_S[number_of_segments];
        double D_starred_HT_P[number_of_segments];
        double D_starred_N_S[number_of_segments];
        double D_starred_N_P[number_of_segments];
        double L[number_of_segments];
        double SwAlpha[number_of_segments];
        double SwAlpha_1[number_of_segments];
        double SwAlpha_2[number_of_segments];
        double observer_position = 10;
        double gamma[number_of_segments];
        double gamma0[number_of_segments];
        double beta[number_of_segments];
        double beta0[number_of_segments];
        double gamma0_gamma_min[number_of_segments];
        double gamma0_gamma_plus[number_of_segments];
        double K1[number_of_segments];
        double K2[number_of_segments];
        double EddyMach_calc[number_of_segments];
        double dist_z[number_of_segments];
        double dist_y[number_of_segments];
        double dist_obs[number_of_segments];
        double D_starred_S[number_of_segments];
        double D_starred_P[number_of_segments];
        double first_term_Dh_S[number_of_segments];
        double first_term_Dl_S[number_of_segments];
        double first_term_Dh_P[number_of_segments];
        double St1[number_of_segments];
        double St2[number_of_segments];
        double b0[number_of_segments];
        double B_min_b0[number_of_segments];
        double B_max_b0[number_of_segments];
        double BR_b0[number_of_segments];
        double RCmod[number_of_segments];
        double ao_Rc[number_of_segments];
        double A_min_ao[number_of_segments];
        double A_max_ao[number_of_segments];
        double K1_3[number_of_segments];
        double AR_ao[number_of_segments];
        double St1_bar[number_of_segments];
        double Re_disp_thick[number_of_segments];
        double delta_K1[number_of_segments];
        double phi_rad[number_of_segments];
        double theta_rad[number_of_segments];
        double b[number_of_segments];
        double a[number_of_segments];
        double XRS[number_of_segments];
        double YRS[number_of_segments];
        double ZRS[number_of_segments];
        double XRT[number_of_segments];
        double YRT[number_of_segments];
        double ZRT[number_of_segments];
        double re[number_of_segments];
        double theta[number_of_segments];
        double phi[number_of_segments];
        double calc_int_a[number_of_segments];
        double psi_e[number_of_segments];
        double theta_e[number_of_segments];
        double r_rt[number_of_segments];
        double r_e[number_of_segments];
        double r_1[number_of_segments];
        double c_1[number_of_segments];
        double r_0[number_of_segments];
        double c_0[number_of_segments];
        double twist[number_of_segments];
        double local_twist[number_of_segments];

        double ri[number_of_segments];
        double ri_1[number_of_segments];
        double ri_2[number_of_segments];
        double A[number_of_segments];
        double phi_lin[number_of_segments];

        double DStarXFoilS[number_of_segments];
        double DStarXFoilP[number_of_segments];

        double sp_OASPL_alpha=0;
        double splog_OASPL_alpha=0;
        double st_OASPL_alpha=0;
        double stlog_OASPL_alpha=0;
        double sp_OASPL_S=0;
        double splog_OASPL_S=0;
        double st_OASPL_S=0;
        double stlog_OASPL_S=0;
        double sp_OASPL_P=0;
        double splog_OASPL_P=0;
        double st_OASPL_P=0;
        double stlog_OASPL_P=0;
        double sp_OASPL=0;
        double splog_OASPL=0;
        double st_OASPL=0;
        double stlog_OASPL=0;
        double sp_dBA=0;
        double splog_dBA=0;
        double st_dBA=0;
        double stlog_dBA=0;
        double sp_dBB=0;
        double splog_dBB=0;
        double st_dBB=0;
        double stlog_dBB=0;
        double sp_dBC=0;
        double splog_dBC=0;
        double st_dBC=0;
        double stlog_dBC=0;
        double sp_LedB=0;
        double splog_LedB=0;
        double st_LedB=0;
        double stlog_LedB=0;
        double sp_LedBAW=0;
        double splog_LedBAW=0;
        double st_LedBAW=0;
        double stlog_LedBAW=0;
        double sp_LedBBW=0;
        double splog_LedBBW=0;
        double st_LedBBW=0;
        double stlog_LedBBW=0;
        double sp_LedBCW=0;
        double splog_LedBCW=0;
        double st_LedBCW=0;
        double stlog_LedBCW=0;
        double SPL_blade_total=0;
        double SPL_blade_total_aux=0;

        double r_R0  =  0.05; double c_R0 = 0.05500;
        double r_R1  =  0.25; double c_R1 = 0.07500;
        double r_R2  =  1.00; double c_R2 = 0.02000;

        bool LE_validation;
        bool BPM_validation;

        sp_OASPL_alpha=0;
        sp_OASPL_S=0;
        sp_OASPL_P=0;
        sp_OASPL=0;
        sp_dBA=0;
        sp_dBB=0;
        sp_dBC=0;
        splog_OASPL_alpha=0;
        splog_OASPL_S=0;
        splog_OASPL_P=0;
        splog_OASPL=0;
        splog_dBA=0;
        splog_dBB=0;
        splog_dBC=0;
        splog_LedB=0;
        splog_LedBAW=0;
        splog_LedBBW=0;
        splog_LedBCW=0;

        double sp_OASPL_alpha_a=0;
        double splog_OASPL_alpha_a=0;
        double st_OASPL_alpha_a=0;
        double stlog_OASPL_alpha_a=0;
        double sp_OASPL_S_a=0;
        double splog_OASPL_S_a=0;
        double st_OASPL_S_a=0;
        double stlog_OASPL_S_a=0;
        double sp_OASPL_P_a=0;
        double splog_OASPL_P_a=0;
        double st_OASPL_P_a=0;
        double stlog_OASPL_P_a=0;
        double sp_OASPL_a=0;
        double splog_OASPL_a=0;
        double st_OASPL_a=0;
        double stlog_OASPL_a=0;
        double sp_dBA_a=0;
        double splog_dBA_a=0;
        double st_dBA_a=0;
        double stlog_dBA_a=0;
        double sp_dBB_a=0;
        double splog_dBB_a=0;
        double st_dBB_a=0;
        double stlog_dBB_a=0;
        double sp_dBC_a=0;
        double splog_dBC_a=0;
        double st_dBC_a=0;
        double stlog_dBC_a=0;
        double sp_LedB_a=0;
        double splog_LedB_a=0;
        double st_LedB_a=0;
        double stlog_LedB_a=0;
        double sp_LedBAW_a=0;
        double splog_LedBAW_a=0;
        double st_LedBAW_a=0;
        double stlog_LedBAW_a=0;
        double sp_LedBBW_a=0;
        double splog_LedBBW_a=0;
        double st_LedBBW_a=0;
        double stlog_LedBBW_a=0;
        double sp_LedBCW_a=0;
        double splog_LedBCW_a=0;
        double st_LedBCW_a=0;
        double stlog_LedBCW_a=0;
        double SPL_blade_total_a=0;
        double SPL_blade_total_aux_a=0;

        for (int i = 0; i < number_of_segments; ++i) {

            // definitions
            axial_ind_fact[i] = bdata->m_a_axial.value(i);
            axial_ind_fact_n[i] = bdata->m_a_axial.value(i+1);

            if (i<number_of_segments/2) {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact[i]);}
            else {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact_n[i]);}

            tangential_speed[i] = omega*bdata->m_pos.value(i)*(1.+bdata->m_a_tangential.value(i));
            resultant_local_speed[i] = qSqrt(pow(axial_velocity[i],2)+pow(tangential_speed[i],2));
            chord[i] = bdata->m_c_local.value(i);
            Reynolds[i] = bdata->m_Reynolds.value(i);

            CPolar *pCPolar = (CPolar *) g_mainFrame->m_pctrlPolar;
            Reynolds_BEM[i]=bdata->m_Reynolds.value(i);
            Reynolds_polar[i]=noiseOpPoints[i]->getReynolds();
            Reynolds_error[i]=qFabs(Reynolds_polar[i]-Reynolds_BEM[i])/Reynolds_BEM[i]*100.;
            Mach_polar[i]=noiseOpPoints[i]->getMach();
            Mach[i]=bdata->m_Mach.value(i);
            Mach_BEM[i] = bdata->m_Mach.value(i);
            Mach_error[i]=qFabs(Mach_polar[i]-Mach_BEM[i])/Mach_BEM[i]*100.;
            alpha_BEM[i] = bdata->m_alpha.value(i);
            alpha[i]=alpha_BEM[i];
            phi_BEM[i] = bdata->m_phi.value(i);
            theta_BEM[i] = bdata->m_theta.value(i);
            cl_cd[i] =  bdata->m_LD.value(i);
            r_R[i] = bdata->m_pos.value(i)/finalradius;

            D_starred_N_S[i]=0;
            D_starred_HT_S[i]=0;

//BPM method p 25 manual and spreadsheet
double aux_alpha_polar[noiseOpPoints.size()];
for (int k=0;k<noiseOpPoints.size();++k){
aux_alpha_polar[k]=qFabs(alpha_BEM[i]-noiseOpPoints[k]->getAlphaDegreeAbsolute());
}

double aux_alpha_polar_set=aux_alpha_polar[0];

for (int k=1;k<noiseOpPoints.size();++k){
    if(aux_alpha_polar_set>aux_alpha_polar[k])
   {aux_alpha_polar_set=aux_alpha_polar[k];
alpha_polar[i]=noiseOpPoints[k]->getAlphaDegreeAbsolute();
    }
}
            if (r_R[i] <= r_R0) {c_Rx[i] = c_R0;}
            if ((r_R[i] > r_R0) & (r_R[i] < r_R1)) {c_Rx[i] = (r_R[i]-r_R0)*(c_R1-c_R0)/(0.25-r_R0)+c_R0;}
            if ((r_R[i] <= r_R1) & (r_R[i] >= r_R1)) {c_Rx[i] = c_R1;}
            if ((r_R[i] > r_R1) & (r_R[i] < r_R2)) {c_Rx[i] = (r_R[i]-r_R1)*(c_R2-c_R1)/(r_R2-r_R1)+c_R1;}
            if (r_R[i] >= r_R2) {c_Rx[i] = c_R2;}

            QString c_R= QString::number(c_Rx[i], 'f', 5);
            QString Mach_error_x= QString::number(Mach_error[i], 'f', 2);
            QString Reynolds_error_x= QString::number(Reynolds_error[i], 'f', 2);
            QString alpha_error_x= QString::number(alpha_error[i], 'f', 2);

//heavy tripping

if (Reynolds[i]>300000){
    D_starred_C_HT[i]=pow(10,(3.411-1.5397*log10(Reynolds[i])+0.1059*pow(log10(Reynolds[i]),2)));
}
else {D_starred_C_HT[i]=0.0601*(pow(Reynolds[i],(-0.114)));}

D_starred_HT[i]=chord[i]*D_starred_C_HT[i];

//natural transition
    D_starred_C_N[i]=pow(10,(3.0187-1.5397*log10(Reynolds[i])+0.1059*pow(log10(Reynolds[i]),2)));

D_starred_N[i]=D_starred_C_N[i]*chord[i];

if ((alpha[i]<=0) & (alpha[i]>=0)){
D_starred_HT_S[i]=D_starred_HT[i];
D_starred_HT_P[i]=D_starred_HT[i];
D_starred_N_S[i]=D_starred_N[i];
D_starred_N_P[i]=D_starred_N[i];
}
else{
//alpha !=0 pressure side
if ((alpha[i]<0) || (alpha[i]>0)){
corr_fact[i]=pow(10,(-0.0432*alpha[i]+0.00113*pow(alpha[i],2)));
D_starred_HT_P[i]=D_starred_HT[i]*corr_fact[i];
D_starred_N_P[i]=D_starred_N[i]*corr_fact[i];
}

//alpha !=0 suction side heavy tripping
if ((alpha[i]>0) & (alpha[i]<=5)){
corr_fact[i]=pow(10,(0.0679*alpha[i]));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

if ((alpha[i]>5) & (alpha[i]<=12.5)){
corr_fact[i]=0.381*(pow(10,(0.1516*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

if ((alpha[i]>12.5) & (alpha[i]<=25)){
corr_fact[i]=14.296*(pow(10,(0.0258*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

//alpha !=0 suction side natural transition
if ((alpha[i]>0) & (alpha[i]<=7.5)){
corr_fact[i]=pow(10.,(0.0679*alpha[i]));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

if ((alpha[i]>7.5) & (alpha[i]<=12.5)){
corr_fact[i]=0.0162*(pow(10.,(0.3066*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

if ((alpha[i]>12.5) & (alpha[i]<=25)){
corr_fact[i]=54.42*(pow(10.,(0.0258*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}}

//For D* Xfoil

DStarXFoilS[i]=m_calculation.m_DStarInterpolatedS3d[i];
DStarXFoilP[i]=m_calculation.m_DStarInterpolatedP3d[i];

//Length of Wetted Trailing Edge
L[i]=bdata->m_pos.value(i+1)-bdata->m_pos.value(i);

//Calculate the Switching Angle
SwAlpha_1[i]=23.43*Mach[i]+4.651;
SwAlpha_2[i]=12.5;

if (SwAlpha_1[i]<SwAlpha_2[i]){SwAlpha[i]=SwAlpha_1[i];}
else {SwAlpha[i]=SwAlpha_2[i];}

double EddyMach = m_parameter.eddyConvectionMach;

if (m_parameter.lowFreq) {
Dh[i]=(2.*pow(sin(qDegreesToRadians(theta[i]/2.)),2)*pow(sin(qDegreesToRadians(phi[i])),2))/pow(1+Mach[i]*cos(qDegreesToRadians(theta[i]))*(1.+(Mach[i]-Mach[i]*EddyMach)*cos(qDegreesToRadians(phi[i]))),2);
}
else{Dh[i]=1;}

if (m_parameter.lowFreq) {
Dl[i]=(2.*pow(sin(qDegreesToRadians(theta[i])),2)*pow(sin(qDegreesToRadians(phi[i])),2))/pow((1+(Mach[i]*cos(qDegreesToRadians(theta[i])))),4);
}
else{Dl[i]=1;}

gamma[i]=27.094*Mach[i]+3.32;

gamma0[i]=SwAlpha_1[i];

beta[i]=72.65*Mach[i]+10.74;

beta0[i]=-34.19*Mach[i]-13.82;

gamma0_gamma_min[i]=gamma0[i]-gamma[i];

gamma0_gamma_plus[i]=gamma0[i]+gamma[i];

if (bdata->m_Reynolds.value(i)<247000)
{K1[i]=-4.31*log10(bdata->m_Reynolds.value(i))+156.3;}
else if (bdata->m_Reynolds.value(i)>800000)
{K1[i]=128.5;}
else {K1[i]=-9.*log10(bdata->m_Reynolds.value(i))+181.6;}

if (alpha[i]<gamma0_gamma_min[i])
{K2[i]=K1[i]-1000.;}
else if (alpha[i]>gamma0_gamma_plus[i])
{K2[i]=K1[i]-12.;}
else
{K2[i]=K1[i]+(sqrt(pow(beta[i],2)-pow((beta[i]/gamma[i]),2)*pow((alpha[i]-gamma0[i]),2)))+beta0[i];}

double EddyMach_perc=EddyMach;

EddyMach_calc[i]=Mach[i]*EddyMach_perc;

//delta starred type, if natural transition or heavy-tripping
if (m_parameter.dstar_type==0){
    D_starred_S[i]=D_starred_N_S[i];
    D_starred_P[i]=D_starred_N_P[i];
}

else if (m_parameter.dstar_type==1){
FoilPolarDlg *pFoilPolarDlg = (FoilPolarDlg *) g_mainFrame->m_pctrlXDirectWidget;

    double TopTrip=pFoilPolarDlg->m_XTopTr;
    double BotTrip=pFoilPolarDlg->m_XBotTr;

if((TopTrip<=0 & TopTrip>=0) & (BotTrip<=0 & BotTrip>=0)) {
//        natural transition
    D_starred_S[i]=D_starred_HT_S[i];
    D_starred_P[i]=D_starred_HT_P[i];
}
else {
//heavy tripping
    D_starred_S[i]=DStarXFoilS[i];
    D_starred_P[i]=DStarXFoilP[i];
}}
else if (m_parameter.dstar_type==2){
//user
    NoiseParameter *pNoiseParameter = (NoiseParameter *) g_mainFrame->m_pSimuWidget;
    D_starred_S[i]=pNoiseParameter->D_starred_S_user[i];
    D_starred_P[i]=pNoiseParameter->D_starred_P_user[i];
}

double B=0;
double XB=0;
double YB=0;
double ZB=0;

//phi type, fixed 90º or free
if (m_parameter.phi_type==0){
    //by the quasi-3D spreadsheet
    ri[i]=bdata->m_pos.value(i);
    ri_1[i]=bdata->m_pos.value(i+1);
    ri_2[i]=bdata->m_pos.value(i+2);
    B=bdata->m_pos.value(number_of_segments-1)/2.;
    A[i]=ri_1[i]+(ri_2[i]-ri_1[i])/2.;
    re[i]=sqrt(pow((A[i]-B),2)+pow(m_parameter.distanceObsever,2));
    phi_lin[i]=qRadiansToDegrees(qAsin(m_parameter.distanceObsever/re[i]));
    phi[i]=180.-phi_lin[i];
    theta[i]=90;
    dist_obs[i]=ri_1[i];
}
else if (m_parameter.phi_type==1){
    //    re phi and theta calculation p 77 C_Project_Log_Text_Jan_16.pdf
    //    Input X e , Y e , Z e
    //    Attribute their respective values to X B , Y B , Z B
        XB=m_parameter.obs_x_pos;
        YB=m_parameter.obs_y_pos;
        ZB=m_parameter.obs_z_pos;

    //For the blade, get the Pitch angle, θp - blade pitch

    //For each blade segment, get: C r i+1 , C r i , r i+1 , r i , β
    //Cri+1 = c_1
    //r+1 = r_1
    //Cri = c_0
    //ri = r_0
    //beta = local_twist

    if((i<=number_of_segments) & (i>=number_of_segments)){
c_1[i]=0;
r_1[i]=0;
    }
    else
    {
c_1[i]=bdata->m_c_local.value(i+1);
r_1[i]=bdata->m_pos.value(i+1);
    }}

c_0[i]=bdata->m_c_local.value(i);
r_0[i]=bdata->m_pos.value(i);

local_twist[i]=theta_BEM[i];

    b[i]=qRadiansToDegrees(qAtan((c_1[i]-c_0[i])/(r_1[i]-r_0[i])));

//    the angle a is the total angle between the Y B Z B blade reference system plane and the local midsection chord line
    a[i]=local_twist[i]+blade_pitch;

    XRS[i]=XB*cos(qDegreesToRadians(a[i]))+YB*sin(qDegreesToRadians(a[i]));
    YRS[i]=-XB*sin(qDegreesToRadians(a[i]))+YB*cos(qDegreesToRadians(a[i]));
    ZRS[i]=ZB-(bdata->m_pos.value(i)-bdata->m_pos.value(i-1))/2.;

//    Calculate Y RS − 0.75 ∗ (C r i+1 − C r i )/2
    calc_int_a[i]=(YRS[i]-0.75*(c_1[i]-c_0[i])/2.);

    XRT[i]=XRS[i];
    YRT[i]=cos(qDegreesToRadians(b[i]))*calc_int_a[i]+sin(qDegreesToRadians(b[i]))*ZRS[i];
    ZRT[i]=-sin(qDegreesToRadians(b[i]))*calc_int_a[i]+cos(qDegreesToRadians(b[i]))*ZRS[i];

    r_e[i]=sqrt(pow(XRT[i],2)+pow(YRT[i],2)+pow(ZRT[i],2));
    r_rt[i]=r_e[i];
    theta_e[i]=qRadiansToDegrees(qAtan(ZRT[i]/YRT[i]));

    psi_e[i]=qRadiansToDegrees(qAtan(XRT[i]/ZRT[i]));

    dist_obs[i]=re[i];

   phi_rad[i]=qDegreesToRadians(phi[i]);
   theta_rad[i]=qDegreesToRadians(theta[i]);

   alpha_error[i]=qFabs(alpha_polar[i]-alpha_BEM[i])/alpha_BEM[i]*100.;

first_term_Dh_S[i]=10.*log10(pow(Mach[i],5)*L[i]*Dh[i]*D_starred_S[i]/pow(dist_obs[i],2));

first_term_Dl_S[i]=10.*log10(pow(Mach[i],5)*L[i]*Dl[i]*D_starred_S[i]/pow(dist_obs[i],2));

first_term_Dh_P[i]=10.*log10(pow(Mach[i],5)*L[i]*Dh[i]*D_starred_P[i]/pow(dist_obs[i],2));

St1[i]=0.02*(pow(Mach[i],-0.6));

if(alpha[i]<1.33){St2[i]=St1[i];}
else if(alpha[i]>12.5){St2[i]=St1[i]*4.72;}
else {St2[i]=St1[i]*pow(10.,(0.0054*pow((alpha[i]-1.33),2)));}

if (bdata->m_Reynolds.value(i)<95200)
{b0[i]= 0.3;}
else if (bdata->m_Reynolds.value(i)>857000)
{b0[i]= 0.56;}
else {b0[i]=-4.48*pow(10,-13)*(pow((bdata->m_Reynolds.value(i)-857000.),2)+0.56);}

if (b0[i]<0.13)
{B_min_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
else if (b0[i]>0.145)
{B_min_b0[i]=-817.81*pow(b0[i],3)+335.21*pow(b0[i],2)-135.024*b0[i]+10.619;}
else{B_min_b0[i]=-83.607*b0[i]+8.138;}

if (b0[i]<0.1)
{B_max_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
else if (b0[i]>0.187)
{B_max_b0[i]=-80.541*pow(b0[i],3)+44.174*pow(b0[i],2)-39.381*b0[i]+2.344;}
else {B_max_b0[i]=-31.33*b0[i]+1.854;}

BR_b0[i]=(-20-B_min_b0[i])/(B_max_b0[i]-B_min_b0[i]);

RCmod[i]=3*Reynolds[i];

if (RCmod[i]<95200){ao_Rc[i]=0.57;}
else if (RCmod[i]>857000){ao_Rc[i]=1.13;}
else {ao_Rc[i]=-9.57*pow(10,-13)*(pow((RCmod[i]-857000),2)+1.13);}

if(ao_Rc[i]<0.204){A_min_ao[i]=sqrt(67.552-886.788*pow(ao_Rc[i],2))-8.219;}
else if (ao_Rc[i]>0.244){A_min_ao[i]=-142.795*pow(ao_Rc[i],3)+103.656*pow(ao_Rc[i],2)-57.757*ao_Rc[i]+6.006;}
else {A_min_ao[i]=-32.665*ao_Rc[i]+3.981;}

if (ao_Rc[i]<0.13){A_max_ao[i]=sqrt(67.552-886.788*pow(ao_Rc[i],2))-8.219;}
else if(ao_Rc[i]>0.321){A_max_ao[i]=-4.669*pow(ao_Rc[i],3)+3.491*pow(ao_Rc[i],2)-16.699*ao_Rc[i]+1.149;}
else {A_max_ao[i]=-15.901*ao_Rc[i]+1.098;}

K1_3[i]=K1[i]-3.;

AR_ao[i]=(-20-A_min_ao[i])/(A_max_ao[i]-A_min_ao[i]);

St1_bar[i]=(St1[i]+St2[i])/2.;

Re_disp_thick[i]=1.225*bdata->m_Windspeed.value(i)*D_starred_P[i]/(0.0000178);

if (Re_disp_thick[i]>5000){delta_K1[i]=0;}
else {delta_K1[i]=alpha[i]*(1.43*log10(Re_disp_thick[i])-5.29);}

    double b0[number_of_segments];
    if (bdata->m_Reynolds.value(i)<95200)
    {b0[i]= 0.3;}
    else if (bdata->m_Reynolds.value(i)>857000)
    {b0[i]= 0.56;}
    else {b0[i]=-4.48*pow(10.,-13)*(pow((bdata->m_Reynolds.value(i)-857000.),2)+0.56);}

    double B_min_b0[number_of_segments];
    if (b0[i]<0.13)
    {B_min_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
    else if (b0[i]>0.145)
    {B_min_b0[i]=-817.81*pow(b0[i],3)+335.21*pow(b0[i],2)-135.024*b0[i]+10.619;}
    else{B_min_b0[i]=-83.607*b0[i]+8.138;}

    double B_max_b0[number_of_segments];
    if (b0[i]<0.1)
    {B_max_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
    else if (b0[i]>0.187)
    {B_max_b0[i]=-80.541*pow(b0[i],3)+44.174*pow(b0[i],2)-39.381*b0[i]+2.344;}
    else {B_max_b0[i]=-31.33*b0[i]+1.854;}

    double BR_b0[number_of_segments];
    BR_b0[i]=(-20-B_min_b0[i])/(B_max_b0[i]-B_min_b0[i]);

//    if((z<=m_parameter.TSRtd) & (z>=m_parameter.TSRtd)){
    QString str= QString::number(z, 'f', 1);


//        stream << "Angles less than the switching angle: "  << endl;
        stream << "Tip Speed Ratio: " << str << endl;
        stream << "Section: " << (i+1)<<"/"<<number_of_segments << endl;
        stream << "Reynolds: " << Reynolds[i] << endl;
        stream << "alpha " << alpha[i] << endl;
        stream << endl;

if(m_parameter.Lowson_type!=0){
    stream << qSetFieldWidth(14)  <<
              "Freq[Hz]"  << ";" <<
//              "Sts"  << ";" <<
//              "b"  << ";" <<
//              "B_min(b)" << ";" <<
//              "B_max(b)" << ";" <<
//              "B(b)"  << ";" <<
//              "a_alpha"  << ";" <<
//              "A_min(a)_alpha" << ";" <<
//              "A_max(a)_alpha" << ";" <<
//              "A'(a)"  << ";" <<
//              "A-Weighting"  << ";" <<
//              "B-Weighting" << ";" <<
//              "C-Weighting"  << ";" <<
//              "SPL_alpha < 0"  << ";" <<
//              "SPL_alpha > 0"  << ";" <<
//              "dB_A_alpha S"  << ";" <<
//              "dB_A_alpha P"  << ";" <<
//              "dB_B_alpha S"  << ";" <<
//              "dB_B_alpha P"  << ";" <<
//              "dB_C_alpha S"  << ";" <<
//              "dB_C_alpha P"  << ";" <<
//              "St1_bar"  << ";" <<
//              "a_S"  << ";" <<
//              "A_min_a_S"  << ";" <<
//              "A_max_a_S"  << ";" <<
//              "A_a_S"  << ";" <<
//              "SPL_dB_S"  << ";" <<
//              "dB_A_S"  << ";" <<
//              "dB_B_S"  << ";" <<
//              "dB_C_S"  << ";" <<
//              "Sts_St1_bar"  << ";" <<
//              "Stp"  << ";" <<
//              "a_P"  << ";" <<
//              "A_min_a_P"  << ";" <<
//              "A_max_a_P"  << ";" <<
//              "A_a_P"  << ";" <<
//              "SPL_dB_P"  << ";" <<
//              "dB_A_P"  << ";" <<
//              "dB_B_P"  << ";" <<
//              "dB_C_P"  << ";" <<
              "SPL (dB)"  << ";" <<
              "SPLa"  << ";" <<
              "SPLs"  << ";" <<
              "SPLp"  << ";" <<
              "SPL (dB(A))"  << ";" <<
              "SPL (dB(B))"  << ";" <<
              "SPL (dB(C))"  << ";" <<
              "SPL_LE (dB) "  << ";" <<
              "SPL_LE ((dB(A)) "  << ";" <<
              "SPL_LE ((dB(B)) "  << ";" <<
              "SPL_LE (dB(C)) "  << ";" <<
//              "s_log_SPL_alpha"  << ";" <<
//              "s_log_SPL_S"  << ";" <<
//              "s_log_SPL_P"  << ";" <<
//              "s_log_SPL"  << ";" <<
//              "s_log_dBA"  << ";" <<
//              "s_log_dBB"  << ";" <<
//              "s_log_dBC"  << ";" <<
//              "s_log_LE_dB"  << ";" <<
              endl;}
else{
        stream << qSetFieldWidth(14)  <<
                  "Freq[Hz]"  << ";" <<
//                  "Sts"  << ";" <<
//                  "b"  << ";" <<
//                  "B_min(b)" << ";" <<
//                  "B_max(b)" << ";" <<
//                  "B(b)"  << ";" <<
//                  "a_alpha"  << ";" <<
//                  "A_min(a)_alpha" << ";" <<
//                  "A_max(a)_alpha" << ";" <<
//                  "A'(a)"  << ";" <<
//                  "A-Weighting"  << ";" <<
//                  "B-Weighting" << ";" <<
//                  "C-Weighting"  << ";" <<
//                  "SPL_alpha < 0"  << ";" <<
//                  "SPL_alpha > 0"  << ";" <<
//                  "dB_A_alpha S"  << ";" <<
//                  "dB_A_alpha P"  << ";" <<
//                  "dB_B_alpha S"  << ";" <<
//                  "dB_B_alpha P"  << ";" <<
//                  "dB_C_alpha S"  << ";" <<
//                  "dB_C_alpha P"  << ";" <<
//                  "St1_bar"  << ";" <<
//                  "a_S"  << ";" <<
//                  "A_min_a_S"  << ";" <<
//                  "A_max_a_S"  << ";" <<
//                  "A_a_S"  << ";" <<
//                  "SPL_dB_S"  << ";" <<
//                  "dB_A_S"  << ";" <<
//                  "dB_B_S"  << ";" <<
//                  "dB_C_S"  << ";" <<
//                  "Sts_St1_bar"  << ";" <<
//                  "Stp"  << ";" <<
//                  "a_P"  << ";" <<
//                  "A_min_a_P"  << ";" <<
//                  "A_max_a_P"  << ";" <<
//                  "A_a_P"  << ";" <<
//                  "SPL_dB_P"  << ";" <<
//                  "dB_A_P"  << ";" <<
//                  "dB_B_P"  << ";" <<
//                  "dB_C_P"  << ";" <<
                  "SPL (dB)"  << ";" <<
                  "SPLa"  << ";" <<
                  "SPLs"  << ";" <<
                  "SPLp"  << ";" <<
                  "SPL (dB(A))"  << ";" <<
                  "SPL (dB(B))"  << ";" <<
                  "SPL (dB(C))"  << ";" <<
//                  "s_log_SPL_alpha"  << ";" <<
//                  "s_log_SPL_S"  << ";" <<
//                  "s_log_SPL_P"  << ";" <<
//                  "s_log_SPL"  << ";" <<
//                  "s_log_dBA"  << ";" <<
//                  "s_log_dBB"  << ";" <<
//                  "s_log_dBC"  << ";" <<
                  endl;}

    int w=30;

    double slog_SPL_alpha[w];
    double slog_SPL_S[w];
    double slog_SPL_P[w];
    double slog_SPL[w];
    double slog_dBA[w];
    double slog_dBB[w];
    double slog_dBC[w];
    double slog_LedB[w];
    double slog_LedBAW[w];
    double slog_LedBBW[w];
    double slog_LedBCW[w];

     for (int j = 0; j < w; ++j) {

         double Sts[w];
         Sts[j]=Frequency[j]*D_starred_S[i]/bdata->m_Windspeed.value(i);

         double b_alpha[w];
         b_alpha[j]=qFabs(log10(Sts[j]/St2[i]));

         double B_min[w];
    if (b_alpha[j]<0.13)
    {B_min[j]=sqrt(16.888-886.788*pow(b_alpha[j],2));}
    else if(b_alpha[j]>0.145)
    {B_min[j]=-817.81*pow(b_alpha[j],3)+335.21*pow(b_alpha[j],2)-135.024*b_alpha[j]+10.619;}
    else {B_min[j]=-83.607*b_alpha[j]+8.138;}

    double B_max[w];
    if (b_alpha[j]<0.1)
    {B_max[j]=sqrt(16.888-886.788*pow(b_alpha[j],2))-4,109;}
    else if(b_alpha[j]>0.187)
    {B_max[j]=-80.541*pow(b_alpha[j],3)+44.174*pow(b_alpha[j],2)-39.381*b_alpha[j]+2.344;}
    else {B_max[j]=-31.33*b_alpha[j]+1.854;}

    double B_b[w];
    B_b[j]=B_min[j]+BR_b0[i]*(B_max[j]-B_min[j]);
    if(qIsInf(B_b[j])){B_b[j]=0;}

    double a_alpha[w];
    a_alpha[j]=qFabs(log10(Sts[j]/St2[i]));

    double A_min_alpha[w];
    if (a_alpha[j]<0.204)
    {A_min_alpha[j]=sqrt(67.552-886.788*pow(a_alpha[j],2))-8.219;}
    else if(a_alpha[j]>0.244)
    {A_min_alpha[j]=-142.795*pow(a_alpha[j],3)+103.656*pow(a_alpha[j],2)-57.757*a_alpha[j]+6.006;}
    else {A_min_alpha[j]=-32.665*a_alpha[j]+3.981;}

    double A_max_alpha[w];
    if (a_alpha[j]<0.13)
    {A_max_alpha[j]=sqrt(67.552-886.788*pow(a_alpha[j],2))-8.219;}
    else if(a_alpha[j]>0.321)
    {A_max_alpha[j]=-4.669*pow(a_alpha[j],3)+3.491*pow(a_alpha[j],2)-16.699*a_alpha[j]+1.149;}
    else {A_max_alpha[j]=-15.901*a_alpha[j]+1.098;}

    double Alin_a[w];
    Alin_a[j]=A_min_alpha[j]+AR_ao[i]*(A_max_alpha[j]-A_min_alpha[j]);

    double SPL_alpha_min0[w];
    SPL_alpha_min0[j]=first_term_Dh_S[i]+K2[i]+B_b[j];

    double SPL_alpha_big0[w];
    SPL_alpha_big0[j]=first_term_Dl_S[i]+K2[i]+Alin_a[j];

    double dBA_alpha_min0[w];
    dBA_alpha_min0[j]=SPL_alpha_min0[j]+AWeighting[j];

    double dBA_alpha_big0[w];
    dBA_alpha_big0[j]=SPL_alpha_big0[j]+AWeighting[j];

    double dBB_alpha_min0[w];
    dBB_alpha_min0[j]=SPL_alpha_min0[j]+BWeighting[j];

    double dBC_alpha_big0[w];
    dBC_alpha_big0[j]=SPL_alpha_big0[j]+CWeighting[j];

    double dBC_alpha_min0[w];
    dBC_alpha_min0[j]=SPL_alpha_min0[j]+CWeighting[j];

    double dBB_alpha_big0[w];
    dBB_alpha_big0[j]=SPL_alpha_big0[j]+BWeighting[j];

    double St1_bar[w];
    St1_bar[j]=(St1[i]+St2[i])/2.;

    double a_S[w];
    a_S[j]=qFabs(log10(Sts[j]/St1_bar[j]));

    double A_min_S[w];
    if (a_S[j]<0.204){A_min_S[j]=sqrt(67.552-886.788*pow(a_S[j],2))-8.219;}
    else if (a_S[j]>0.244){A_min_S[j]=-142.795*pow(a_S[j],3)+103.656*pow(a_S[j],2)-57.757*a_S[j]+6.006;}
    else {A_min_S[j]=-32.665*a_S[j]+3.981;}

    double A_max_S[w];
    if (a_S[j]<0.13){A_max_S[j]=sqrt(67.552-886.788*pow(a_S[j],2))-8.219;}
    else if (a_S[j]>0.321){A_max_S[j]=-4.669*pow(a_S[j],3)+3.491*pow(a_S[j],2)-16.699*a_S[j]+1.149;}
    else {A_max_S[j]=-15.901*a_S[j]+1.098;}

    double A_a_S[w];
    A_a_S[j]=A_min_S[j]+AR_ao[i]*(A_max_S[j]-A_min_S[j]);

    double SPL_dB_S[w];
    SPL_dB_S[j]=first_term_Dh_S[i]+A_a_S[j]+K1_3[i];

    double dBA_S[w];
    dBA_S[j]=SPL_dB_S[j]+AWeighting[j];

    double dBB_S[w];
    dBB_S[j]=SPL_dB_S[j]+BWeighting[j];

    double dBC_S[w];
    dBC_S[j]=SPL_dB_S[j]+CWeighting[j];

    double Sts_St1_bar[w];
    Sts_St1_bar[j]=Sts[j]/St1_bar[i];

    double Stp_P[w];
    Stp_P[j]=Frequency[j]*D_starred_P[i]/bdata->m_Windspeed.value(i);

    double a_P[w];
    a_P[j]=qFabs(log10(Stp_P[j]/St1[i]));

    double A_min_P[w];
    if(a_P[j]<0.204){A_min_P[j]=sqrt(67.552-886.788*pow(a_P[j],2))-8.219;}
    else if (a_P[j]>0.244){A_min_P[j]=-142.795*pow(a_P[j],3)+103.656*pow(a_P[j],2)-57.757*a_P[j]+6.006;}
    else {A_min_P[j]=-32.665*a_P[j]+3.981;}

    double A_max_P[w];
    if(a_P[j]<0.13){A_max_P[j]=sqrt(67.552-886.788*pow(a_P[j],2))-8.219;}
    else if (a_P[j]>0.321){A_max_P[j]=-4.669*pow(a_P[j],3)+3.491*pow(a_P[j],2)-16.699*a_P[j]+1.149;}
    else {A_max_P[j]=-15.901*a_P[j]+1.098;}

    double A_a_P[w];
    A_a_P[j]=A_min_P[j]+AR_ao[i]*(A_max_P[j]-A_min_P[j]);

    double SPL_dB_P[w];
    SPL_dB_P[j]=delta_K1[i]+A_a_P[j]+K1_3[i]+first_term_Dh_P[i];

    double dBA_P[w];
    dBA_P[j]=SPL_dB_P[j]+AWeighting[j];

    double dBB_P[w];
    dBB_P[j]=SPL_dB_P[j]+BWeighting[j];

    double dBC_P[w];
    dBC_P[j]=SPL_dB_P[j]+CWeighting[j];

    double SPL_alpha[w];
    if (alpha[i]<SwAlpha[i]){SPL_alpha[j]=SPL_alpha_min0[j];}
    else {SPL_alpha[j]=SPL_alpha_big0[j];}

    double SPL_S[w];
    if (alpha[i]<SwAlpha[i]){SPL_S[j]=SPL_dB_S[j];}
    else {SPL_S[j]=-999999999999.;}

    double SPL_P[w];
    if (alpha[i]<SwAlpha[i]){SPL_P[j]=SPL_dB_P[j];}
    else {SPL_P[j]=-999999999999.;}

    double SPL_dB[w];
    SPL_dB[j]=10*log10(pow(10.,(SPL_alpha[j]/10.))+pow(10.,(SPL_S[j]/10.))+pow(10.,(SPL_P[j]/10.)));

    double SPL_A[w];
    SPL_A[j]=SPL_dB[j]+AWeighting[j];

    double SPL_B[w];
    SPL_B[j]=SPL_dB[j]+BWeighting[j];

    double SPL_C[w];
    SPL_C[j]=SPL_dB[j]+CWeighting[j];

    double SPL_LedB[w];
    SPL_LedB[w]=0;

    double SPL_LedBAW[w];
    SPL_LedBAW[w]=0;

    double SPL_LedBBW[w];
    SPL_LedBBW[w]=0;

    double SPL_LedBCW[w];
    SPL_LedBCW[w]=0;

    u_le=m_parameter.u_wind_speed*100.;
    c_le=100.*chord[i];
    I_le=m_parameter.TurbulenceIntensity;
    lambda_le=100.*m_parameter.IntegralLengthScale;
    r_e_le=100.*dist_obs[i];
    beta_le=sqrt(1-pow(Mach[i],2));
    L_le=100.*L[i];
    D_L_le=0.5*Dl[i];
    aux_le=0.5*(lambda_le*L_le*pow(rho, 2)*pow(c_0_le, 2)*pow(u_le, 2)*pow(Mach[i], 3)*pow(I_le, 2)*D_L_le)/(pow(r_e_le, 2));
    K_le=M_PI*Frequency[i]*c_le/u_le;
    S_le=sqrt(pow((2.*M_PI*K_le/(pow(beta_le, 2)))+(pow((1+(2.4*K_le/pow(beta_le,2))),-1)),-1));
    LFC_le = 10.*Mach[i]*pow(S_le*K_le/beta_le,2);

    if(m_parameter.Lowson_type==2){
        c_const_le = c_const_rd_le;
        d_const_le = d_const_rd_le;
     }
    else if(m_parameter.Lowson_type==1){
            c_const_le = c_const_vk_le;
            d_const_le = d_const_vk_le;
         }
        else {
            c_const_le=0;
            d_const_le = 0;
        }

aux0_le=0.5*(lambda_le*L_le*pow(rho, 2)*pow(c_0_le, 2)*pow(u_le, 2)*pow(Mach[i], 3)*pow(I_le, 2)*Dl[i])/(pow(r_e_le, 2));
aux1_le=10.*log10(pow(LFC_le/(1+LFC_le), 2))+d_const_le;
aux4_le=pow(K_le,3)/pow(1+(pow(K_le,2)),c_const_le);
aux5_le=10.*log10(aux0_le*aux4_le);

//Validation:

//Lowson validation:
if (((((m_parameter.Lowson_type!=0) & (Mach[i]<=0.18)) & (Mach[i]>0) & (Reynolds[i]<=(6.*pow(10,5)))) & (Reynolds[i]>0))){
SPL_LedB[j]=10.*log10(pow(10,(aux1_le+aux5_le)/10.));
SPL_LedBAW[j]=SPL_LedB[j]+AWeighting[j];
SPL_LedBBW[j]=SPL_LedB[j]+BWeighting[j];
SPL_LedBCW[j]=SPL_LedB[j]+CWeighting[j];
LE_validation=true;
}
else{
SPL_LedB[j]=-999999999999.;
SPL_LedBAW[j]=-999999999999.;
SPL_LedBBW[j]=-999999999999.;
SPL_LedBCW[j]=-999999999999.;
LE_validation=false;
}

//BPM validation:
//p 17 C_Project_Log_Text_15_jan_16
if(((((alpha[i]<=19.8 & Reynolds[i]<3*pow(10,6)) & (Mach[i]<0.21)) & (Reynolds[i]>0)) & (Mach[i]>0))){
BPM_validation=true;
}
else{
SPL_alpha[j]=-999999999999.;
SPL_S[j]=-999999999999.;
SPL_P[j]=-999999999999.;
SPL_dB[j]=-999999999999.;
SPL_A[j]=-999999999999.;
SPL_B[j]=-999999999999.;
SPL_C[j]=-999999999999.;
BPM_validation=false;
}

    slog_SPL_alpha[j]=pow(10.,(SPL_alpha[j]/10.));
    slog_SPL_S[j]=pow(10.,(SPL_S[j]/10.));
    slog_SPL_P[j]=pow(10.,(SPL_P[j]/10.));
    slog_SPL[j]=pow(10.,(SPL_dB[j]/10.));
    slog_dBA[j]=pow(10.,(SPL_A[j]/10.));
    slog_dBB[j]=pow(10.,(SPL_B[j]/10.));
    slog_dBC[j]=pow(10.,(SPL_C[j]/10.));
    slog_LedB[j]=pow(10.,(SPL_LedB[j]/10.));
    slog_LedBAW[j]=pow(10.,(SPL_LedBAW[j]/10.));
    slog_LedBBW[j]=pow(10.,(SPL_LedBBW[j]/10.));
    slog_LedBCW[j]=pow(10.,(SPL_LedBCW[j]/10.));

    if(qIsNaN(slog_SPL_alpha[j])){slog_SPL_alpha[j]=0;}
    if(qIsNaN(slog_SPL_S[j])){slog_SPL_S[j]=0;}
    if(qIsNaN(slog_SPL_P[j])){slog_SPL_P[j]=0;}
    if(qIsNaN(slog_SPL[j])){slog_SPL[j]=0;}
    if(qIsNaN(slog_dBA[j])){slog_dBA[j]=0;}
    if(qIsNaN(slog_dBB[j])){slog_dBB[j]=0;}
    if(qIsNaN(slog_dBC[j])){slog_dBC[j]=0;}
    if(qIsNaN(slog_LedB[j])){slog_LedB[j]=0;}
    if(qIsNaN(slog_LedBAW[j])){slog_LedBAW[j]=0;}
    if(qIsNaN(slog_LedBBW[j])){slog_LedBBW[j]=0;}
    if(qIsNaN(slog_LedBCW[j])){slog_LedBCW[j]=0;}

    splog_OASPL_alpha+=slog_SPL_alpha[j];
    splog_OASPL_S+=slog_SPL_S[j];
    splog_OASPL_P+=slog_SPL_P[j];
    splog_OASPL+=slog_SPL[j];
    splog_dBA+=slog_dBA[j];
    splog_dBB+=slog_dBB[j];
    splog_dBC+=slog_dBC[j];
    splog_LedB+=slog_LedB[j];
    splog_LedBAW+=slog_LedBAW[j];
    splog_LedBBW+=slog_LedBBW[j];
    splog_LedBCW+=slog_LedBCW[j];

    splog_OASPL_alpha_a+=slog_SPL_alpha[j];
    splog_OASPL_S_a+=slog_SPL_S[j];
    splog_OASPL_P_a+=slog_SPL_P[j];
    splog_OASPL_a+=slog_SPL[j];
    splog_dBA_a+=slog_dBA[j];
    splog_dBB_a+=slog_dBB[j];
    splog_dBC_a+=slog_dBC[j];
    splog_LedB_a+=slog_LedB[j];
    splog_LedBAW_a+=slog_LedBAW[j];
    splog_LedBBW_a+=slog_LedBBW[j];
    splog_LedBCW_a+=slog_LedBCW[j];

    sp_OASPL_alpha=10*log10(splog_OASPL_alpha);
    sp_OASPL_S=10*log10(splog_OASPL_S);
    sp_OASPL_P=10*log10(splog_OASPL_P);
    sp_OASPL=10*log10(splog_OASPL);
    sp_dBA=10*log10(splog_dBA);
    sp_dBB=10*log10(splog_dBB);
    sp_dBC=10*log10(splog_dBC);
    sp_LedB=10*log10(splog_LedB);
    sp_LedBAW=10*log10(splog_LedBAW);
    sp_LedBBW=10*log10(splog_LedBBW);
    sp_LedBCW=10*log10(splog_LedBCW);

    sp_OASPL_alpha_a=10*log10(splog_OASPL_alpha_a);
    sp_OASPL_S_a=10*log10(splog_OASPL_S_a);
    sp_OASPL_P_a=10*log10(splog_OASPL_P_a);
    sp_OASPL_a=10*log10(splog_OASPL_a);
    sp_dBA_a=10*log10(splog_dBA_a);
    sp_dBB_a=10*log10(splog_dBB_a);
    sp_dBC_a=10*log10(splog_dBC_a);
    sp_LedB_a=10*log10(splog_LedB_a);
    sp_LedBAW_a=10*log10(splog_LedBAW_a);
    sp_LedBBW_a=10*log10(splog_LedBBW_a);
    sp_LedBCW_a=10*log10(splog_LedBCW_a);

    if(qIsInf(sp_OASPL_alpha)){sp_OASPL_alpha=0;}
    if(qIsInf(sp_OASPL_S)){sp_OASPL_S=0;}
    if(qIsInf(sp_OASPL_P)){sp_OASPL_P=0;}
    if(qIsInf(sp_OASPL)){sp_OASPL=0;}
    if(qIsInf(sp_dBA)){sp_dBA=0;}
    if(qIsInf(sp_dBB)){sp_dBB=0;}
    if(qIsInf(sp_dBC)){sp_dBC=0;}
    if(qIsInf(sp_LedB)){sp_LedB=0;}
    if(qIsInf(sp_LedBAW)){sp_LedBAW=0;}
    if(qIsInf(sp_LedBBW)){sp_LedBBW=0;}
    if(qIsInf(sp_LedBCW)){sp_LedBCW=0;}

    if(qIsInf(sp_OASPL_alpha)){sp_OASPL_alpha=0;}
    if(qIsInf(splog_OASPL_S)){splog_OASPL_S=0;}
    if(qIsInf(splog_OASPL_P)){splog_OASPL_P=0;}
    if(qIsInf(splog_OASPL)){splog_OASPL=0;}
    if(qIsInf(splog_dBA)){splog_dBA=0;}
    if(qIsInf(splog_dBB)){splog_dBB=0;}
    if(qIsInf(splog_dBC)){splog_dBC=0;}
    if(qIsInf(splog_LedB)){splog_LedB=0;}
    if(qIsInf(splog_LedBAW)){splog_LedBAW=0;}
    if(qIsInf(splog_LedBBW)){splog_LedBBW=0;}
    if(qIsInf(splog_LedBCW)){splog_LedBCW=0;}

    if(qIsInf(sp_OASPL_alpha_a)){sp_OASPL_alpha_a=0;}
    if(qIsInf(sp_OASPL_S_a)){sp_OASPL_S_a=0;}
    if(qIsInf(sp_OASPL_P_a)){sp_OASPL_P_a=0;}
    if(qIsInf(sp_OASPL_a)){sp_OASPL_a=0;}
    if(qIsInf(sp_dBA_a)){sp_dBA_a=0;}
    if(qIsInf(sp_dBB_a)){sp_dBB_a=0;}
    if(qIsInf(sp_dBC_a)){sp_dBC_a=0;}
    if(qIsInf(sp_LedB_a)){sp_LedB_a=0;}
    if(qIsInf(sp_LedBAW_a)){sp_LedBAW_a=0;}
    if(qIsInf(sp_LedBBW_a)){sp_LedBBW_a=0;}
    if(qIsInf(sp_LedBCW_a)){sp_LedBCW_a=0;}

    if(qIsInf(sp_OASPL_alpha_a)){sp_OASPL_alpha_a=0;}
    if(qIsInf(splog_OASPL_S_a)){splog_OASPL_S_a=0;}
    if(qIsInf(splog_OASPL_P_a)){splog_OASPL_P_a=0;}
    if(qIsInf(splog_OASPL_a)){splog_OASPL_a=0;}
    if(qIsInf(splog_dBA_a)){splog_dBA_a=0;}
    if(qIsInf(splog_dBB_a)){splog_dBB_a=0;}
    if(qIsInf(splog_dBC_a)){splog_dBC_a=0;}
    if(qIsInf(splog_LedB_a)){splog_LedB_a=0;}
    if(qIsInf(splog_LedBAW_a)){splog_LedBAW_a=0;}
    if(qIsInf(splog_LedBBW_a)){splog_LedBBW_a=0;}
    if(qIsInf(splog_LedBCW_a)){splog_LedBCW_a=0;}

    stlog_OASPL_alpha_a+=splog_OASPL_alpha_a;
    stlog_OASPL_S_a+=splog_OASPL_S_a;
    stlog_OASPL_P_a+=splog_OASPL_P_a;
    stlog_OASPL_a+=splog_OASPL_a;
    stlog_dBA_a+=splog_dBA_a;
    stlog_dBB_a+=splog_dBB_a;
    stlog_dBC_a+=splog_dBC_a;
    stlog_LedB_a+=splog_LedB_a;
    stlog_LedBAW_a+=splog_LedBAW_a;
    stlog_LedBBW_a+=splog_LedBBW_a;
    stlog_LedBCW_a+=splog_LedBCW_a;

    st_OASPL_alpha_a=10*log10(stlog_OASPL_alpha_a);
    st_OASPL_S_a=10*log10(stlog_OASPL_S_a);
    st_OASPL_P_a=10*log10(stlog_OASPL_P_a);
    st_OASPL_a=10*log10(stlog_OASPL_a);
    st_dBA_a=10*log10(stlog_dBA_a);
    st_dBB_a=10*log10(stlog_dBB_a);
    st_dBC_a=10*log10(stlog_dBC_a);
    st_LedB_a=10*log10(stlog_LedB_a);
    st_LedBAW_a=10*log10(stlog_LedBAW_a);
    st_LedBBW_a=10*log10(stlog_LedBBW_a);
    st_LedBCW_a=10*log10(stlog_LedBAW_a);

    if(qIsInf(st_OASPL_alpha_a)){st_OASPL_alpha_a=0;}
    if(qIsInf(st_OASPL_S_a)){st_OASPL_S_a=0;}
    if(qIsInf(st_OASPL_P_a)){st_OASPL_P_a=0;}
    if(qIsInf(st_OASPL_a)){st_OASPL_a=0;}
    if(qIsInf(st_dBA_a)){st_dBA_a=0;}
    if(qIsInf(st_dBB_a)){st_dBB_a=0;}
    if(qIsInf(st_dBC_a)){st_dBC_a=0;}
    if(qIsInf(st_LedB_a)){st_LedB_a=0;}
    if(qIsInf(st_LedBAW_a)){st_LedBAW_a=0;}
    if(qIsInf(st_LedBBW_a)){st_LedBBW_a=0;}
    if(qIsInf(st_LedBCW_a)){st_LedBCW_a=0;}

    stlog_OASPL_alpha+=splog_OASPL_alpha;
    stlog_OASPL_S+=splog_OASPL_S;
    stlog_OASPL_P+=splog_OASPL_P;
    stlog_OASPL+=splog_OASPL;
    stlog_dBA+=splog_dBA;
    stlog_dBB+=splog_dBB;
    stlog_dBC+=splog_dBC;
    stlog_LedB+=splog_LedB;
    stlog_LedBAW+=splog_LedBAW;
    stlog_LedBBW+=splog_LedBBW;
    stlog_LedBCW+=splog_LedBCW;

    st_OASPL_alpha=10*log10(stlog_OASPL_alpha);
    st_OASPL_S=10*log10(stlog_OASPL_S);
    st_OASPL_P=10*log10(stlog_OASPL_P);
    st_OASPL=10*log10(stlog_OASPL);
    st_dBA=10*log10(stlog_dBA);
    st_dBB=10*log10(stlog_dBB);
    st_dBC=10*log10(stlog_dBC);
    st_LedB=10*log10(stlog_LedB);
    st_LedBAW=10*log10(stlog_LedBAW);
    st_LedBBW=10*log10(stlog_LedBBW);
    st_LedBCW=10*log10(stlog_LedBAW);

    if(qIsInf(st_OASPL_alpha)){st_OASPL_alpha=0;}
    if(qIsInf(st_OASPL_S)){st_OASPL_S=0;}
    if(qIsInf(st_OASPL_P)){st_OASPL_P=0;}
    if(qIsInf(st_OASPL)){st_OASPL=0;}
    if(qIsInf(st_dBA)){st_dBA=0;}
    if(qIsInf(st_dBB)){st_dBB=0;}
    if(qIsInf(st_dBC)){st_dBC=0;}
    if(qIsInf(st_LedB)){st_LedB=0;}
    if(qIsInf(st_LedBAW)){st_LedBAW=0;}
    if(qIsInf(st_LedBBW)){st_LedBBW=0;}
    if(qIsInf(st_LedBCW)){st_LedBCW=0;}

    if (m_parameter.Lowson_type==2)   {
SPL_blade_total_aux+=st_OASPL_alpha+st_OASPL_S+st_OASPL_P+st_OASPL+st_dBA+st_dBB+st_dBC+st_LedB+st_LedBAW+st_LedBBW+st_LedBCW;

SPL_blade_total_aux_a+=st_OASPL_alpha_a+st_OASPL_S_a+st_OASPL_P_a+st_OASPL_a+st_dBA_a+st_dBB_a+st_dBC_a+st_LedB_a+st_LedBAW_a+st_LedBBW_a+st_LedBCW_a;
    }
    else {
SPL_blade_total_aux=st_OASPL_alpha+st_OASPL_S+st_OASPL_P+st_OASPL+st_dBA+st_dBB+st_dBC+st_LedB+st_LedBAW+st_LedBBW;

SPL_blade_total_aux_a+=st_OASPL_alpha_a+st_OASPL_S_a+st_OASPL_P_a+st_OASPL_a+st_dBA_a+st_dBB_a+st_dBC_a+st_LedB_a+st_LedBAW_a+st_LedBBW_a;
            }

if (j==w){
SPL_blade_total = 10.*log10(1./number_of_segments*(SPL_blade_total_aux));

SPL_blade_total_a = 10.*log10(1./number_of_segments*(SPL_blade_total_aux_a));
}

    QString observations_x("");
    if (!((Reynolds[i] >9.5*pow(10,5)) & (Reynolds[i]<2.5*pow(10,6))) || Mach[i]>=0.19 || !((alpha[i]<=0) & (alpha[i]>=0)) || !((alpha[i]<=5) & (alpha[i]>=5)) || !((alpha[i]<=10) & (alpha[i]>=10))){
    observations_x.append("1 ");}
    if((Reynolds[i] <=4.8*pow(10,4)) || (Reynolds[i]>=2.5*pow(10,6)) || Mach[i]<0.208 || !((alpha[i]<=0) & (alpha[i]>=0))){
    observations_x.append("2 ");}
    if (Reynolds[i] >=3*pow(10,6) || Mach[i]<0.208 || alpha[i]<=19.8){
    observations_x.append("3 ");}
    if (Reynolds[i] >1.5*pow(10,6) || Mach[i]<0.208 || alpha[i]<=19.8){
    observations_x.append('4');}

    SPL_LedB_val.clear();
    SPL_LedB_valAW.clear();
    SPL_LedB_valBW.clear();
    SPL_LedB_valCW.clear();
    SPL_alpha_val.clear();
    SPL_dB_S_val.clear();
    SPL_dB_P_val.clear();
    SPL_S_val.clear();
    SPL_P_val.clear();
    SPL_dB_val.clear();
    SPL_A_val.clear();
    SPL_B_val.clear();
    SPL_C_val.clear();

        if(BPM_validation){
    SPL_alpha_val.append(QString::number(SPL_alpha[j], 'f', 2));
    SPL_dB_S_val.append(QString::number(SPL_dB_S[j], 'f', 2));
    SPL_dB_P_val.append(QString::number(SPL_dB_P[j], 'f', 2));
    SPL_dB_val.append(QString::number(SPL_dB_S[j], 'f', 2));
    SPL_A_val.append(QString::number(SPL_A[j], 'f', 2));
    SPL_B_val.append(QString::number(SPL_B[j], 'f', 2));
    SPL_C_val.append(QString::number(SPL_C[j], 'f', 2));

    if(SPL_S[j]==-999999999999.){SPL_S_val.append("N/A");} else {SPL_S_val.append(QString::number(SPL_S[j], 'f', 2));}
    if(SPL_P[j]==-999999999999.){SPL_P_val.append("N/A");} else {SPL_P_val.append(QString::number(SPL_P[j], 'f', 2));}
       }
        else{
    SPL_alpha_val.append("N/A");
    SPL_dB_S_val.append("N/A");
    SPL_dB_P_val.append("N/A");
    SPL_S_val.append("N/A");
    SPL_P_val.append("N/A");
    SPL_dB_val.append("N/A");
    SPL_A_val.append("N/A");
    SPL_B_val.append("N/A");
    SPL_C_val.append("N/A");
        }

    if (LE_validation){
    SPL_LedB_val.append(QString::number(SPL_LedB[j], 'f', 2));
    SPL_LedB_valAW.append(QString::number(SPL_LedBAW[j], 'f', 2));
    SPL_LedB_valBW.append(QString::number(SPL_LedBBW[j], 'f', 2));
    SPL_LedB_valCW.append(QString::number(SPL_LedBCW[j], 'f', 2));
    }
    else{
    SPL_LedB_val.append("N/A");
    SPL_LedB_valAW.append("N/A");
    SPL_LedB_valBW.append("N/A");
    SPL_LedB_valCW.append("N/A");
    }

if (m_parameter.Lowson_type!=0){
        stream << qSetFieldWidth(14)  <<
                  Frequency[j] << ";" <<
//                  Sts[j]<< ";" <<
//                  b_alpha[j]  << ";" <<
//                  B_min[j] << ";" <<
//                  B_max[j] << ";" <<
//                  B_b[j]  << ";" <<
//                  a_alpha[j] << ";" <<
//                  A_min_alpha[j]  << ";" <<
//                  A_max_alpha[j]  << ";" <<
//                  Alin_a[j]  << ";" <<
//                  AWeighting[j]  << ";" <<
//                  BWeighting[j]  << ";" <<
//                  CWeighting[j]  << ";" <<
//                  SPL_alpha_min0[j]  << ";" <<
//                  SPL_alpha_big0[j]  << ";" <<
//                  dBA_alpha_min0[j]  << ";" <<
//                  dBA_alpha_big0[j]  << ";" <<
//                  dBB_alpha_min0[j]  << ";" <<
//                  dBB_alpha_big0[j]  << ";" <<
//                  dBC_alpha_min0[j]  << ";" <<
//                  dBC_alpha_big0[j]  << ";" <<
//                  St1_bar[j]  << ";" <<
//                  a_S[j]  << ";" <<
//                  A_min_S[j]  << ";" <<
//                  A_max_S[j]  << ";" <<
//                  A_a_S[j]  << ";" <<
//                  SPL_dB_S_val  << ";" <<
//                  dBA_S[j]  << ";" <<
//                  dBB_S[j]  << ";" <<
//                  dBC_S[j]  << ";" <<
//                  Sts_St1_bar[j]  << ";" <<
//                  Stp_P[j]  << ";" <<
//                  a_P[j]  << ";" <<
//                  A_min_P[j]  << ";" <<
//                  A_max_P[j]  << ";" <<
//                  A_a_P[j]  << ";" <<
//                  SPL_dB_P_val  << ";" <<
//                  dBA_P[j]  << ";" <<
//                  dBB_P[j]  << ";" <<
//                  dBC_P[j]  << ";" <<
                  SPL_dB_val  << ";" <<
                  SPL_alpha_val  << ";" <<
                  SPL_S_val << ";" <<
                  SPL_P_val  << ";" <<
                  SPL_A_val  << ";" <<
                  SPL_B_val  << ";" <<
                  SPL_C_val  << ";" <<
                  SPL_LedB_val  << ";" <<
                  SPL_LedB_valAW  << ";" <<
                  SPL_LedB_valBW  << ";" <<
                  SPL_LedB_valCW  << ";" <<
//                  slog_SPL_alpha[j]  << ";" <<
//                  slog_SPL_S[j]  << ";" <<
//                  slog_SPL_P[j]  << ";" <<
//                  slog_SPL[j]  << ";" <<
//                  slog_dBA[j]  << ";" <<
//                  slog_dBB[j]  << ";" <<
//                  slog_dBC[j]  << ";" <<
//                  slog_LedB[j]  << ";" <<
                  endl;
}
else
{
        stream << qSetFieldWidth(14)  <<
                  Frequency[j] << ";" <<
//                  Sts[j]<< ";" <<
//                  b_alpha[j]  << ";" <<
//                  B_min[j] << ";" <<
//                  B_max[j] << ";" <<
//                  B_b[j]  << ";" <<
//                  a_alpha[j] << ";" <<
//                  A_min_alpha[j]  << ";" <<
//                  A_max_alpha[j]  << ";" <<
//                  Alin_a[j]  << ";" <<
//                  AWeighting[j]  << ";" <<
//                  BWeighting[j]  << ";" <<
//                  CWeighting[j]  << ";" <<
//                  SPL_alpha_min0[j]  << ";" <<
//                  SPL_alpha_big0[j]  << ";" <<
//                  dBA_alpha_min0[j]  << ";" <<
//                  dBA_alpha_big0[j]  << ";" <<
//                  dBB_alpha_min0[j]  << ";" <<
//                  dBB_alpha_big0[j]  << ";" <<
//                  dBC_alpha_min0[j]  << ";" <<
//                  dBC_alpha_big0[j]  << ";" <<
//                  St1_bar[j]  << ";" <<
//                  a_S[j]  << ";" <<
//                  A_min_S[j]  << ";" <<
//                  A_max_S[j]  << ";" <<
//                  A_a_S[j]  << ";" <<
//                  SPL_dB_S_val  << ";" <<
//                  dBA_S[j]  << ";" <<
//                  dBB_S[j]  << ";" <<
//                  dBC_S[j]  << ";" <<
//                  Sts_St1_bar[j]  << ";" <<
//                  Stp_P[j]  << ";" <<
//                  a_P[j]  << ";" <<
//                  A_min_P[j]  << ";" <<
//                  A_max_P[j]  << ";" <<
//                  A_a_P[j]  << ";" <<
//                  SPL_dB_P_val  << ";" <<
//                  dBA_P[j]  << ";" <<
//                  dBB_P[j]  << ";" <<
//                  dBC_P[j]  << ";" <<
                  SPL_dB_val  << ";" <<
                  SPL_alpha_val  << ";" <<
                  SPL_S_val << ";" <<
                  SPL_P_val  << ";" <<
                  SPL_A_val  << ";" <<
                  SPL_B_val  << ";" <<
                  SPL_C_val  << ";" <<
//                  slog_SPL_alpha[j]  << ";" <<
//                  slog_SPL_S[j]  << ";" <<
//                  slog_SPL_P[j]  << ";" <<
//                  slog_SPL[j]  << ";" <<
//                  slog_dBA[j]  << ";" <<
//                  slog_dBB[j]  << ";" <<
//                  slog_dBC[j]  << ";" <<
                  endl;
}

if(j==(w-1)){
        stream << endl;
        if (sp_OASPL<=0 & sp_OASPL>=0){stream << "OASPL: "  << "N/A"<< endl;} else {stream << "SPL: "  << sp_OASPL<< "dB" << endl;}
        if (sp_dBA<=0 & sp_dBA>=0){stream << "OASPL(A): "  << "N/A"<< endl;} else {stream << "OASPL(A): "  << sp_dBA<< "dB" << endl;}
        if (sp_dBB<=0 & sp_dBB>=0){stream << "OASPL(B): "  << "N/A"<< endl;} else {stream << "OASPL(B): "  << sp_dBB<< "dB" << endl;}
        if (sp_dBC<=0 & sp_dBC>=0){stream << "OASPL(C): "  << "N/A"<< endl;} else {stream << "OASPL(C): "  << sp_dBC<< "dB" << endl;}
        if (sp_OASPL_alpha<=0 & sp_OASPL_alpha>=0){stream << "SPL_a: "  << "N/A"<< endl;} else {stream << "SPL_a: "  << sp_OASPL_alpha<< endl;}
        if (sp_OASPL_S<=0 & sp_OASPL_S>=0){stream << "SPL_s: "  << "N/A"<< endl;} else {stream << "SPL_s: "  << sp_OASPL_S<< endl;}
         if (sp_OASPL_P<=0 & sp_OASPL_P>=0){stream << "SPL_p: "  << "N/A"<< endl;} else {stream << "SPL_p: "  << sp_OASPL_P<< endl;}
          if (m_parameter.Lowson_type!=0){
        if (sp_LedB<=0 & sp_LedB>=0){stream << "SPL_LE: "  << "N/A"<< endl;} else {stream << "SPL_LE: "  << sp_LedB<< "dB" << endl;}
        if (sp_LedBAW<=0 & sp_LedBAW>=0){stream << "SPL_LE(A): "  << "N/A"<< endl;} else {stream << "SPL_LE(A): "  << sp_LedBAW<< "dB" << endl;}
        if (sp_LedBBW<=0 & sp_LedBAW>=0){stream << "SPL_LE(A): "  << "N/A"<< endl;} else {stream << "SPL_LE(A): "  << sp_LedBBW<< "dB" << endl;}
        if (sp_LedBCW<=0 & sp_LedBAW>=0){stream << "SPL_LE(A): "  << "N/A"<< endl;} else {stream << "SPL_LE(A): "  << sp_LedBCW<< "dB" << endl;}
          }
                stream << endl;
if (i>=(number_of_segments-1)){
stream << "***********************************************************"<<endl;
stream << "Total Blade: "  << endl;
stream << "OASPL: "  << st_OASPL_a<< "dB" << endl;
stream << "OASPL (A): "  << st_dBA_a<< "dB" <<  endl;
stream << "OASPL (B): "  << st_dBB_a<< "dB" <<  endl;
stream << "OASPL (C): "  << st_dBC_a<< "dB" <<  endl;
stream << "SPL_a: "  << st_OASPL_alpha_a<< endl;
stream << "SPL_s: "  << st_OASPL_S_a<< endl;
stream << "SPL_p: "  << st_OASPL_P_a<< endl;
if (m_parameter.Lowson_type!=0){
stream << "SPL_LE: "  << st_LedB_a<< "dB" <<  endl;
stream << "SPL_LE (A): "  << st_LedBAW_a<< "dB" <<  endl;
stream << "SPL_LE (B): "  << st_LedBBW_a<< "dB" <<  endl;
stream << "SPL_LE (C): "  << st_LedBCW_a<< "dB" <<  endl;}
stream << "SPL_Blade_total: " << SPL_blade_total_a << "dB" <<  endl;
stream << "***********************************************************"<<endl;
stream << endl;}

if(z>=lend & i>=(number_of_segments-1)){
stream << "***********************************************************"<<endl;
stream << "TOTAL COMPLETE" << endl;
stream <<"***********************************************************"<<endl;
                stream << "Total: "  << endl;
                stream << "OASPL: "  << st_OASPL<< "dB" << endl;
                stream << "OASPL (A): "  << st_dBA<< "dB" <<  endl;
                stream << "OASPL (B): "  << st_dBB<< "dB" <<  endl;
                stream << "OASPL (C): "  << st_dBC<< "dB" <<  endl;
                stream << "SPL_a: "  << st_OASPL_alpha<< endl;
                stream << "SPL_s: "  << st_OASPL_S<< endl;
                stream << "SPL_p: "  << st_OASPL_P<< endl;
          if (m_parameter.Lowson_type!=0){
                stream << "SPL_LE: "  << st_LedB<< "dB" <<  endl;
          stream << "SPL_LE (A): "  << st_LedBAW<< "dB" <<  endl;
        stream << "SPL_LE (B): "  << st_LedBBW<< "dB" <<  endl;
stream << "SPL_LE (C): "  << st_LedBCW<< "dB" <<  endl;}
                stream << "SPL_Blade_total: " << SPL_blade_total << "dB" <<  endl;
                stream << "***********************************************************"<<endl;
stream << endl;
}
}}}
z=z+ldelta;
        }}

void NoiseSimulation::exportqs3DCalculation(QTextStream &stream)
{
        stream.setRealNumberNotation(QTextStream::FixedNotation);
        stream.setRealNumberPrecision(5);

        const double AWeighting[] = {-44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1,
                                                                           -13.4, -10.9,  -8.6,  -6.6,  -4.8,  -3.2,  -1.9,  -0.8,
                                                                             0.0,   0.6,   1.0,   1.2,   1.3,   1.2,   1.0,   0.5,
                                                                            -0.1,  -1.1,  -2.5,  -4.3,  -6.6,-  9.3};

        const double BWeighting[] = {-20.4, -17.1, -14.2, -11.6,  -9.3,  -7.4,  -5.6,  -4.2,
                                                                            -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.1,   0.0,
                                                                             0.0,   0.0,   0.0,  -0.1,  -0.2,  -0.4,  -0.7,  -1.2,
                                                                            -1.9,  -2.9,  -4.3,  -6.1,  -8.4, -11.1};

        const double CWeighting[] = {-4.4,  -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.2,
                                                                            -0.1,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                                                                             0.0,   0.0,  -0.1,  -0.2,  -0.3,  -0.5,  -0.8,  -1.3,
                                                                            -2.0,  -3.0,  -4.4,  -6.2,  -8.5, -11.2};

        const double Frequency[]= {25, 31.5, 40, 50, 63, 80, 100, 125, 160,200, 250, 315, 400,500, 630, 800, 1000,1250, 1600, 2000, 2500, 3150, 4000, 5000,6300, 8000, 10000, 12500, 16000, 20000};

        stream << "Quasi 3D Noise Log" << endl;
        stream << endl;

NoiseCreatorDialog *pnoisecreatordialog = (NoiseCreatorDialog *) g_mainFrame->m_pSimuWidget;

    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

    SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
    double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
    double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
    double lend  =   pSimuWidget->m_pctrlLELineEdit->getValue();
    double z=lstart;
    double approaxing_wind_speed = m_parameter.originalVelocity;

    QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
    double blade_pitch=pbem->m_pctrlFixedPitch->getValue();

        foreach(BData * bdata, pbem->m_pBEMData->GetBData()){
        int number_of_segments = bdata->m_pos.size();
        double rho = pbem->dlg_rho;
        double dynamic_visc = pbem->dlg_visc;
        double cin_visc = dynamic_visc/rho;
        double K_air = 1.4;
        double R_air = 286.9;
        double T_std_cond = pbem->dlg_temp;
        double P_std_cond = 101300;
        double lambda = pbem->dlg_lambda;
        int mpos_size = bdata->m_pos.size(); //total number of segments
        double finalradius = bdata->m_pos.value(mpos_size-1);
        double nom_tg_speed = bdata->windspeed*lambda;
        double omega = nom_tg_speed/finalradius;
        double rotation = 60/(M_PI*100/nom_tg_speed);
        double c_0_le = 34000;
        double c_const_rd_le = 19./6.;
        double d_const_rd_le = 65.95;
        double c_const_vk_le = 7./3.;
        double d_const_vk_le = 58.4;
        double c_const_le=0;
        double d_const_le=0;
        double aux_le;
        double aux_1_le;
        double aux_4_le;
        double aux_5_le;
        double u_le;
        double c_le;
        double I_le;
        double lambda_le;
        double r_e_le;
        double beta_le;
        double L_le;
        double D_L_le;
        double K_le;
        double S_le;
        double LFC_le;
        double aux0_le;
        double aux1_le;
        double aux4_le;
        double aux5_le;

        QString SPL_alpha_val("");
        QString SPL_dB_S_val("");
        QString SPL_dB_P_val("");
        QString SPL_S_val("");
        QString SPL_P_val("");
        QString SPL_dB_val("");
        QString SPL_A_val("");
        QString SPL_B_val("");
        QString SPL_C_val("");
        QString SPL_LedB_val("");
        QString SPL_LedB_valAW("");
        QString SPL_LedB_valBW("");
        QString SPL_LedB_valCW("");

        //definitions
        double axial_ind_fact[number_of_segments];
        double axial_ind_fact_n[number_of_segments];
        double axial_velocity[number_of_segments];
        double tangential_speed[number_of_segments];
        double resultant_local_speed[number_of_segments];
        double chord[number_of_segments];
        double Reynolds[number_of_segments];
        double Reynolds_BEM[number_of_segments];
        double Reynolds_polar[number_of_segments];
        double Reynolds_error[number_of_segments];
        double Mach_BEM[number_of_segments];
        double Mach_polar[number_of_segments];
        double Mach_error[number_of_segments];
        double alpha[number_of_segments];
        double alpha_BEM[number_of_segments];
        double alpha_polar[number_of_segments];
        double alpha_error[number_of_segments];
        double phi_BEM[number_of_segments];
        double theta_BEM[number_of_segments];
        double cl_cd[number_of_segments];
        double r_R[number_of_segments];
        double c_Rx[number_of_segments];
        double D_starred_C_HT[number_of_segments];
        double D_starred_HT[number_of_segments];
        double D_starred_C_N[number_of_segments];
        double D_starred_N[number_of_segments];
        double Dh[number_of_segments];
        double Dl[number_of_segments];
        double Mach[number_of_segments];
        double corr_fact[number_of_segments];
        double D_starred_HT_S[number_of_segments];
        double D_starred_HT_P[number_of_segments];
        double D_starred_N_S[number_of_segments];
        double D_starred_N_P[number_of_segments];
        double L[number_of_segments];
        double SwAlpha[number_of_segments];
        double SwAlpha_1[number_of_segments];
        double SwAlpha_2[number_of_segments];
        double observer_position = 10;
        double gamma[number_of_segments];
        double gamma0[number_of_segments];
        double beta[number_of_segments];
        double beta0[number_of_segments];
        double gamma0_gamma_min[number_of_segments];
        double gamma0_gamma_plus[number_of_segments];
        double K1[number_of_segments];
        double K2[number_of_segments];
        double EddyMach_calc[number_of_segments];
        double dist_z[number_of_segments];
        double dist_y[number_of_segments];
        double dist_obs[number_of_segments];
        double D_starred_S[number_of_segments];
        double D_starred_P[number_of_segments];
        double first_term_Dh_S[number_of_segments];
        double first_term_Dl_S[number_of_segments];
        double first_term_Dh_P[number_of_segments];
        double St1[number_of_segments];
        double St2[number_of_segments];
        double b0[number_of_segments];
        double B_min_b0[number_of_segments];
        double B_max_b0[number_of_segments];
        double BR_b0[number_of_segments];
        double RCmod[number_of_segments];
        double ao_Rc[number_of_segments];
        double A_min_ao[number_of_segments];
        double A_max_ao[number_of_segments];
        double K1_3[number_of_segments];
        double AR_ao[number_of_segments];
        double St1_bar[number_of_segments];
        double Re_disp_thick[number_of_segments];
        double delta_K1[number_of_segments];
        double phi_rad[number_of_segments];
        double theta_rad[number_of_segments];
        double b[number_of_segments];
        double a[number_of_segments];
        double XRS[number_of_segments];
        double YRS[number_of_segments];
        double ZRS[number_of_segments];
        double XRT[number_of_segments];
        double YRT[number_of_segments];
        double ZRT[number_of_segments];
        double re[number_of_segments];
        double theta[number_of_segments];
        double phi[number_of_segments];
        double calc_int_a[number_of_segments];
        double psi_e[number_of_segments];
        double theta_e[number_of_segments];
        double r_rt[number_of_segments];
        double r_e[number_of_segments];
        double r_1[number_of_segments];
        double c_1[number_of_segments];
        double r_0[number_of_segments];
        double c_0[number_of_segments];
        double twist[number_of_segments];
        double local_twist[number_of_segments];

        double ri[number_of_segments];
        double ri_1[number_of_segments];
        double ri_2[number_of_segments];
        double A[number_of_segments];
        double phi_lin[number_of_segments];

        double DStarXFoilS[number_of_segments];
        double DStarXFoilP[number_of_segments];

        double sp_OASPL_alpha=0;
        double splog_OASPL_alpha=0;
        double st_OASPL_alpha=0;
        double stlog_OASPL_alpha=0;
        double sp_OASPL_S=0;
        double splog_OASPL_S=0;
        double st_OASPL_S=0;
        double stlog_OASPL_S=0;
        double sp_OASPL_P=0;
        double splog_OASPL_P=0;
        double st_OASPL_P=0;
        double stlog_OASPL_P=0;
        double sp_OASPL=0;
        double splog_OASPL=0;
        double st_OASPL=0;
        double stlog_OASPL=0;
        double sp_dBA=0;
        double splog_dBA=0;
        double st_dBA=0;
        double stlog_dBA=0;
        double sp_dBB=0;
        double splog_dBB=0;
        double st_dBB=0;
        double stlog_dBB=0;
        double sp_dBC=0;
        double splog_dBC=0;
        double st_dBC=0;
        double stlog_dBC=0;
        double sp_LedB=0;
        double splog_LedB=0;
        double st_LedB=0;
        double stlog_LedB=0;
        double sp_LedBAW=0;
        double splog_LedBAW=0;
        double st_LedBAW=0;
        double stlog_LedBAW=0;
        double sp_LedBBW=0;
        double splog_LedBBW=0;
        double st_LedBBW=0;
        double stlog_LedBBW=0;
        double sp_LedBCW=0;
        double splog_LedBCW=0;
        double st_LedBCW=0;
        double stlog_LedBCW=0;
        double SPL_blade_total=0;
        double SPL_blade_total_aux=0;

        double r_R0  =  0.05; double c_R0 = 0.05500;
        double r_R1  =  0.25; double c_R1 = 0.07500;
        double r_R2  =  1.00; double c_R2 = 0.02000;

        bool LE_validation;
        bool BPM_validation;

        sp_OASPL_alpha=0;
        sp_OASPL_S=0;
        sp_OASPL_P=0;
        sp_OASPL=0;
        sp_dBA=0;
        sp_dBB=0;
        sp_dBC=0;
        splog_OASPL_alpha=0;
        splog_OASPL_S=0;
        splog_OASPL_P=0;
        splog_OASPL=0;
        splog_dBA=0;
        splog_dBB=0;
        splog_dBC=0;
        splog_LedB=0;
        splog_LedBAW=0;
        splog_LedBBW=0;
        splog_LedBCW=0;

        double sp_OASPL_alpha_a=0;
        double splog_OASPL_alpha_a=0;
        double st_OASPL_alpha_a=0;
        double stlog_OASPL_alpha_a=0;
        double sp_OASPL_S_a=0;
        double splog_OASPL_S_a=0;
        double st_OASPL_S_a=0;
        double stlog_OASPL_S_a=0;
        double sp_OASPL_P_a=0;
        double splog_OASPL_P_a=0;
        double st_OASPL_P_a=0;
        double stlog_OASPL_P_a=0;
        double sp_OASPL_a=0;
        double splog_OASPL_a=0;
        double st_OASPL_a=0;
        double stlog_OASPL_a=0;
        double sp_dBA_a=0;
        double splog_dBA_a=0;
        double st_dBA_a=0;
        double stlog_dBA_a=0;
        double sp_dBB_a=0;
        double splog_dBB_a=0;
        double st_dBB_a=0;
        double stlog_dBB_a=0;
        double sp_dBC_a=0;
        double splog_dBC_a=0;
        double st_dBC_a=0;
        double stlog_dBC_a=0;
        double sp_LedB_a=0;
        double splog_LedB_a=0;
        double st_LedB_a=0;
        double stlog_LedB_a=0;
        double sp_LedBAW_a=0;
        double splog_LedBAW_a=0;
        double st_LedBAW_a=0;
        double stlog_LedBAW_a=0;
        double sp_LedBBW_a=0;
        double splog_LedBBW_a=0;
        double st_LedBBW_a=0;
        double stlog_LedBBW_a=0;
        double sp_LedBCW_a=0;
        double splog_LedBCW_a=0;
        double st_LedBCW_a=0;
        double stlog_LedBCW_a=0;
        double SPL_blade_total_a=0;
        double SPL_blade_total_aux_a=0;

        for (int i = 0; i < number_of_segments; ++i) {

            // definitions
            axial_ind_fact[i] = bdata->m_a_axial.value(i);
            axial_ind_fact_n[i] = bdata->m_a_axial.value(i+1);

            if (i<number_of_segments/2) {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact[i]);}
            else {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact_n[i]);}

            tangential_speed[i] = omega*bdata->m_pos.value(i)*(1.+bdata->m_a_tangential.value(i));
            resultant_local_speed[i] = qSqrt(pow(axial_velocity[i],2)+pow(tangential_speed[i],2));
            chord[i] = bdata->m_c_local.value(i);
            Reynolds[i] = bdata->m_Reynolds.value(i);

            CPolar *pCPolar = (CPolar *) g_mainFrame->m_pctrlPolar;
            Reynolds_BEM[i]=bdata->m_Reynolds.value(i);
            Reynolds_polar[i]=noiseOpPoints[i]->getReynolds();
            Reynolds_error[i]=qFabs(Reynolds_polar[i]-Reynolds_BEM[i])/Reynolds_BEM[i]*100.;
            Mach_polar[i]=noiseOpPoints[i]->getMach();
            Mach[i]=bdata->m_Mach.value(i);
            Mach_BEM[i] = bdata->m_Mach.value(i);
            Mach_error[i]=qFabs(Mach_polar[i]-Mach_BEM[i])/Mach_BEM[i]*100.;
            alpha_BEM[i] = bdata->m_alpha.value(i);
            alpha[i]=alpha_BEM[i];
            phi_BEM[i] = bdata->m_phi.value(i);
            theta_BEM[i] = bdata->m_theta.value(i);
            cl_cd[i] =  bdata->m_LD.value(i);
            r_R[i] = bdata->m_pos.value(i)/finalradius;

            D_starred_N_S[i]=0;
            D_starred_HT_S[i]=0;

//BPM method p 25 manual and spreadsheet
double aux_alpha_polar[noiseOpPoints.size()];
for (int k=0;k<noiseOpPoints.size();++k){
aux_alpha_polar[k]=qFabs(alpha_BEM[i]-noiseOpPoints[k]->getAlphaDegreeAbsolute());
}

double aux_alpha_polar_set=aux_alpha_polar[0];

for (int k=1;k<noiseOpPoints.size();++k){
    if(aux_alpha_polar_set>aux_alpha_polar[k])
   {aux_alpha_polar_set=aux_alpha_polar[k];
alpha_polar[i]=noiseOpPoints[k]->getAlphaDegreeAbsolute();
    }
}
            if (r_R[i] <= r_R0) {c_Rx[i] = c_R0;}
            if ((r_R[i] > r_R0) & (r_R[i] < r_R1)) {c_Rx[i] = (r_R[i]-r_R0)*(c_R1-c_R0)/(0.25-r_R0)+c_R0;}
            if ((r_R[i] <= r_R1) & (r_R[i] >= r_R1)) {c_Rx[i] = c_R1;}
            if ((r_R[i] > r_R1) & (r_R[i] < r_R2)) {c_Rx[i] = (r_R[i]-r_R1)*(c_R2-c_R1)/(r_R2-r_R1)+c_R1;}
            if (r_R[i] >= r_R2) {c_Rx[i] = c_R2;}

            QString c_R= QString::number(c_Rx[i], 'f', 5);
            QString Mach_error_x= QString::number(Mach_error[i], 'f', 2);
            QString Reynolds_error_x= QString::number(Reynolds_error[i], 'f', 2);
            QString alpha_error_x= QString::number(alpha_error[i], 'f', 2);

//heavy tripping

if (Reynolds[i]>300000){
    D_starred_C_HT[i]=pow(10,(3.411-1.5397*log10(Reynolds[i])+0.1059*pow(log10(Reynolds[i]),2)));
}
else {D_starred_C_HT[i]=0.0601*(pow(Reynolds[i],(-0.114)));}

D_starred_HT[i]=chord[i]*D_starred_C_HT[i];

//natural transition
    D_starred_C_N[i]=pow(10,(3.0187-1.5397*log10(Reynolds[i])+0.1059*pow(log10(Reynolds[i]),2)));

D_starred_N[i]=D_starred_C_N[i]*chord[i];

if ((alpha[i]<=0) & (alpha[i]>=0)){
D_starred_HT_S[i]=D_starred_HT[i];
D_starred_HT_P[i]=D_starred_HT[i];
D_starred_N_S[i]=D_starred_N[i];
D_starred_N_P[i]=D_starred_N[i];
}
else{
//alpha !=0 pressure side
if ((alpha[i]<0) || (alpha[i]>0)){
corr_fact[i]=pow(10,(-0.0432*alpha[i]+0.00113*pow(alpha[i],2)));
D_starred_HT_P[i]=D_starred_HT[i]*corr_fact[i];
D_starred_N_P[i]=D_starred_N[i]*corr_fact[i];
}

//alpha !=0 suction side heavy tripping
if ((alpha[i]>0) & (alpha[i]<=5)){
corr_fact[i]=pow(10,(0.0679*alpha[i]));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

if ((alpha[i]>5) & (alpha[i]<=12.5)){
corr_fact[i]=0.381*(pow(10,(0.1516*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

if ((alpha[i]>12.5) & (alpha[i]<=25)){
corr_fact[i]=14.296*(pow(10,(0.0258*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

//alpha !=0 suction side natural transition
if ((alpha[i]>0) & (alpha[i]<=7.5)){
corr_fact[i]=pow(10.,(0.0679*alpha[i]));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

if ((alpha[i]>7.5) & (alpha[i]<=12.5)){
corr_fact[i]=0.0162*(pow(10.,(0.3066*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

if ((alpha[i]>12.5) & (alpha[i]<=25)){
corr_fact[i]=54.42*(pow(10.,(0.0258*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}}

//For D* Xfoil

DStarXFoilS[i]=m_calculation.m_DStarInterpolatedS3d[i];
DStarXFoilP[i]=m_calculation.m_DStarInterpolatedP3d[i];

//Length of Wetted Trailing Edge
L[i]=bdata->m_pos.value(i+1)-bdata->m_pos.value(i);

//Calculate the Switching Angle
SwAlpha_1[i]=23.43*Mach[i]+4.651;
SwAlpha_2[i]=12.5;

if (SwAlpha_1[i]<SwAlpha_2[i]){SwAlpha[i]=SwAlpha_1[i];}
else {SwAlpha[i]=SwAlpha_2[i];}

double EddyMach = m_parameter.eddyConvectionMach;

if (m_parameter.lowFreq) {
Dh[i]=(2.*pow(sin(qDegreesToRadians(theta[i]/2.)),2)*pow(sin(qDegreesToRadians(phi[i])),2))/pow(1+Mach[i]*cos(qDegreesToRadians(theta[i]))*(1.+(Mach[i]-Mach[i]*EddyMach)*cos(qDegreesToRadians(phi[i]))),2);
}
else{Dh[i]=1;}

if (m_parameter.lowFreq) {
Dl[i]=(2.*pow(sin(qDegreesToRadians(theta[i])),2)*pow(sin(qDegreesToRadians(phi[i])),2))/pow((1+(Mach[i]*cos(qDegreesToRadians(theta[i])))),4);
}
else{Dl[i]=1;}

gamma[i]=27.094*Mach[i]+3.32;

gamma0[i]=SwAlpha_1[i];

beta[i]=72.65*Mach[i]+10.74;

beta0[i]=-34.19*Mach[i]-13.82;

gamma0_gamma_min[i]=gamma0[i]-gamma[i];

gamma0_gamma_plus[i]=gamma0[i]+gamma[i];

if (bdata->m_Reynolds.value(i)<247000)
{K1[i]=-4.31*log10(bdata->m_Reynolds.value(i))+156.3;}
else if (bdata->m_Reynolds.value(i)>800000)
{K1[i]=128.5;}
else {K1[i]=-9.*log10(bdata->m_Reynolds.value(i))+181.6;}

if (alpha[i]<gamma0_gamma_min[i])
{K2[i]=K1[i]-1000.;}
else if (alpha[i]>gamma0_gamma_plus[i])
{K2[i]=K1[i]-12.;}
else
{K2[i]=K1[i]+(sqrt(pow(beta[i],2)-pow((beta[i]/gamma[i]),2)*pow((alpha[i]-gamma0[i]),2)))+beta0[i];}

double EddyMach_perc=EddyMach;

EddyMach_calc[i]=Mach[i]*EddyMach_perc;

//delta starred type, if natural transition or heavy-tripping
if (m_parameter.dstar_type==0){
    D_starred_S[i]=D_starred_N_S[i];
    D_starred_P[i]=D_starred_N_P[i];
}

else if (m_parameter.dstar_type==1){
FoilPolarDlg *pFoilPolarDlg = (FoilPolarDlg *) g_mainFrame->m_pctrlXDirectWidget;

    double TopTrip=pFoilPolarDlg->m_XTopTr;
    double BotTrip=pFoilPolarDlg->m_XBotTr;

if((TopTrip<=0 & TopTrip>=0) & (BotTrip<=0 & BotTrip>=0)) {
//        natural transition
    D_starred_S[i]=D_starred_HT_S[i];
    D_starred_P[i]=D_starred_HT_P[i];
}
else {
//heavy tripping
    D_starred_S[i]=DStarXFoilS[i];
    D_starred_P[i]=DStarXFoilP[i];
}}
else if (m_parameter.dstar_type==2){
//user
    NoiseParameter *pNoiseParameter = (NoiseParameter *) g_mainFrame->m_pSimuWidget;
    D_starred_S[i]=pNoiseParameter->D_starred_S_user[i];
    D_starred_P[i]=pNoiseParameter->D_starred_P_user[i];
}

double B=0;
double XB=0;
double YB=0;
double ZB=0;

//phi type, fixed 90º or free
if (m_parameter.phi_type==0){
    //by the quasi-3D spreadsheet
    ri[i]=bdata->m_pos.value(i);
    ri_1[i]=bdata->m_pos.value(i+1);
    ri_2[i]=bdata->m_pos.value(i+2);
    B=bdata->m_pos.value(number_of_segments-1)/2.;
    A[i]=ri_1[i]+(ri_2[i]-ri_1[i])/2.;
    re[i]=sqrt(pow((A[i]-B),2)+pow(m_parameter.distanceObsever,2));
    phi_lin[i]=qRadiansToDegrees(qAsin(m_parameter.distanceObsever/re[i]));
    phi[i]=180.-phi_lin[i];
    theta[i]=90;
    dist_obs[i]=ri_1[i];
}
else if (m_parameter.phi_type==1){
    //    re phi and theta calculation p 77 C_Project_Log_Text_Jan_16.pdf
    //    Input X e , Y e , Z e
    //    Attribute their respective values to X B , Y B , Z B
        XB=m_parameter.obs_x_pos;
        YB=m_parameter.obs_y_pos;
        ZB=m_parameter.obs_z_pos;

    //For the blade, get the Pitch angle, θp - blade pitch

    //For each blade segment, get: C r i+1 , C r i , r i+1 , r i , β
    //Cri+1 = c_1
    //r+1 = r_1
    //Cri = c_0
    //ri = r_0
    //beta = local_twist

    if((i<=number_of_segments) & (i>=number_of_segments)){
c_1[i]=0;
r_1[i]=0;
    }
    else
    {
c_1[i]=bdata->m_c_local.value(i+1);
r_1[i]=bdata->m_pos.value(i+1);
    }}

c_0[i]=bdata->m_c_local.value(i);
r_0[i]=bdata->m_pos.value(i);

local_twist[i]=theta_BEM[i];

    b[i]=qRadiansToDegrees(qAtan((c_1[i]-c_0[i])/(r_1[i]-r_0[i])));

//    the angle a is the total angle between the Y B Z B blade reference system plane and the local midsection chord line
    a[i]=local_twist[i]+blade_pitch;

    XRS[i]=XB*cos(qDegreesToRadians(a[i]))+YB*sin(qDegreesToRadians(a[i]));
    YRS[i]=-XB*sin(qDegreesToRadians(a[i]))+YB*cos(qDegreesToRadians(a[i]));
    ZRS[i]=ZB-(bdata->m_pos.value(i)-bdata->m_pos.value(i-1))/2.;

//    Calculate Y RS − 0.75 ∗ (C r i+1 − C r i )/2
    calc_int_a[i]=(YRS[i]-0.75*(c_1[i]-c_0[i])/2.);

    XRT[i]=XRS[i];
    YRT[i]=cos(qDegreesToRadians(b[i]))*calc_int_a[i]+sin(qDegreesToRadians(b[i]))*ZRS[i];
    ZRT[i]=-sin(qDegreesToRadians(b[i]))*calc_int_a[i]+cos(qDegreesToRadians(b[i]))*ZRS[i];

    r_e[i]=sqrt(pow(XRT[i],2)+pow(YRT[i],2)+pow(ZRT[i],2));
    r_rt[i]=r_e[i];
    theta_e[i]=qRadiansToDegrees(qAtan(ZRT[i]/YRT[i]));

    psi_e[i]=qRadiansToDegrees(qAtan(XRT[i]/ZRT[i]));

    dist_obs[i]=re[i];

   phi_rad[i]=qDegreesToRadians(phi[i]);
   theta_rad[i]=qDegreesToRadians(theta[i]);

   alpha_error[i]=qFabs(alpha_polar[i]-alpha_BEM[i])/alpha_BEM[i]*100.;

first_term_Dh_S[i]=10.*log10(pow(Mach[i],5)*L[i]*Dh[i]*D_starred_S[i]/pow(dist_obs[i],2));

first_term_Dl_S[i]=10.*log10(pow(Mach[i],5)*L[i]*Dl[i]*D_starred_S[i]/pow(dist_obs[i],2));

first_term_Dh_P[i]=10.*log10(pow(Mach[i],5)*L[i]*Dh[i]*D_starred_P[i]/pow(dist_obs[i],2));

St1[i]=0.02*(pow(Mach[i],-0.6));

if(alpha[i]<1.33){St2[i]=St1[i];}
else if(alpha[i]>12.5){St2[i]=St1[i]*4.72;}
else {St2[i]=St1[i]*pow(10.,(0.0054*pow((alpha[i]-1.33),2)));}

if (bdata->m_Reynolds.value(i)<95200)
{b0[i]= 0.3;}
else if (bdata->m_Reynolds.value(i)>857000)
{b0[i]= 0.56;}
else {b0[i]=-4.48*pow(10,-13)*(pow((bdata->m_Reynolds.value(i)-857000.),2)+0.56);}

if (b0[i]<0.13)
{B_min_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
else if (b0[i]>0.145)
{B_min_b0[i]=-817.81*pow(b0[i],3)+335.21*pow(b0[i],2)-135.024*b0[i]+10.619;}
else{B_min_b0[i]=-83.607*b0[i]+8.138;}

if (b0[i]<0.1)
{B_max_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
else if (b0[i]>0.187)
{B_max_b0[i]=-80.541*pow(b0[i],3)+44.174*pow(b0[i],2)-39.381*b0[i]+2.344;}
else {B_max_b0[i]=-31.33*b0[i]+1.854;}

BR_b0[i]=(-20-B_min_b0[i])/(B_max_b0[i]-B_min_b0[i]);

RCmod[i]=3*Reynolds[i];

if (RCmod[i]<95200){ao_Rc[i]=0.57;}
else if (RCmod[i]>857000){ao_Rc[i]=1.13;}
else {ao_Rc[i]=-9.57*pow(10,-13)*(pow((RCmod[i]-857000),2)+1.13);}

if(ao_Rc[i]<0.204){A_min_ao[i]=sqrt(67.552-886.788*pow(ao_Rc[i],2))-8.219;}
else if (ao_Rc[i]>0.244){A_min_ao[i]=-142.795*pow(ao_Rc[i],3)+103.656*pow(ao_Rc[i],2)-57.757*ao_Rc[i]+6.006;}
else {A_min_ao[i]=-32.665*ao_Rc[i]+3.981;}

if (ao_Rc[i]<0.13){A_max_ao[i]=sqrt(67.552-886.788*pow(ao_Rc[i],2))-8.219;}
else if(ao_Rc[i]>0.321){A_max_ao[i]=-4.669*pow(ao_Rc[i],3)+3.491*pow(ao_Rc[i],2)-16.699*ao_Rc[i]+1.149;}
else {A_max_ao[i]=-15.901*ao_Rc[i]+1.098;}

K1_3[i]=K1[i]-3.;

AR_ao[i]=(-20-A_min_ao[i])/(A_max_ao[i]-A_min_ao[i]);

St1_bar[i]=(St1[i]+St2[i])/2.;

Re_disp_thick[i]=1.225*bdata->m_Windspeed.value(i)*D_starred_P[i]/(0.0000178);

if (Re_disp_thick[i]>5000){delta_K1[i]=0;}
else {delta_K1[i]=alpha[i]*(1.43*log10(Re_disp_thick[i])-5.29);}

    double b0[number_of_segments];
    if (bdata->m_Reynolds.value(i)<95200)
    {b0[i]= 0.3;}
    else if (bdata->m_Reynolds.value(i)>857000)
    {b0[i]= 0.56;}
    else {b0[i]=-4.48*pow(10.,-13)*(pow((bdata->m_Reynolds.value(i)-857000.),2)+0.56);}

    double B_min_b0[number_of_segments];
    if (b0[i]<0.13)
    {B_min_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
    else if (b0[i]>0.145)
    {B_min_b0[i]=-817.81*pow(b0[i],3)+335.21*pow(b0[i],2)-135.024*b0[i]+10.619;}
    else{B_min_b0[i]=-83.607*b0[i]+8.138;}

    double B_max_b0[number_of_segments];
    if (b0[i]<0.1)
    {B_max_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
    else if (b0[i]>0.187)
    {B_max_b0[i]=-80.541*pow(b0[i],3)+44.174*pow(b0[i],2)-39.381*b0[i]+2.344;}
    else {B_max_b0[i]=-31.33*b0[i]+1.854;}

    double BR_b0[number_of_segments];
    BR_b0[i]=(-20-B_min_b0[i])/(B_max_b0[i]-B_min_b0[i]);

    if((z<=m_parameter.TSRtd) & (z>=m_parameter.TSRtd)){
    QString str= QString::number(z, 'f', 1);


//        stream << "Angles less than the switching angle: "  << endl;
        stream << "Tip Speed Ratio: " << str << endl;
        stream << "Section: " << (i+1)<<"/"<<number_of_segments << endl;
        stream << "Reynolds: " << Reynolds[i] << endl;
        stream << "alpha " << alpha[i] << endl;
        stream << endl;

if(m_parameter.Lowson_type!=0){
    stream << qSetFieldWidth(14)  <<
              "Freq[Hz]"  << ";" <<
//              "Sts"  << ";" <<
//              "b"  << ";" <<
//              "B_min(b)" << ";" <<
//              "B_max(b)" << ";" <<
//              "B(b)"  << ";" <<
//              "a_alpha"  << ";" <<
//              "A_min(a)_alpha" << ";" <<
//              "A_max(a)_alpha" << ";" <<
//              "A'(a)"  << ";" <<
//              "A-Weighting"  << ";" <<
//              "B-Weighting" << ";" <<
//              "C-Weighting"  << ";" <<
//              "SPL_alpha < 0"  << ";" <<
//              "SPL_alpha > 0"  << ";" <<
//              "dB_A_alpha S"  << ";" <<
//              "dB_A_alpha P"  << ";" <<
//              "dB_B_alpha S"  << ";" <<
//              "dB_B_alpha P"  << ";" <<
//              "dB_C_alpha S"  << ";" <<
//              "dB_C_alpha P"  << ";" <<
//              "St1_bar"  << ";" <<
//              "a_S"  << ";" <<
//              "A_min_a_S"  << ";" <<
//              "A_max_a_S"  << ";" <<
//              "A_a_S"  << ";" <<
//              "SPL_dB_S"  << ";" <<
//              "dB_A_S"  << ";" <<
//              "dB_B_S"  << ";" <<
//              "dB_C_S"  << ";" <<
//              "Sts_St1_bar"  << ";" <<
//              "Stp"  << ";" <<
//              "a_P"  << ";" <<
//              "A_min_a_P"  << ";" <<
//              "A_max_a_P"  << ";" <<
//              "A_a_P"  << ";" <<
//              "SPL_dB_P"  << ";" <<
//              "dB_A_P"  << ";" <<
//              "dB_B_P"  << ";" <<
//              "dB_C_P"  << ";" <<
              "SPL (dB)"  << ";" <<
              "SPLa"  << ";" <<
              "SPLs"  << ";" <<
              "SPLp"  << ";" <<
              "SPL (dB(A))"  << ";" <<
              "SPL (dB(B))"  << ";" <<
              "SPL (dB(C))"  << ";" <<
              "SPL_LE (dB) "  << ";" <<
              "SPL_LE ((dB(A)) "  << ";" <<
              "SPL_LE ((dB(B)) "  << ";" <<
              "SPL_LE (dB(C)) "  << ";" <<
//              "s_log_SPL_alpha"  << ";" <<
//              "s_log_SPL_S"  << ";" <<
//              "s_log_SPL_P"  << ";" <<
//              "s_log_SPL"  << ";" <<
//              "s_log_dBA"  << ";" <<
//              "s_log_dBB"  << ";" <<
//              "s_log_dBC"  << ";" <<
//              "s_log_LE_dB"  << ";" <<
              endl;}
else{
        stream << qSetFieldWidth(14)  <<
                  "Freq[Hz]"  << ";" <<
//                  "Sts"  << ";" <<
//                  "b"  << ";" <<
//                  "B_min(b)" << ";" <<
//                  "B_max(b)" << ";" <<
//                  "B(b)"  << ";" <<
//                  "a_alpha"  << ";" <<
//                  "A_min(a)_alpha" << ";" <<
//                  "A_max(a)_alpha" << ";" <<
//                  "A'(a)"  << ";" <<
//                  "A-Weighting"  << ";" <<
//                  "B-Weighting" << ";" <<
//                  "C-Weighting"  << ";" <<
//                  "SPL_alpha < 0"  << ";" <<
//                  "SPL_alpha > 0"  << ";" <<
//                  "dB_A_alpha S"  << ";" <<
//                  "dB_A_alpha P"  << ";" <<
//                  "dB_B_alpha S"  << ";" <<
//                  "dB_B_alpha P"  << ";" <<
//                  "dB_C_alpha S"  << ";" <<
//                  "dB_C_alpha P"  << ";" <<
//                  "St1_bar"  << ";" <<
//                  "a_S"  << ";" <<
//                  "A_min_a_S"  << ";" <<
//                  "A_max_a_S"  << ";" <<
//                  "A_a_S"  << ";" <<
//                  "SPL_dB_S"  << ";" <<
//                  "dB_A_S"  << ";" <<
//                  "dB_B_S"  << ";" <<
//                  "dB_C_S"  << ";" <<
//                  "Sts_St1_bar"  << ";" <<
//                  "Stp"  << ";" <<
//                  "a_P"  << ";" <<
//                  "A_min_a_P"  << ";" <<
//                  "A_max_a_P"  << ";" <<
//                  "A_a_P"  << ";" <<
//                  "SPL_dB_P"  << ";" <<
//                  "dB_A_P"  << ";" <<
//                  "dB_B_P"  << ";" <<
//                  "dB_C_P"  << ";" <<
                  "SPL (dB)"  << ";" <<
                  "SPLa"  << ";" <<
                  "SPLs"  << ";" <<
                  "SPLp"  << ";" <<
                  "SPL (dB(A))"  << ";" <<
                  "SPL (dB(B))"  << ";" <<
                  "SPL (dB(C))"  << ";" <<
//                  "s_log_SPL_alpha"  << ";" <<
//                  "s_log_SPL_S"  << ";" <<
//                  "s_log_SPL_P"  << ";" <<
//                  "s_log_SPL"  << ";" <<
//                  "s_log_dBA"  << ";" <<
//                  "s_log_dBB"  << ";" <<
//                  "s_log_dBC"  << ";" <<
                  endl;}

    int w=30;

    double slog_SPL_alpha[w];
    double slog_SPL_S[w];
    double slog_SPL_P[w];
    double slog_SPL[w];
    double slog_dBA[w];
    double slog_dBB[w];
    double slog_dBC[w];
    double slog_LedB[w];
    double slog_LedBAW[w];
    double slog_LedBBW[w];
    double slog_LedBCW[w];

     for (int j = 0; j < w; ++j) {

         double Sts[w];
         Sts[j]=Frequency[j]*D_starred_S[i]/bdata->m_Windspeed.value(i);

         double b_alpha[w];
         b_alpha[j]=qFabs(log10(Sts[j]/St2[i]));

         double B_min[w];
    if (b_alpha[j]<0.13)
    {B_min[j]=sqrt(16.888-886.788*pow(b_alpha[j],2));}
    else if(b_alpha[j]>0.145)
    {B_min[j]=-817.81*pow(b_alpha[j],3)+335.21*pow(b_alpha[j],2)-135.024*b_alpha[j]+10.619;}
    else {B_min[j]=-83.607*b_alpha[j]+8.138;}

    double B_max[w];
    if (b_alpha[j]<0.1)
    {B_max[j]=sqrt(16.888-886.788*pow(b_alpha[j],2))-4,109;}
    else if(b_alpha[j]>0.187)
    {B_max[j]=-80.541*pow(b_alpha[j],3)+44.174*pow(b_alpha[j],2)-39.381*b_alpha[j]+2.344;}
    else {B_max[j]=-31.33*b_alpha[j]+1.854;}

    double B_b[w];
    B_b[j]=B_min[j]+BR_b0[i]*(B_max[j]-B_min[j]);
    if(qIsInf(B_b[j])){B_b[j]=0;}

    double a_alpha[w];
    a_alpha[j]=qFabs(log10(Sts[j]/St2[i]));

    double A_min_alpha[w];
    if (a_alpha[j]<0.204)
    {A_min_alpha[j]=sqrt(67.552-886.788*pow(a_alpha[j],2))-8.219;}
    else if(a_alpha[j]>0.244)
    {A_min_alpha[j]=-142.795*pow(a_alpha[j],3)+103.656*pow(a_alpha[j],2)-57.757*a_alpha[j]+6.006;}
    else {A_min_alpha[j]=-32.665*a_alpha[j]+3.981;}

    double A_max_alpha[w];
    if (a_alpha[j]<0.13)
    {A_max_alpha[j]=sqrt(67.552-886.788*pow(a_alpha[j],2))-8.219;}
    else if(a_alpha[j]>0.321)
    {A_max_alpha[j]=-4.669*pow(a_alpha[j],3)+3.491*pow(a_alpha[j],2)-16.699*a_alpha[j]+1.149;}
    else {A_max_alpha[j]=-15.901*a_alpha[j]+1.098;}

    double Alin_a[w];
    Alin_a[j]=A_min_alpha[j]+AR_ao[i]*(A_max_alpha[j]-A_min_alpha[j]);

    double SPL_alpha_min0[w];
    SPL_alpha_min0[j]=first_term_Dh_S[i]+K2[i]+B_b[j];

    double SPL_alpha_big0[w];
    SPL_alpha_big0[j]=first_term_Dl_S[i]+K2[i]+Alin_a[j];

    double dBA_alpha_min0[w];
    dBA_alpha_min0[j]=SPL_alpha_min0[j]+AWeighting[j];

    double dBA_alpha_big0[w];
    dBA_alpha_big0[j]=SPL_alpha_big0[j]+AWeighting[j];

    double dBB_alpha_min0[w];
    dBB_alpha_min0[j]=SPL_alpha_min0[j]+BWeighting[j];

    double dBC_alpha_big0[w];
    dBC_alpha_big0[j]=SPL_alpha_big0[j]+CWeighting[j];

    double dBC_alpha_min0[w];
    dBC_alpha_min0[j]=SPL_alpha_min0[j]+CWeighting[j];

    double dBB_alpha_big0[w];
    dBB_alpha_big0[j]=SPL_alpha_big0[j]+BWeighting[j];

    double St1_bar[w];
    St1_bar[j]=(St1[i]+St2[i])/2.;

    double a_S[w];
    a_S[j]=qFabs(log10(Sts[j]/St1_bar[j]));

    double A_min_S[w];
    if (a_S[j]<0.204){A_min_S[j]=sqrt(67.552-886.788*pow(a_S[j],2))-8.219;}
    else if (a_S[j]>0.244){A_min_S[j]=-142.795*pow(a_S[j],3)+103.656*pow(a_S[j],2)-57.757*a_S[j]+6.006;}
    else {A_min_S[j]=-32.665*a_S[j]+3.981;}

    double A_max_S[w];
    if (a_S[j]<0.13){A_max_S[j]=sqrt(67.552-886.788*pow(a_S[j],2))-8.219;}
    else if (a_S[j]>0.321){A_max_S[j]=-4.669*pow(a_S[j],3)+3.491*pow(a_S[j],2)-16.699*a_S[j]+1.149;}
    else {A_max_S[j]=-15.901*a_S[j]+1.098;}

    double A_a_S[w];
    A_a_S[j]=A_min_S[j]+AR_ao[i]*(A_max_S[j]-A_min_S[j]);

    double SPL_dB_S[w];
    SPL_dB_S[j]=first_term_Dh_S[i]+A_a_S[j]+K1_3[i];

    double dBA_S[w];
    dBA_S[j]=SPL_dB_S[j]+AWeighting[j];

    double dBB_S[w];
    dBB_S[j]=SPL_dB_S[j]+BWeighting[j];

    double dBC_S[w];
    dBC_S[j]=SPL_dB_S[j]+CWeighting[j];

    double Sts_St1_bar[w];
    Sts_St1_bar[j]=Sts[j]/St1_bar[i];

    double Stp_P[w];
    Stp_P[j]=Frequency[j]*D_starred_P[i]/bdata->m_Windspeed.value(i);

    double a_P[w];
    a_P[j]=qFabs(log10(Stp_P[j]/St1[i]));

    double A_min_P[w];
    if(a_P[j]<0.204){A_min_P[j]=sqrt(67.552-886.788*pow(a_P[j],2))-8.219;}
    else if (a_P[j]>0.244){A_min_P[j]=-142.795*pow(a_P[j],3)+103.656*pow(a_P[j],2)-57.757*a_P[j]+6.006;}
    else {A_min_P[j]=-32.665*a_P[j]+3.981;}

    double A_max_P[w];
    if(a_P[j]<0.13){A_max_P[j]=sqrt(67.552-886.788*pow(a_P[j],2))-8.219;}
    else if (a_P[j]>0.321){A_max_P[j]=-4.669*pow(a_P[j],3)+3.491*pow(a_P[j],2)-16.699*a_P[j]+1.149;}
    else {A_max_P[j]=-15.901*a_P[j]+1.098;}

    double A_a_P[w];
    A_a_P[j]=A_min_P[j]+AR_ao[i]*(A_max_P[j]-A_min_P[j]);

    double SPL_dB_P[w];
    SPL_dB_P[j]=delta_K1[i]+A_a_P[j]+K1_3[i]+first_term_Dh_P[i];

    double dBA_P[w];
    dBA_P[j]=SPL_dB_P[j]+AWeighting[j];

    double dBB_P[w];
    dBB_P[j]=SPL_dB_P[j]+BWeighting[j];

    double dBC_P[w];
    dBC_P[j]=SPL_dB_P[j]+CWeighting[j];

    double SPL_alpha[w];
    if (alpha[i]<SwAlpha[i]){SPL_alpha[j]=SPL_alpha_min0[j];}
    else {SPL_alpha[j]=SPL_alpha_big0[j];}

    double SPL_S[w];
    if (alpha[i]<SwAlpha[i]){SPL_S[j]=SPL_dB_S[j];}
    else {SPL_S[j]=-999999999999.;}

    double SPL_P[w];
    if (alpha[i]<SwAlpha[i]){SPL_P[j]=SPL_dB_P[j];}
    else {SPL_P[j]=-999999999999.;}

    double SPL_dB[w];
    SPL_dB[j]=10*log10(pow(10.,(SPL_alpha[j]/10.))+pow(10.,(SPL_S[j]/10.))+pow(10.,(SPL_P[j]/10.)));

    double SPL_A[w];
    SPL_A[j]=SPL_dB[j]+AWeighting[j];

    double SPL_B[w];
    SPL_B[j]=SPL_dB[j]+BWeighting[j];

    double SPL_C[w];
    SPL_C[j]=SPL_dB[j]+CWeighting[j];

    double SPL_LedB[w];
    SPL_LedB[w]=0;

    double SPL_LedBAW[w];
    SPL_LedBAW[w]=0;

    double SPL_LedBBW[w];
    SPL_LedBBW[w]=0;

    double SPL_LedBCW[w];
    SPL_LedBCW[w]=0;

    u_le=m_parameter.u_wind_speed*100.;
    c_le=100.*chord[i];
    I_le=m_parameter.TurbulenceIntensity;
    lambda_le=100.*m_parameter.IntegralLengthScale;
    r_e_le=100.*dist_obs[i];
    beta_le=sqrt(1-pow(Mach[i],2));
    L_le=100.*L[i];
    D_L_le=0.5*Dl[i];
    aux_le=0.5*(lambda_le*L_le*pow(rho, 2)*pow(c_0_le, 2)*pow(u_le, 2)*pow(Mach[i], 3)*pow(I_le, 2)*D_L_le)/(pow(r_e_le, 2));
    K_le=M_PI*Frequency[i]*c_le/u_le;
    S_le=sqrt(pow((2.*M_PI*K_le/(pow(beta_le, 2)))+(pow((1+(2.4*K_le/pow(beta_le,2))),-1)),-1));
    LFC_le = 10.*Mach[i]*pow(S_le*K_le/beta_le,2);

    if(m_parameter.Lowson_type==2){
        c_const_le = c_const_rd_le;
        d_const_le = d_const_rd_le;
     }
    else if(m_parameter.Lowson_type==1){
            c_const_le = c_const_vk_le;
            d_const_le = d_const_vk_le;
         }
        else {
            c_const_le=0;
            d_const_le = 0;
        }

aux0_le=0.5*(lambda_le*L_le*pow(rho, 2)*pow(c_0_le, 2)*pow(u_le, 2)*pow(Mach[i], 3)*pow(I_le, 2)*Dl[i])/(pow(r_e_le, 2));
aux1_le=10.*log10(pow(LFC_le/(1+LFC_le), 2))+d_const_le;
aux4_le=pow(K_le,3)/pow(1+(pow(K_le,2)),c_const_le);
aux5_le=10.*log10(aux0_le*aux4_le);

//Validation:

//Lowson validation:
if (((((m_parameter.Lowson_type!=0) & (Mach[i]<=0.18)) & (Mach[i]>0) & (Reynolds[i]<=(6.*pow(10,5)))) & (Reynolds[i]>0))){
SPL_LedB[j]=10.*log10(pow(10,(aux1_le+aux5_le)/10.));
SPL_LedBAW[j]=SPL_LedB[j]+AWeighting[j];
SPL_LedBBW[j]=SPL_LedB[j]+BWeighting[j];
SPL_LedBCW[j]=SPL_LedB[j]+CWeighting[j];
LE_validation=true;
}
else{
SPL_LedB[j]=-999999999999.;
SPL_LedBAW[j]=-999999999999.;
SPL_LedBBW[j]=-999999999999.;
SPL_LedBCW[j]=-999999999999.;
LE_validation=false;
}

//BPM validation:
//p 17 C_Project_Log_Text_15_jan_16
if(((((alpha[i]<=19.8 & Reynolds[i]<3*pow(10,6)) & (Mach[i]<0.21)) & (Reynolds[i]>0)) & (Mach[i]>0))){
BPM_validation=true;
}
else{
SPL_alpha[j]=-999999999999.;
SPL_S[j]=-999999999999.;
SPL_P[j]=-999999999999.;
SPL_dB[j]=-999999999999.;
SPL_A[j]=-999999999999.;
SPL_B[j]=-999999999999.;
SPL_C[j]=-999999999999.;
BPM_validation=false;
}

    slog_SPL_alpha[j]=pow(10.,(SPL_alpha[j]/10.));
    slog_SPL_S[j]=pow(10.,(SPL_S[j]/10.));
    slog_SPL_P[j]=pow(10.,(SPL_P[j]/10.));
    slog_SPL[j]=pow(10.,(SPL_dB[j]/10.));
    slog_dBA[j]=pow(10.,(SPL_A[j]/10.));
    slog_dBB[j]=pow(10.,(SPL_B[j]/10.));
    slog_dBC[j]=pow(10.,(SPL_C[j]/10.));
    slog_LedB[j]=pow(10.,(SPL_LedB[j]/10.));
    slog_LedBAW[j]=pow(10.,(SPL_LedBAW[j]/10.));
    slog_LedBBW[j]=pow(10.,(SPL_LedBBW[j]/10.));
    slog_LedBCW[j]=pow(10.,(SPL_LedBCW[j]/10.));

    if(qIsNaN(slog_SPL_alpha[j])){slog_SPL_alpha[j]=0;}
    if(qIsNaN(slog_SPL_S[j])){slog_SPL_S[j]=0;}
    if(qIsNaN(slog_SPL_P[j])){slog_SPL_P[j]=0;}
    if(qIsNaN(slog_SPL[j])){slog_SPL[j]=0;}
    if(qIsNaN(slog_dBA[j])){slog_dBA[j]=0;}
    if(qIsNaN(slog_dBB[j])){slog_dBB[j]=0;}
    if(qIsNaN(slog_dBC[j])){slog_dBC[j]=0;}
    if(qIsNaN(slog_LedB[j])){slog_LedB[j]=0;}
    if(qIsNaN(slog_LedBAW[j])){slog_LedBAW[j]=0;}
    if(qIsNaN(slog_LedBBW[j])){slog_LedBBW[j]=0;}
    if(qIsNaN(slog_LedBCW[j])){slog_LedBCW[j]=0;}

    splog_OASPL_alpha+=slog_SPL_alpha[j];
    splog_OASPL_S+=slog_SPL_S[j];
    splog_OASPL_P+=slog_SPL_P[j];
    splog_OASPL+=slog_SPL[j];
    splog_dBA+=slog_dBA[j];
    splog_dBB+=slog_dBB[j];
    splog_dBC+=slog_dBC[j];
    splog_LedB+=slog_LedB[j];
    splog_LedBAW+=slog_LedBAW[j];
    splog_LedBBW+=slog_LedBBW[j];
    splog_LedBCW+=slog_LedBCW[j];

    splog_OASPL_alpha_a+=slog_SPL_alpha[j];
    splog_OASPL_S_a+=slog_SPL_S[j];
    splog_OASPL_P_a+=slog_SPL_P[j];
    splog_OASPL_a+=slog_SPL[j];
    splog_dBA_a+=slog_dBA[j];
    splog_dBB_a+=slog_dBB[j];
    splog_dBC_a+=slog_dBC[j];
    splog_LedB_a+=slog_LedB[j];
    splog_LedBAW_a+=slog_LedBAW[j];
    splog_LedBBW_a+=slog_LedBBW[j];
    splog_LedBCW_a+=slog_LedBCW[j];

    sp_OASPL_alpha=10*log10(splog_OASPL_alpha);
    sp_OASPL_S=10*log10(splog_OASPL_S);
    sp_OASPL_P=10*log10(splog_OASPL_P);
    sp_OASPL=10*log10(splog_OASPL);
    sp_dBA=10*log10(splog_dBA);
    sp_dBB=10*log10(splog_dBB);
    sp_dBC=10*log10(splog_dBC);
    sp_LedB=10*log10(splog_LedB);
    sp_LedBAW=10*log10(splog_LedBAW);
    sp_LedBBW=10*log10(splog_LedBBW);
    sp_LedBCW=10*log10(splog_LedBCW);

    sp_OASPL_alpha_a=10*log10(splog_OASPL_alpha_a);
    sp_OASPL_S_a=10*log10(splog_OASPL_S_a);
    sp_OASPL_P_a=10*log10(splog_OASPL_P_a);
    sp_OASPL_a=10*log10(splog_OASPL_a);
    sp_dBA_a=10*log10(splog_dBA_a);
    sp_dBB_a=10*log10(splog_dBB_a);
    sp_dBC_a=10*log10(splog_dBC_a);
    sp_LedB_a=10*log10(splog_LedB_a);
    sp_LedBAW_a=10*log10(splog_LedBAW_a);
    sp_LedBBW_a=10*log10(splog_LedBBW_a);
    sp_LedBCW_a=10*log10(splog_LedBCW_a);

    if(qIsInf(sp_OASPL_alpha)){sp_OASPL_alpha=0;}
    if(qIsInf(sp_OASPL_S)){sp_OASPL_S=0;}
    if(qIsInf(sp_OASPL_P)){sp_OASPL_P=0;}
    if(qIsInf(sp_OASPL)){sp_OASPL=0;}
    if(qIsInf(sp_dBA)){sp_dBA=0;}
    if(qIsInf(sp_dBB)){sp_dBB=0;}
    if(qIsInf(sp_dBC)){sp_dBC=0;}
    if(qIsInf(sp_LedB)){sp_LedB=0;}
    if(qIsInf(sp_LedBAW)){sp_LedBAW=0;}
    if(qIsInf(sp_LedBBW)){sp_LedBBW=0;}
    if(qIsInf(sp_LedBCW)){sp_LedBCW=0;}

    if(qIsInf(sp_OASPL_alpha)){sp_OASPL_alpha=0;}
    if(qIsInf(splog_OASPL_S)){splog_OASPL_S=0;}
    if(qIsInf(splog_OASPL_P)){splog_OASPL_P=0;}
    if(qIsInf(splog_OASPL)){splog_OASPL=0;}
    if(qIsInf(splog_dBA)){splog_dBA=0;}
    if(qIsInf(splog_dBB)){splog_dBB=0;}
    if(qIsInf(splog_dBC)){splog_dBC=0;}
    if(qIsInf(splog_LedB)){splog_LedB=0;}
    if(qIsInf(splog_LedBAW)){splog_LedBAW=0;}
    if(qIsInf(splog_LedBBW)){splog_LedBBW=0;}
    if(qIsInf(splog_LedBCW)){splog_LedBCW=0;}

    if(qIsInf(sp_OASPL_alpha_a)){sp_OASPL_alpha_a=0;}
    if(qIsInf(sp_OASPL_S_a)){sp_OASPL_S_a=0;}
    if(qIsInf(sp_OASPL_P_a)){sp_OASPL_P_a=0;}
    if(qIsInf(sp_OASPL_a)){sp_OASPL_a=0;}
    if(qIsInf(sp_dBA_a)){sp_dBA_a=0;}
    if(qIsInf(sp_dBB_a)){sp_dBB_a=0;}
    if(qIsInf(sp_dBC_a)){sp_dBC_a=0;}
    if(qIsInf(sp_LedB_a)){sp_LedB_a=0;}
    if(qIsInf(sp_LedBAW_a)){sp_LedBAW_a=0;}
    if(qIsInf(sp_LedBBW_a)){sp_LedBBW_a=0;}
    if(qIsInf(sp_LedBCW_a)){sp_LedBCW_a=0;}

    if(qIsInf(sp_OASPL_alpha_a)){sp_OASPL_alpha_a=0;}
    if(qIsInf(splog_OASPL_S_a)){splog_OASPL_S_a=0;}
    if(qIsInf(splog_OASPL_P_a)){splog_OASPL_P_a=0;}
    if(qIsInf(splog_OASPL_a)){splog_OASPL_a=0;}
    if(qIsInf(splog_dBA_a)){splog_dBA_a=0;}
    if(qIsInf(splog_dBB_a)){splog_dBB_a=0;}
    if(qIsInf(splog_dBC_a)){splog_dBC_a=0;}
    if(qIsInf(splog_LedB_a)){splog_LedB_a=0;}
    if(qIsInf(splog_LedBAW_a)){splog_LedBAW_a=0;}
    if(qIsInf(splog_LedBBW_a)){splog_LedBBW_a=0;}
    if(qIsInf(splog_LedBCW_a)){splog_LedBCW_a=0;}

    stlog_OASPL_alpha_a+=splog_OASPL_alpha_a;
    stlog_OASPL_S_a+=splog_OASPL_S_a;
    stlog_OASPL_P_a+=splog_OASPL_P_a;
    stlog_OASPL_a+=splog_OASPL_a;
    stlog_dBA_a+=splog_dBA_a;
    stlog_dBB_a+=splog_dBB_a;
    stlog_dBC_a+=splog_dBC_a;
    stlog_LedB_a+=splog_LedB_a;
    stlog_LedBAW_a+=splog_LedBAW_a;
    stlog_LedBBW_a+=splog_LedBBW_a;
    stlog_LedBCW_a+=splog_LedBCW_a;

    st_OASPL_alpha_a=10*log10(stlog_OASPL_alpha_a);
    st_OASPL_S_a=10*log10(stlog_OASPL_S_a);
    st_OASPL_P_a=10*log10(stlog_OASPL_P_a);
    st_OASPL_a=10*log10(stlog_OASPL_a);
    st_dBA_a=10*log10(stlog_dBA_a);
    st_dBB_a=10*log10(stlog_dBB_a);
    st_dBC_a=10*log10(stlog_dBC_a);
    st_LedB_a=10*log10(stlog_LedB_a);
    st_LedBAW_a=10*log10(stlog_LedBAW_a);
    st_LedBBW_a=10*log10(stlog_LedBBW_a);
    st_LedBCW_a=10*log10(stlog_LedBAW_a);

    if(qIsInf(st_OASPL_alpha_a)){st_OASPL_alpha_a=0;}
    if(qIsInf(st_OASPL_S_a)){st_OASPL_S_a=0;}
    if(qIsInf(st_OASPL_P_a)){st_OASPL_P_a=0;}
    if(qIsInf(st_OASPL_a)){st_OASPL_a=0;}
    if(qIsInf(st_dBA_a)){st_dBA_a=0;}
    if(qIsInf(st_dBB_a)){st_dBB_a=0;}
    if(qIsInf(st_dBC_a)){st_dBC_a=0;}
    if(qIsInf(st_LedB_a)){st_LedB_a=0;}
    if(qIsInf(st_LedBAW_a)){st_LedBAW_a=0;}
    if(qIsInf(st_LedBBW_a)){st_LedBBW_a=0;}
    if(qIsInf(st_LedBCW_a)){st_LedBCW_a=0;}

    stlog_OASPL_alpha+=splog_OASPL_alpha;
    stlog_OASPL_S+=splog_OASPL_S;
    stlog_OASPL_P+=splog_OASPL_P;
    stlog_OASPL+=splog_OASPL;
    stlog_dBA+=splog_dBA;
    stlog_dBB+=splog_dBB;
    stlog_dBC+=splog_dBC;
    stlog_LedB+=splog_LedB;
    stlog_LedBAW+=splog_LedBAW;
    stlog_LedBBW+=splog_LedBBW;
    stlog_LedBCW+=splog_LedBCW;

    st_OASPL_alpha=10*log10(stlog_OASPL_alpha);
    st_OASPL_S=10*log10(stlog_OASPL_S);
    st_OASPL_P=10*log10(stlog_OASPL_P);
    st_OASPL=10*log10(stlog_OASPL);
    st_dBA=10*log10(stlog_dBA);
    st_dBB=10*log10(stlog_dBB);
    st_dBC=10*log10(stlog_dBC);
    st_LedB=10*log10(stlog_LedB);
    st_LedBAW=10*log10(stlog_LedBAW);
    st_LedBBW=10*log10(stlog_LedBBW);
    st_LedBCW=10*log10(stlog_LedBAW);

    if(qIsInf(st_OASPL_alpha)){st_OASPL_alpha=0;}
    if(qIsInf(st_OASPL_S)){st_OASPL_S=0;}
    if(qIsInf(st_OASPL_P)){st_OASPL_P=0;}
    if(qIsInf(st_OASPL)){st_OASPL=0;}
    if(qIsInf(st_dBA)){st_dBA=0;}
    if(qIsInf(st_dBB)){st_dBB=0;}
    if(qIsInf(st_dBC)){st_dBC=0;}
    if(qIsInf(st_LedB)){st_LedB=0;}
    if(qIsInf(st_LedBAW)){st_LedBAW=0;}
    if(qIsInf(st_LedBBW)){st_LedBBW=0;}
    if(qIsInf(st_LedBCW)){st_LedBCW=0;}

    if (m_parameter.Lowson_type==2)   {
SPL_blade_total_aux+=st_OASPL_alpha+st_OASPL_S+st_OASPL_P+st_OASPL+st_dBA+st_dBB+st_dBC+st_LedB+st_LedBAW+st_LedBBW+st_LedBCW;

SPL_blade_total_aux_a+=st_OASPL_alpha_a+st_OASPL_S_a+st_OASPL_P_a+st_OASPL_a+st_dBA_a+st_dBB_a+st_dBC_a+st_LedB_a+st_LedBAW_a+st_LedBBW_a+st_LedBCW_a;
    }
    else {
SPL_blade_total_aux=st_OASPL_alpha+st_OASPL_S+st_OASPL_P+st_OASPL+st_dBA+st_dBB+st_dBC+st_LedB+st_LedBAW+st_LedBBW;

SPL_blade_total_aux_a+=st_OASPL_alpha_a+st_OASPL_S_a+st_OASPL_P_a+st_OASPL_a+st_dBA_a+st_dBB_a+st_dBC_a+st_LedB_a+st_LedBAW_a+st_LedBBW_a;
            }

if (j==w){
SPL_blade_total = 10.*log10(1./number_of_segments*(SPL_blade_total_aux));

SPL_blade_total_a = 10.*log10(1./number_of_segments*(SPL_blade_total_aux_a));
}

    QString observations_x("");
    if (!((Reynolds[i] >9.5*pow(10,5)) & (Reynolds[i]<2.5*pow(10,6))) || Mach[i]>=0.19 || !((alpha[i]<=0) & (alpha[i]>=0)) || !((alpha[i]<=5) & (alpha[i]>=5)) || !((alpha[i]<=10) & (alpha[i]>=10))){
    observations_x.append("1 ");}
    if((Reynolds[i] <=4.8*pow(10,4)) || (Reynolds[i]>=2.5*pow(10,6)) || Mach[i]<0.208 || !((alpha[i]<=0) & (alpha[i]>=0))){
    observations_x.append("2 ");}
    if (Reynolds[i] >=3*pow(10,6) || Mach[i]<0.208 || alpha[i]<=19.8){
    observations_x.append("3 ");}
    if (Reynolds[i] >1.5*pow(10,6) || Mach[i]<0.208 || alpha[i]<=19.8){
    observations_x.append('4');}

    SPL_LedB_val.clear();
    SPL_LedB_valAW.clear();
    SPL_LedB_valBW.clear();
    SPL_LedB_valCW.clear();
    SPL_alpha_val.clear();
    SPL_dB_S_val.clear();
    SPL_dB_P_val.clear();
    SPL_S_val.clear();
    SPL_P_val.clear();
    SPL_dB_val.clear();
    SPL_A_val.clear();
    SPL_B_val.clear();
    SPL_C_val.clear();

        if(BPM_validation){
    SPL_alpha_val.append(QString::number(SPL_alpha[j], 'f', 2));
    SPL_dB_S_val.append(QString::number(SPL_dB_S[j], 'f', 2));
    SPL_dB_P_val.append(QString::number(SPL_dB_P[j], 'f', 2));
    SPL_dB_val.append(QString::number(SPL_dB_S[j], 'f', 2));
    SPL_A_val.append(QString::number(SPL_A[j], 'f', 2));
    SPL_B_val.append(QString::number(SPL_B[j], 'f', 2));
    SPL_C_val.append(QString::number(SPL_C[j], 'f', 2));

    if(SPL_S[j]==-999999999999.){SPL_S_val.append("N/A");} else {SPL_S_val.append(QString::number(SPL_S[j], 'f', 2));}
    if(SPL_P[j]==-999999999999.){SPL_P_val.append("N/A");} else {SPL_P_val.append(QString::number(SPL_P[j], 'f', 2));}
       }
        else{
    SPL_alpha_val.append("N/A");
    SPL_dB_S_val.append("N/A");
    SPL_dB_P_val.append("N/A");
    SPL_S_val.append("N/A");
    SPL_P_val.append("N/A");
    SPL_dB_val.append("N/A");
    SPL_A_val.append("N/A");
    SPL_B_val.append("N/A");
    SPL_C_val.append("N/A");
        }

    if (LE_validation){
    SPL_LedB_val.append(QString::number(SPL_LedB[j], 'f', 2));
    SPL_LedB_valAW.append(QString::number(SPL_LedBAW[j], 'f', 2));
    SPL_LedB_valBW.append(QString::number(SPL_LedBBW[j], 'f', 2));
    SPL_LedB_valCW.append(QString::number(SPL_LedBCW[j], 'f', 2));
    }
    else{
    SPL_LedB_val.append("N/A");
    SPL_LedB_valAW.append("N/A");
    SPL_LedB_valBW.append("N/A");
    SPL_LedB_valCW.append("N/A");
    }

if (m_parameter.Lowson_type!=0){
        stream << qSetFieldWidth(14)  <<
                  Frequency[j] << ";" <<
//                  Sts[j]<< ";" <<
//                  b_alpha[j]  << ";" <<
//                  B_min[j] << ";" <<
//                  B_max[j] << ";" <<
//                  B_b[j]  << ";" <<
//                  a_alpha[j] << ";" <<
//                  A_min_alpha[j]  << ";" <<
//                  A_max_alpha[j]  << ";" <<
//                  Alin_a[j]  << ";" <<
//                  AWeighting[j]  << ";" <<
//                  BWeighting[j]  << ";" <<
//                  CWeighting[j]  << ";" <<
//                  SPL_alpha_min0[j]  << ";" <<
//                  SPL_alpha_big0[j]  << ";" <<
//                  dBA_alpha_min0[j]  << ";" <<
//                  dBA_alpha_big0[j]  << ";" <<
//                  dBB_alpha_min0[j]  << ";" <<
//                  dBB_alpha_big0[j]  << ";" <<
//                  dBC_alpha_min0[j]  << ";" <<
//                  dBC_alpha_big0[j]  << ";" <<
//                  St1_bar[j]  << ";" <<
//                  a_S[j]  << ";" <<
//                  A_min_S[j]  << ";" <<
//                  A_max_S[j]  << ";" <<
//                  A_a_S[j]  << ";" <<
//                  SPL_dB_S_val  << ";" <<
//                  dBA_S[j]  << ";" <<
//                  dBB_S[j]  << ";" <<
//                  dBC_S[j]  << ";" <<
//                  Sts_St1_bar[j]  << ";" <<
//                  Stp_P[j]  << ";" <<
//                  a_P[j]  << ";" <<
//                  A_min_P[j]  << ";" <<
//                  A_max_P[j]  << ";" <<
//                  A_a_P[j]  << ";" <<
//                  SPL_dB_P_val  << ";" <<
//                  dBA_P[j]  << ";" <<
//                  dBB_P[j]  << ";" <<
//                  dBC_P[j]  << ";" <<
                  SPL_dB_val  << ";" <<
                  SPL_alpha_val  << ";" <<
                  SPL_S_val << ";" <<
                  SPL_P_val  << ";" <<
                  SPL_A_val  << ";" <<
                  SPL_B_val  << ";" <<
                  SPL_C_val  << ";" <<
                  SPL_LedB_val  << ";" <<
                  SPL_LedB_valAW  << ";" <<
                  SPL_LedB_valBW  << ";" <<
                  SPL_LedB_valCW  << ";" <<
//                  slog_SPL_alpha[j]  << ";" <<
//                  slog_SPL_S[j]  << ";" <<
//                  slog_SPL_P[j]  << ";" <<
//                  slog_SPL[j]  << ";" <<
//                  slog_dBA[j]  << ";" <<
//                  slog_dBB[j]  << ";" <<
//                  slog_dBC[j]  << ";" <<
//                  slog_LedB[j]  << ";" <<
                  endl;
}
else
{
        stream << qSetFieldWidth(14)  <<
                  Frequency[j] << ";" <<
//                  Sts[j]<< ";" <<
//                  b_alpha[j]  << ";" <<
//                  B_min[j] << ";" <<
//                  B_max[j] << ";" <<
//                  B_b[j]  << ";" <<
//                  a_alpha[j] << ";" <<
//                  A_min_alpha[j]  << ";" <<
//                  A_max_alpha[j]  << ";" <<
//                  Alin_a[j]  << ";" <<
//                  AWeighting[j]  << ";" <<
//                  BWeighting[j]  << ";" <<
//                  CWeighting[j]  << ";" <<
//                  SPL_alpha_min0[j]  << ";" <<
//                  SPL_alpha_big0[j]  << ";" <<
//                  dBA_alpha_min0[j]  << ";" <<
//                  dBA_alpha_big0[j]  << ";" <<
//                  dBB_alpha_min0[j]  << ";" <<
//                  dBB_alpha_big0[j]  << ";" <<
//                  dBC_alpha_min0[j]  << ";" <<
//                  dBC_alpha_big0[j]  << ";" <<
//                  St1_bar[j]  << ";" <<
//                  a_S[j]  << ";" <<
//                  A_min_S[j]  << ";" <<
//                  A_max_S[j]  << ";" <<
//                  A_a_S[j]  << ";" <<
//                  SPL_dB_S_val  << ";" <<
//                  dBA_S[j]  << ";" <<
//                  dBB_S[j]  << ";" <<
//                  dBC_S[j]  << ";" <<
//                  Sts_St1_bar[j]  << ";" <<
//                  Stp_P[j]  << ";" <<
//                  a_P[j]  << ";" <<
//                  A_min_P[j]  << ";" <<
//                  A_max_P[j]  << ";" <<
//                  A_a_P[j]  << ";" <<
//                  SPL_dB_P_val  << ";" <<
//                  dBA_P[j]  << ";" <<
//                  dBB_P[j]  << ";" <<
//                  dBC_P[j]  << ";" <<
                  SPL_dB_val  << ";" <<
                  SPL_alpha_val  << ";" <<
                  SPL_S_val << ";" <<
                  SPL_P_val  << ";" <<
                  SPL_A_val  << ";" <<
                  SPL_B_val  << ";" <<
                  SPL_C_val  << ";" <<
//                  slog_SPL_alpha[j]  << ";" <<
//                  slog_SPL_S[j]  << ";" <<
//                  slog_SPL_P[j]  << ";" <<
//                  slog_SPL[j]  << ";" <<
//                  slog_dBA[j]  << ";" <<
//                  slog_dBB[j]  << ";" <<
//                  slog_dBC[j]  << ";" <<
                  endl;
}

if(j==(w-1)){
        stream << endl;
        if (sp_OASPL<=0 & sp_OASPL>=0){stream << "OASPL: "  << "N/A"<< endl;} else {stream << "SPL: "  << sp_OASPL<< "dB" << endl;}
        if (sp_dBA<=0 & sp_dBA>=0){stream << "OASPL(A): "  << "N/A"<< endl;} else {stream << "OASPL(A): "  << sp_dBA<< "dB" << endl;}
        if (sp_dBB<=0 & sp_dBB>=0){stream << "OASPL(B): "  << "N/A"<< endl;} else {stream << "OASPL(B): "  << sp_dBB<< "dB" << endl;}
        if (sp_dBC<=0 & sp_dBC>=0){stream << "OASPL(C): "  << "N/A"<< endl;} else {stream << "OASPL(C): "  << sp_dBC<< "dB" << endl;}
        if (sp_OASPL_alpha<=0 & sp_OASPL_alpha>=0){stream << "SPL_a: "  << "N/A"<< endl;} else {stream << "SPL_a: "  << sp_OASPL_alpha<< endl;}
        if (sp_OASPL_S<=0 & sp_OASPL_S>=0){stream << "SPL_s: "  << "N/A"<< endl;} else {stream << "SPL_s: "  << sp_OASPL_S<< endl;}
         if (sp_OASPL_P<=0 & sp_OASPL_P>=0){stream << "SPL_p: "  << "N/A"<< endl;} else {stream << "SPL_p: "  << sp_OASPL_P<< endl;}
          if (m_parameter.Lowson_type!=0){
        if (sp_LedB<=0 & sp_LedB>=0){stream << "SPL_LE: "  << "N/A"<< endl;} else {stream << "SPL_LE: "  << sp_LedB<< "dB" << endl;}
        if (sp_LedBAW<=0 & sp_LedBAW>=0){stream << "SPL_LE(A): "  << "N/A"<< endl;} else {stream << "SPL_LE(A): "  << sp_LedBAW<< "dB" << endl;}
        if (sp_LedBBW<=0 & sp_LedBAW>=0){stream << "SPL_LE(A): "  << "N/A"<< endl;} else {stream << "SPL_LE(A): "  << sp_LedBBW<< "dB" << endl;}
        if (sp_LedBCW<=0 & sp_LedBAW>=0){stream << "SPL_LE(A): "  << "N/A"<< endl;} else {stream << "SPL_LE(A): "  << sp_LedBCW<< "dB" << endl;}
          }
                stream << endl;
if (i>=(number_of_segments-1)){
stream << "***********************************************************"<<endl;
stream << "Total Blade: "  << endl;
stream << "OASPL: "  << st_OASPL_a<< "dB" << endl;
stream << "OASPL (A): "  << st_dBA_a<< "dB" <<  endl;
stream << "OASPL (B): "  << st_dBB_a<< "dB" <<  endl;
stream << "OASPL (C): "  << st_dBC_a<< "dB" <<  endl;
stream << "SPL_a: "  << st_OASPL_alpha_a<< endl;
stream << "SPL_s: "  << st_OASPL_S_a<< endl;
stream << "SPL_p: "  << st_OASPL_P_a<< endl;
if (m_parameter.Lowson_type!=0){
stream << "SPL_LE: "  << st_LedB_a<< "dB" <<  endl;
stream << "SPL_LE (A): "  << st_LedBAW_a<< "dB" <<  endl;
stream << "SPL_LE (B): "  << st_LedBBW_a<< "dB" <<  endl;
stream << "SPL_LE (C): "  << st_LedBCW_a<< "dB" <<  endl;}
stream << "SPL_Blade_total: " << SPL_blade_total_a << "dB" <<  endl;
stream << "***********************************************************"<<endl;
stream << endl;}

if(z>=lend & i>=(number_of_segments-1)){
stream << "***********************************************************"<<endl;
stream << "TOTAL COMPLETE" << endl;
stream <<"***********************************************************"<<endl;
                stream << "Total: "  << endl;
                stream << "OASPL: "  << st_OASPL<< "dB" << endl;
                stream << "OASPL (A): "  << st_dBA<< "dB" <<  endl;
                stream << "OASPL (B): "  << st_dBB<< "dB" <<  endl;
                stream << "OASPL (C): "  << st_dBC<< "dB" <<  endl;
                stream << "SPL_a: "  << st_OASPL_alpha<< endl;
                stream << "SPL_s: "  << st_OASPL_S<< endl;
                stream << "SPL_p: "  << st_OASPL_P<< endl;
          if (m_parameter.Lowson_type!=0){
                stream << "SPL_LE: "  << st_LedB<< "dB" <<  endl;
          stream << "SPL_LE (A): "  << st_LedBAW<< "dB" <<  endl;
        stream << "SPL_LE (B): "  << st_LedBBW<< "dB" <<  endl;
stream << "SPL_LE (C): "  << st_LedBCW<< "dB" <<  endl;}
                stream << "SPL_Blade_total: " << SPL_blade_total << "dB" <<  endl;
                stream << "***********************************************************"<<endl;
stream << endl;
}
}}}}
z=z+ldelta;
        }}
//Sara

void NoiseSimulation::setAnalyzedOpPoints(QVector<OpPoint *> newList) {
    removeAllParents();
    for (int i = 0; i < newList.size(); ++i) {
        addParent(newList[i]);
    }

    m_parameter.analyzedOpPoints = newList;
}

QVariant NoiseSimulation::accessParameter(Parameter::NoiseSimulation::Key key, QVariant value) {
//Sara
    QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
    double outer_radius=pbem->m_pTData->OuterRadius;
//Sara

    typedef Parameter::NoiseSimulation P;

    const bool set = value.isValid();
    switch (key) {
    case P::Name: if(set) setName(value.toString()); else value = getName(); break;
    case P::WettedLength:
        if(set) m_parameter.wettedLength = value.toDouble();
        else value = m_parameter.wettedLength; break;
    case P::DistanceObsever:
        if(set) m_parameter.distanceObsever = value.toDouble();
        else value = m_parameter.distanceObsever; break;
    case P::OriginalVelocity:
        if(set) m_parameter.originalVelocity = value.toDouble();
        else value = m_parameter.originalVelocity; break;
    case P::OriginalChordLength:
        if(set) m_parameter.originalChordLength = value.toDouble();
        else value = m_parameter.originalChordLength; break;
    case P::OriginalMach:
        if(set) m_parameter.originalMach = value.toDouble();
        else value = m_parameter.originalMach; break;
    case P::DStarChordStation:
        if(set) m_parameter.dStarChordStation = value.toDouble();
        else value = m_parameter.dStarChordStation; break;
    case P::DStarScalingFactor:
        if(set) m_parameter.dStarScalingFactor = value.toDouble();
        else value = m_parameter.dStarScalingFactor; break;
    case P::EddyConvectionMach:
        if(set) m_parameter.eddyConvectionMach = value.toDouble();
        else value = m_parameter.eddyConvectionMach; break;
    case P::DirectivityTheta:
        if(set) m_parameter.directivityGreek = value.toDouble();
        else value = m_parameter.directivityGreek; break;
    case P::DirectivityPhi:
        if(set) m_parameter.directivityPhi = value.toDouble();
        else value = m_parameter.directivityPhi; break;
    case P::SeparatedFlow:
        if(set) m_parameter.separatedFlow = value.toBool();
        else value = m_parameter.separatedFlow; break;
    case P::SuctionSide:
        if(set) m_parameter.suctionSide = value.toBool();
        else value = m_parameter.suctionSide; break;
    case P::PressureSide:
        if(set) m_parameter.pressureSide = value.toBool();
        else value = m_parameter.pressureSide; break;
    case P::Aoa:
        if(set) m_parameter.aoa = value.toDouble();
        else value = m_parameter.aoa; break;
    case P::ChordBasedReynolds:
        if(set) m_parameter.chordBasedReynolds = value.toDouble();
        else value = m_parameter.chordBasedReynolds; break;
    case P::Transition:
        if(set) m_parameter.transition = static_cast<NoiseParameter::TransitionType>(value.toInt());
        else value = static_cast<int>(m_parameter.transition); break;
        //Alexandre MOD
        case P::IntegralLengthScale:
            if(set) m_parameter.IntegralLengthScale = value.toDouble();
            else value = m_parameter.IntegralLengthScale; break;
        case P::TurbulenceIntensity:
        if(set) m_parameter.TurbulenceIntensity = value.toDouble();
            else value = m_parameter.TurbulenceIntensity; break;

        //Sara
    case P::sects:
        if(set) m_parameter.sects=value.toDouble();
else {m_parameter.sects = pbem->dlg_elements; value=m_parameter.sects;}
        break;

    case P::TSRtd:
        if(set) m_parameter.TSRtd = value.toDouble();
        else {value=m_parameter.TSRtd;}
break;

    case P::u_wind_speed:
        if(set) m_parameter.u_wind_speed = value.toDouble();
        else {value=m_parameter.u_wind_speed;}
break;

    case P::rot_speed:
        if(set) {m_parameter.rot_speed = value.toDouble();}
        else{value=m_parameter.rot_speed;}
        break;

    case P::TSR_check:
        if(set) m_parameter.TSR_check = value.toBool();
        else {value=m_parameter.TSR_check;}
break;

    case P::u_wind_speed_check:
        if(set) m_parameter.u_wind_speed_check = value.toBool();
        else {value=m_parameter.u_wind_speed_check;}
break;

    case P::rot_speed_check:
        if(set) {m_parameter.rot_speed_check = value.toBool();}
else {value=m_parameter.rot_speed_check;}
        break;

    case P::dstar_type:
        if(set) m_parameter.dstar_type = value.toInt();
        else {value = m_parameter.dstar_type;}break;

    case P::phi_type:
        if(set) m_parameter.phi_type = value.toInt();
        else {value = m_parameter.phi_type;}break;

    case P::obs_x_pos:
//        if (m_parameter.obs_x_pos==1){m_parameter.obs_x_pos=10;}
        if(set) m_parameter.obs_x_pos = value.toDouble();
        else {value = m_parameter.obs_x_pos;}break;

    case P::obs_y_pos:
//        if (m_parameter.obs_y_pos==1){m_parameter.obs_y_pos=10;}
        if(set) m_parameter.obs_y_pos = value.toDouble();
        else {value = m_parameter.obs_y_pos;}break;

    case P::obs_z_pos:
//        if (m_parameter.obs_z_pos==1){
//            double hub_radius=pbem->m_pBlade->m_HubRadius;
//            double blade_radius=(outer_radius-hub_radius);
//            m_parameter.obs_z_pos=blade_radius/2.;
//        }
        if(set) m_parameter.obs_z_pos = value.toDouble();
        else {value = m_parameter.obs_z_pos;}break;

    case P::Lowson_type:
        if(set) {m_parameter.Lowson_type = value.toInt();}
else {value=m_parameter.Lowson_type;}
        break;

    }
// Sara

SimuWidget *pSimuWidget = (SimuWidget *) pbem->m_pSimuWidget;

//cálculos TSR w e u
m_parameter.TSR_calc=2.*PI*m_parameter.rot_speed/60.*outer_radius/m_parameter.u_wind_speed;

m_parameter.rot_speed_calc=m_parameter.TSRtd*m_parameter.u_wind_speed*60./(2.*PI*outer_radius);

m_parameter.u_wind_speed_calc=2.*PI*m_parameter.rot_speed/60.*outer_radius/m_parameter.TSRtd;

//cálculo para não sets
if(m_parameter.u_wind_speed_check==false){
m_parameter.u_wind_speed=m_parameter.u_wind_speed_calc;}

if(m_parameter.TSR_check==false){
SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
double lend  =   pSimuWidget->m_pctrlLELineEdit->getValue();

double v=lstart;
while (qAbs(m_parameter.TSR_calc-v)>ldelta){
v=v+ldelta;
}

int q=0;
if(m_parameter.TSR_calc<v || m_parameter.TSR_calc>v){
m_parameter.TSR_calc=v;
m_parameter.TSR_check=true;
m_parameter.rot_speed=m_parameter.rot_speed_calc;
m_parameter.rot_speed_check=false;
m_parameter.u_wind_speed_check=true;
}
m_parameter.TSRtd=m_parameter.TSR_calc;
}

if(m_parameter.rot_speed_check==false){m_parameter.rot_speed=m_parameter.rot_speed_calc;}

//condição inicial
if((m_parameter.u_wind_speed_check==false) && (m_parameter.rot_speed_check==false) && (m_parameter.TSR_check==false)){
 m_parameter.u_wind_speed_check=true;
 m_parameter.rot_speed_check=false;
 m_parameter.TSR_check=true;
 m_parameter.TSRtd=7;
 m_parameter.u_wind_speed=pSimuWidget->m_pctrlWindspeed->getValue();
 m_parameter.rot_speed=m_parameter.rot_speed_calc;
}

// Sara

    return (set ? QVariant() : value);
}
