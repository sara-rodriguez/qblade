#include "NoiseSimulation.h"
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
#include "NoiseException.h"
#include "NoiseOpPoint.h"
#include <cmath>
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
        default: return NULL;
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
              << "SPL (dB(C))";

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
        stream << endl;

        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" <<
                  "SPL (dB)" <<
                  "SPLa" <<
                  "SPLs" <<
                  "SPLp" <<
                  "SPL (dB(A))" <<
                  "SPL (dB(B))" <<
                  "SPL (dB(C))" << endl;

        for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {
            stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] <<
                      m_calculation.SPLdB()[i][j] <<
                      m_calculation.SPLadB()[i][j] <<
                      m_calculation.SPLsdB()[i][j] <<
                      m_calculation.SPLpdB()[i][j] <<
                      m_calculation.SPLdBAW()[i][j] <<
                      m_calculation.SPLdBBW()[i][j] <<
                      m_calculation.SPLdBCW()[i][j] << endl;
        }

        stream << endl;
        stream << endl;
    }
    qDeleteAll(noiseOpPoints);
}

//Sara
void NoiseSimulation::export3DCalculation(QTextStream &stream) {
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

    const double CWeighting[] = { -4.4,  -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.2,
                                                                        -0.1,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                                                                         0.0,   0.0,  -0.1,  -0.2,  -0.3,  -0.5,  -0.8,  -1.3,
                                                                        -2.0,  -3.0,  -4.4,  -6.2,  -8.5, -11.2};

    const double Frequency[]= {25, 31.5, 40, 50, 63, 80, 100, 125, 160,200, 250, 315, 400,500, 630, 800, 1000,1250, 1600, 2000, 2500, 3150, 4000, 5000,6300, 8000, 10000, 12500, 16000, 20000};

    stream << "3D Noise prediction file export" << endl;
    stream << endl;

QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();


SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
double z=lstart;
double approaxing_wind_speed = m_parameter.originalVelocity;

//nome da pa
//QGroupBox* SimulationCreatorDialog<ParameterGroup>::constructParameterBox(QString defaultName)
//SimulationCreatorDialog *pSimulationCreatorDialog = (SimulationCreatorDialog *) constructParameterBox->;
//SimulationCreatorDialog->ParameterGroup->constructParameterBox->defaultName

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    foreach(BData * bdata, pBEM->m_pBEMData->GetBData()){
    int number_of_segments = bdata->m_pos.size();
    double rho = pBEM->dlg_rho;
    double dynamic_visc = pBEM->dlg_visc;
    double cin_visc = dynamic_visc/rho;
    double K_air = 1.4;
    double R_air = 286.9;
    double T_std_cond = 288.15;
    double P_std_cond = 101300;
    double lambda = pBEM->dlg_lambda;
    int mpos_size = bdata->m_pos.size(); //total number of segments
    double finalradius = bdata->m_pos.value(mpos_size-1);
    double nom_tg_speed = bdata->windspeed*lambda;
    double omega = nom_tg_speed/finalradius;
    double rotation = 60/(M_PI*100/nom_tg_speed);

//    qDebug() << "tamanho de alpha: " << bdata->m_alpha.size();
//    qDebug() << "tamanho de cl/cd: " << bdata->m_LD.size();
//    qDebug() << "segmentos: " << number_of_segments;

    QString str= QString::number(z, 'f', 1);

    stream << "Tip Speed Ratio: " << str << endl;
    stream << endl;

    stream << qSetFieldWidth(14)  <<
              "Sect"  << ";" <<
              "Radius [m]"  << ";" <<
              "r/R"  << ";" <<
              "Chord [m]" << ";" <<
              "Theta [deg]" << ";" <<
              "Phi [deg]" << ";" <<
              "Axial Ind. Fact. (a)"  << ";" <<
              "Axial Velocity [m/s]" << ";" <<
              "Tg. Ind. Fact. [a']"  << ";" <<
              "Tg. Speed [m/s]" << ";" <<
//              "Res. Local Speed Calc [m/s]" << ";" <<
              "Res. Local Speed BEM [m/s]" << ";" <<
              "Re BEM"  << ";" <<
//              "Re calc"  << ";" <<
//              "Mach BEM"  << ";" <<
              "Mach calc"  << ";" <<
              "(cl/cd) max"   << ";" <<
              "cl"   << ";" <<
              "cd"   << ";" <<
              "c/R"    <<";"   <<
              "Alpha [deg]"   <<";"   <<
              "D* nat.trans. S"   <<";"   <<
              "D* nat.trans. P" <<";"   <<
              "D* heavy trip. S"<<";"   <<
              "D* heavy trip. P"   <<";"   <<
              endl;

    //definitions
    double axial_ind_fact[number_of_segments];
    double axial_ind_fact_n[number_of_segments];
    double axial_velocity[number_of_segments];
    double tangential_speed[number_of_segments];
    double resultant_local_speed[number_of_segments];
    double chord[number_of_segments];
//    double Reynolds_calc[number_of_segments];
    double Reynolds[number_of_segments];
    double Mach_calc[number_of_segments];
    double alpha[number_of_segments];
    double phi[number_of_segments];
    double theta[number_of_segments];
    double cl_cd[number_of_segments];
    double r_R[number_of_segments];
    double c_Rx[number_of_segments];
    double D_starred_C_HT[number_of_segments];
    double D_starred_HT[number_of_segments];
    double D_starred_C_N[number_of_segments];
    double D_starred_N[number_of_segments];

    double r_R0  =  0.05; double c_R0 = 0.05500;
    double r_R1  =  0.25; double c_R1 = 0.07500;
    double r_R2  =  1.00; double c_R2 = 0.02000;

        for (int i = 0; i < number_of_segments; ++i) {

            // definitions
            axial_ind_fact[i] = bdata->m_a_axial.value(i);
            axial_ind_fact_n[i] = bdata->m_a_axial.value(i+1);

            if (i<number_of_segments/2) {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact[i]);}
            else {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact_n[i]);}

            tangential_speed[i] = omega*bdata->m_pos.value(i)*(1.f+bdata->m_a_tangential.value(i));
            resultant_local_speed[i] = qSqrt(pow(axial_velocity[i],2)+pow(tangential_speed[i],2));
            chord[i] = bdata->m_c_local.value(i);
//            Reynolds_calc[i] = resultant_local_speed[i]*chord[i]/dynamic_visc;//*rho
            Reynolds[i] = bdata->m_Reynolds.value(i);
            Mach_calc[i] = resultant_local_speed[i]/sqrt(R_air*K_air*T_std_cond);
            alpha[i] = bdata->m_alpha.value(i);
            phi[i] = bdata->m_phi.value(i);
            theta[i] = bdata->m_theta.value(i);
            cl_cd[i] =  bdata->m_LD.value(i);
            r_R[i] = bdata->m_pos.value(i)/finalradius;

            if (r_R[i] <= r_R0) {c_Rx[i] = c_R0;}
            if (r_R[i] > r_R0 && r_R[i] < r_R1) {c_Rx[i] = (r_R[i]-r_R0)*(c_R1-c_R0)/(0.25-r_R0)+c_R0;}
            if (r_R[i] <= r_R1 && r_R[i] >= r_R1) {c_Rx[i] = c_R1;}
            if (r_R[i] > r_R1 && r_R[i] < r_R2) {c_Rx[i] = (r_R[i]-r_R1)*(c_R2-c_R1)/(r_R2-r_R1)+c_R1;}
            if (r_R[i] >= r_R2) {c_Rx[i] = c_R2;}

            QString c_R= QString::number(c_Rx[i], 'f', 5);
            double Mach[number_of_segments];
            Mach[i]=Mach_calc[i];

            double corr_fact[number_of_segments];

//heavy tripping
if (Reynolds[i]>300000){
    D_starred_C_HT[i]=pow(10,(3.411-1.5397*log10(Reynolds[i])+0.1059*pow(log10(Reynolds[i]),2)));
}
else {D_starred_C_HT[i]=0.0601*(pow(Reynolds[i],(-0.114)));}

D_starred_HT[i]=chord[i]*D_starred_C_HT[i];

//natural transition
    D_starred_C_N[i]=pow(10,(3.0187-1.5397*log10(Reynolds[i])+0.1059*pow(log10(Reynolds[i]),2)));

D_starred_N[i]=D_starred_C_N[i]*alpha[i];

double D_starred_HT_S[number_of_segments];
double D_starred_HT_P[number_of_segments];
double D_starred_N_S[number_of_segments];
double D_starred_N_P[number_of_segments];

D_starred_HT_S[i]=D_starred_HT[i];
D_starred_HT_P[i]=D_starred_HT[i];
D_starred_N_S[i]=D_starred_N[i];
D_starred_N_P[i]=D_starred_N[i];


//alpha !=0 pressure side
if (alpha[i]!=0){
corr_fact[i]=pow(10,(-0.0432*alpha[i]+0.00113*pow(alpha[i],2)));
D_starred_HT_P[i]=D_starred_HT[i]*corr_fact[i];
D_starred_N_P[i]=D_starred_N[i]*corr_fact[i];
}

//alpha !=0 suction side heavy tripping
if (alpha[i]>0 & alpha[i]<=5){
corr_fact[i]=pow(10,(0.0679*alpha[i]));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

if (alpha[i]>5 & alpha[i]<=12.5){
corr_fact[i]=0.381*(pow(10,(0.1516*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

if (alpha[i]>12.5 & alpha[i]<=25){
corr_fact[i]=14.296*(pow(10,(0.0258*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
}

//alpha !=0 suction side natural transition
if (alpha[i]>0 & alpha[i]<=7.5){
corr_fact[i]=pow(10,(0.0679*alpha[i]));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

if (alpha[i]>7.5 & alpha[i]<=12.5){
corr_fact[i]=0.0162*(pow(10,(0.3066*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

if (alpha[i]>12.5 & alpha[i]<=25){
corr_fact[i]=54.42*(pow(10,(0.0258*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
}

//teste inicio

//Rodar de novo o XFoil com os Reynolds e Machs definidos em Mach_calc[i] e bdata->m_Reynolds.value(i)

//Fazer um loop externo e extrair os dados de i=1 a i=40 para cada corda de:  m_DStarInterpolatedS e m_DStarInterpolatedP e dentro com os pontos operacionais de todas as frequencias de 25 a 20000Hz (como era feito antes).

//

//            D_starred[i]=chord[i]*D_starred_C[i];
//            //D_starred[i]=m_opPoint->topDStar.second[i];

//            double m_DStarInterpolatedS[number_of_segments];
//            double m_DStarInterpolatedP[number_of_segments];

//    bool dStarOrder[number_of_segments];
//    for (i=1;i<(number_of_segments+1);i++){
//        if (alpha[i]<0){dStarOrder[i]= true;} else {dStarOrder[i]= false;}
//    }

    //        m_DStarInterpolatedS = getDStarInterpolated(dStarOrder,nop);

//            double m_DStarFinalS[number_of_segments];
//            double m_DStarFinalP[number_of_segments];
    //        m_DStarFinalS = m_DStarInterpolatedS * m_parameter->originalChordLength * m_parameter->dStarScalingFactor;

    //        m_DStarFinalS = m_DStarInterpolatedS * m_parameter->originalChordLength * m_parameter->dStarScalingFactor;
    //        m_DStarFinalP = m_DStarInterpolatedP * m_parameter->originalChordLength * m_parameter->dStarScalingFactor;

// teste fim

        stream << qSetFieldWidth(14)  <<
                      (i+1) << ";" <<
                      bdata->m_pos.value(i) << ";" <<
                      r_R[i] << ";" <<
                      chord[i] << ";" <<
                      bdata->m_theta.value(i) << ";" <<
                      bdata->m_phi.value(i) << ";" <<
                      axial_ind_fact[i] << ";" <<
                      axial_velocity[i] << ";" <<
                      bdata->m_a_tangential.value(i) << ";" <<
                      tangential_speed[i] << ";" <<
//                      resultant_local_speed[i] << ";" <<
                      bdata->m_Windspeed.value(i) << ";" <<
                      bdata->m_Reynolds.value(i) << ";" <<
//                      Reynolds_calc[i] << ";" <<
//                      Mach[i] << ";" <<
                      Mach_calc[i] << ";" <<
                      cl_cd[i] << ";" <<
                      bdata->m_CL.value(i) << ";" <<
                      bdata->m_CD.value(i) << ";" <<
                      c_R  << ";" <<
                      bdata->m_alpha.value(i) << ";" <<
                      D_starred_N_S[i] << ";" <<
                      D_starred_N_P[i] << ";" <<
                      D_starred_HT_S[i] << ";" <<
                      D_starred_HT_P[i] << ";" <<
                      endl;
}
z=z+ldelta;
    }

        //SPL_Alpha ********************************************************

    double r_R0  =  0.05; double c_R0 = 0.05500;
    double r_R1  =  0.25; double c_R1 = 0.07500;
    double r_R2  =  1.00; double c_R2 = 0.02000;

    z=lstart;
    foreach(BData * bdata, pBEM->m_pBEMData->GetBData()){
    int number_of_segments = bdata->m_pos.size();

    double axial_ind_fact[number_of_segments];
    double axial_ind_fact_n[number_of_segments];
    double axial_velocity[number_of_segments];
    double tangential_speed[number_of_segments];
    double resultant_local_speed[number_of_segments];
    double chord[number_of_segments];
    double Reynolds_calc[number_of_segments];
    double Mach_calc[number_of_segments];
    double alpha[number_of_segments];
    double phi[number_of_segments];
    double theta[number_of_segments];
    double cl_cd[number_of_segments];
    double r_R[number_of_segments];
    double c_Rx[number_of_segments];

    double rho = pBEM->dlg_rho;
    double dynamic_visc = pBEM->dlg_visc;
    double cin_visc = dynamic_visc/rho;
    double K_air = 1.4;
    double R_air = 286.9;
    double T_std_cond = 288.15;
    double P_std_cond = 101300;
    double lambda = pBEM->dlg_lambda;
    int mpos_size = bdata->m_pos.size(); //total number of segments
    double finalradius = bdata->m_pos.value(mpos_size-1);
    double nom_tg_speed = bdata->windspeed*lambda;
    double omega = nom_tg_speed/finalradius;
    double rotation = 60/(M_PI*100/nom_tg_speed);
        for (int i = 0; i < number_of_segments; ++i) {

            // definitions
            axial_ind_fact[i] = bdata->m_a_axial.value(i);
            axial_ind_fact_n[i] = bdata->m_a_axial.value(i+1);

            if (i<number_of_segments/2) {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact[i]);}
            else {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact_n[i]);}

            tangential_speed[i] = omega*bdata->m_pos.value(i)*(1.f+bdata->m_a_tangential.value(i));
            resultant_local_speed[i] = qSqrt(pow(axial_velocity[i],2)+pow(tangential_speed[i],2));
            chord[i] = bdata->m_c_local.value(i);
            Reynolds_calc[i] = resultant_local_speed[i]*chord[i]/dynamic_visc;//*rho
            Mach_calc[i] = resultant_local_speed[i]/sqrt(R_air*K_air*T_std_cond);
//            Mach_calc[i] = resultant_local_speed[i]/1235;
            alpha[i] = bdata->m_alpha.value(i);
            phi[i] = bdata->m_phi.value(i);
            theta[i] = bdata->m_theta.value(i);
            cl_cd[i] =  bdata->m_LD.value(i);
            r_R[i] = bdata->m_pos.value(i)/finalradius;

            if (r_R[i] <= r_R0) {c_Rx[i] = c_R0;}
            if (r_R[i] > r_R0 && r_R[i] < r_R1) {c_Rx[i] = (r_R[i]-r_R0)*(c_R1-c_R0)/(0.25-r_R0)+c_R0;}
            if (r_R[i] <= r_R1 && r_R[i] >= r_R1) {c_Rx[i] = c_R1;}
            if (r_R[i] > r_R1 && r_R[i] < r_R2) {c_Rx[i] = (r_R[i]-r_R1)*(c_R2-c_R1)/(r_R2-r_R1)+c_R1;}
            if (r_R[i] >= r_R2) {c_Rx[i] = c_R2;}

            QString c_R= QString::number(c_Rx[i], 'f', 5);
            double Mach[number_of_segments];
            Mach[i]=Mach_calc[i];

        //Calculate the Switching Angle
        double SwAlpha[number_of_segments];
        double SwAlpha_1[number_of_segments];

        SwAlpha_1[i]=23.43*Mach[i]+4.651;
        double SwAlpha_2=12.5;

        if (SwAlpha_1[i]<SwAlpha_2){SwAlpha[i]=SwAlpha_1[i];}
        else {SwAlpha[i]=SwAlpha_2;}

        //Length of Wetted  Trailing Edge
        double L[number_of_segments];
        L[i]=bdata->m_pos.value(i+1)-bdata->m_pos.value(i); //Default = 1 m; 0.4572 m for comparison to BPM data (see page 5 of report).

        double Dh[number_of_segments];
        double Dl[number_of_segments];
        double observer_position = 1.22;

        double EddyMach = m_parameter.eddyConvectionMach;

Dh[i]=(2*pow(sin(qDegreesToRadians(theta[i]/2)),2)*pow(sin(qDegreesToRadians(phi[i])),2))/pow(1+Mach[i]*cos(qDegreesToRadians(theta[i]))*(1+(Mach[i]-Mach[i]*EddyMach)*cos(qDegreesToRadians(phi[i]))),2);

Dl[i]=(2*pow(sin(qDegreesToRadians(theta[i])),2)*pow(sin(qDegreesToRadians(phi[i])),2))/pow((1+(Mach[i]*cos(qDegreesToRadians(theta[i])))),4);

        //For Heavy Tripping
        double D_starred_C[number_of_segments];
        if (bdata->m_Reynolds.value(i)>30000) {D_starred_C[i] =pow(10,(3.411-1.5397*log10(bdata->m_Reynolds.value(i))))+0.1059*pow(log10(bdata->m_Reynolds.value(i)),2);}
        else {D_starred_C[i]=0.0601*pow(bdata->m_Reynolds.value(i),-0.114);}

        //D*
        double D_starred[number_of_segments];
        D_starred[i]=chord[i]*D_starred_C[i];
        //D_starred[i]=m_opPoint->topDStar.second[i];

        double m_DStarInterpolatedS[number_of_segments];
        double m_DStarInterpolatedP[number_of_segments];

bool dStarOrder[number_of_segments];
for (i=1;i<(number_of_segments+1);i++){
    if (alpha[i]<0){dStarOrder[i]= true;} else {dStarOrder[i]= false;}
}

//        m_DStarInterpolatedS = getDStarInterpolated(dStarOrder,nop);

        double m_DStarFinalS[number_of_segments];
        double m_DStarFinalP[number_of_segments];
//        m_DStarFinalS = m_DStarInterpolatedS * m_parameter->originalChordLength * m_parameter->dStarScalingFactor;

//        m_DStarFinalS = m_DStarInterpolatedS * m_parameter->originalChordLength * m_parameter->dStarScalingFactor;
//        m_DStarFinalP = m_DStarInterpolatedP * m_parameter->originalChordLength * m_parameter->dStarScalingFactor;
        //parei aqui

        double first_term[number_of_segments];
        first_term[i]=10*log10(pow(Mach[i],5)*L[i]*Dh[i]*D_starred[i]/pow(observer_position,2));

        double St1[number_of_segments];
        St1[i] = 0.02*(pow(Mach[i],-0.6));

        double St2[number_of_segments];
        if (alpha[i]<1.33){St2[i]=St1[i];}
        else if (alpha[i]>12.5) {St2[i]=St1[i]*4.72;}
        else {St2[i]=St1[i]*pow(10,(0.0054*pow((alpha[i]-1.33),2)));}

        double gamma[number_of_segments];
        gamma[i]=27.094*Mach[i]+3.32;

        double gamma0[number_of_segments];
        gamma0[i]=SwAlpha_1[i];

        double beta[number_of_segments];
        beta[i]=72.65*Mach[i]+10.74;

        double beta0[number_of_segments];
        beta0[i]=-34.19*Mach[i]-13.82;

        double gamma0_gamma_min[number_of_segments];
        gamma0_gamma_min[i]=gamma0[i]-gamma[i];

        double gamma0_gamma_plus[number_of_segments];
        gamma0_gamma_plus[i]=gamma0[i]+gamma[i];

        double K1[number_of_segments];
        if (bdata->m_Reynolds.value(i)<247000)
        {K1[i]=-4.31*log10(bdata->m_Reynolds.value(i))+156.3;}
        else if (bdata->m_Reynolds.value(i)>800000)
        {K1[i]=128.5;}
        else {K1[i]=-9*log10(bdata->m_Reynolds.value(i))+181.6;}

        double K2[number_of_segments];
        if (alpha[i]<gamma0_gamma_min[i])
        {K2[i]=K1[i]-1000;}
        else if (alpha[i]>gamma0_gamma_plus[i])
        {K2[i]=K1[i]-12;}
        else
        {K2[i]=K1[i]+(sqrt(pow(beta[i],2)-pow((beta[i]/gamma[i]),2)*pow((alpha[i]-gamma0[i]),2)))+beta0[i];}

int q_aux=0;
int y_aux=0;

        if (alpha[i]<=SwAlpha[i])
        // if alpha is less than or equal to switching angle
        {
q_aux=q_aux+1;
        double b0[number_of_segments];
        if (bdata->m_Reynolds.value(i)<95200)
        {b0[i]= 0.3;}
        else if (bdata->m_Reynolds.value(i)>857000)
        {b0[i]= 0.56;}
        else {b0[i]=-4.48*pow(10,-13)*(pow((bdata->m_Reynolds.value(i)-857000),2)+0.56);}

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

if (q_aux==1){
        QString str= QString::number(z, 'f', 1);

            stream << "SPL Alpha: " <<  endl;
            stream << "Angles less than the switching angle: "  << endl;
            stream << "Tip Speed Ratio: " << str << endl;
            stream << "Section: " << (i+1)<<"/"<<number_of_segments << endl;
            stream << "D*: " << D_starred[i] << endl;
            stream << endl;

        stream << qSetFieldWidth(14)  <<
                  "Freq[Hz]"  << ";" <<
                  "Sts"  << ";" <<
                  "b"  << ";" <<
                  "B_min(b)" << ";" <<
                  "B_max(b)" << ";" <<
                  "B(b)"  << ";" <<
                  "SPL [dB]" << ";" <<
                  "A-Weighting"  << ";" <<
                  "dB(A)" << ";" <<
                  "B-Weighting" << ";" <<
                  "dB(B)" << ";" <<
                  "C-Weighting"  << ";" <<
                  "dB(C)"  << ";" <<
                  endl;
        }

        int w=30;

        for (int j = 0; j < w; ++j) {

            double Sts[w];
            Sts[j]=Frequency[j]*D_starred[i]/bdata->m_a_tangential.value(i);

            double b[w];
            b[j]=qFabs(log10(Sts[j]/St2[i]));

            double B_min[w];
if (b[j]<0.13)
{B_min[j]=sqrt(16.888-886.788*pow(b[j],2));}
else if(b[j]>0.145)
{B_min[j]=-817.81*pow(b[j],3)+335.21*pow(b[j],2)-135.024*b[j]+10.619;}
else {B_min[j]=-83.607*b[j]+8.138;}

double B_max[w];
if (b[j]<0.1)
{B_max[j]=sqrt(16.888-886.788*pow(b[j],2))-4,109;}
else if(b[j]>0.187)
{B_max[j]=-80.541*pow(b[j],3)+44.174*pow(b[j],2)-39.381*b[j]+2.344;}
else {B_max[j]=-31.33*b[j]+1.854;}

double B_b[w];
B_b[j]=B_min[j]+BR_b0[i]*(B_max[j]-B_min[j]);

double SPL[w];
SPL[j]=first_term[i]+B_b[j]+K2[i];

double dB_A[w];
dB_A[j]=SPL[j]+AWeighting[j];

double dB_B[w];
dB_B[j]=SPL[j]+BWeighting[j];

double dB_C[w];
dB_C[j]=SPL[j]+CWeighting[j];

        stream << qSetFieldWidth(14)  <<
                  Frequency[j]  << ";" <<
                  Sts[j]  << ";" <<
                  b[j]  << ";" <<
                  B_min[j] << ";" <<
                  B_max[j] << ";" <<
                  B_b[j]  << ";" <<
                  SPL[j] << ";" <<
                  AWeighting[j]  << ";" <<
                  dB_A[j] << ";" <<
                  BWeighting[j] << ";" <<
                  dB_B[j] << ";" <<
                  CWeighting[j]  << ";" <<
                  dB_C[j]  << ";" <<
                  endl;

if(q_aux==w){q_aux=0;}
        }
        stream << endl;
        }
else {
//todo
            // if alpha is bigger than the switching angle
double first_term_a[number_of_segments];

first_term_a[i]=10*log10(pow(bdata->m_Mach.value(i),5)*L[i]*Dl[i]*D_starred[i]/pow(observer_position,2));

//first_term_a[i]=10*log10();



    y_aux=y_aux+1;
            double a0[number_of_segments];

            if (bdata->m_Reynolds.value(i)<95200)
            {a0[i]= 0.57;}
            else if (bdata->m_Reynolds.value(i)>857000)
            {a0[i]= 1.13;}
            else {a0[i]=-9.57*pow(10,-13)*(pow((3*bdata->m_Reynolds.value(i)-857000),2))+1.13;}

            double A_min_a0[number_of_segments];
            if (a0[i]<0.204)
            {A_min_a0[i]=sqrt(67.552-886.788*pow(a0[i],2))-8.219;}
            else if (a0[i]>0.244)
            {A_min_a0[i]=-142.795*pow(a0[i],3)+103.656*pow(a0[i],2)-57.757*a0[i]+6.006;}
            else{A_min_a0[i]=-32.665*a0[i]+3.981;}

            double A_max_a0[number_of_segments];
            if (a0[i]<0.13)
            {A_max_a0[i]=sqrt(67.552-886.788*pow(a0[i],2))-8.219;}
            else if (a0[i]>0.321)
            {A_max_a0[i]=-4.669*pow(a0[i],3)+3,491*pow(a0[i],2)-16.699*a0[i]+1.149;}
            else {A_max_a0[i]=-15.901*a0[i]+1.098;}

            double AR_a0[number_of_segments];
            AR_a0[i]=(-20-A_min_a0[i])/(A_max_a0[i]-A_min_a0[i]);

    if (y_aux==1){
            QString str= QString::number(z, 'f', 1);

                stream << "SPL Alpha: " <<  endl;
                stream << "Angles bigger than the switching angle: "  << endl;
                stream << "Tip Speed Ratio: " << str << endl;
                stream << "Section: " << (i+1)<<"/"<<number_of_segments << endl;
                stream << "D*: " << D_starred[i] << endl;
                stream << endl;

            stream << qSetFieldWidth(14)  <<
                      "Freq[Hz]"  << ";" <<
                      "Sts"  << ";" <<
                      "a"  << ";" <<
                      "A_min(b)" << ";" <<
                      "A_max(b)" << ";" <<
                      "A'(b)"  << ";" <<
                      "SPL [dB]" << ";" <<
                      "A-Weighting"  << ";" <<
                      "dB(A)" << ";" <<
                      "B-Weighting" << ";" <<
                      "dB(B)" << ";" <<
                      "C-Weighting"  << ";" <<
                      "dB(C)"  << ";" <<
                      endl;
            }

            int w=30;

            for (int j = 0; j < w; ++j) {

                double Sts[w];
                Sts[j]=Frequency[j]*D_starred[i]/bdata->m_a_tangential.value(i);

                double a[w];
                a[j]=qFabs(log10(Sts[j]/St2[i]));

                double A_min[w];

    if (a[j]<0.204)
    {A_min[j]=sqrt(67.552-886.788*pow(a[j],2))-8.219;}
    else if(a[j]>0.244)
    {A_min[j]=-142.795*pow(a[j],3)+103.656*pow(a[j],2)-57.757*a[j]+6.006;}
    else {A_min[j]=-32.665*a[j]+3.981;}

    double A_max[w];
    if (a[j]<0.13)
    {A_max[j]=sqrt(67.552-886.788*pow(a[j],2))-8.219;}
    else if(a[j]>0.321)
    {A_max[j]=-4.669*pow(a[j],3)+3.491*pow(a[j],2)-16.699*a[j]+1.149;}
    else {A_max[j]=-15.901*a[j]+1.098;}

    double Alin_a[w];
    Alin_a[j]=A_min[j]+AR_a0[i]*(A_max[j]-A_min[j]);

    double SPL[w];
    SPL[j]=first_term_a[i]+Alin_a[j]+K2[i];

    double dB_A[w];
    dB_A[j]=SPL[j]+AWeighting[j];

    double dB_B[w];
    dB_B[j]=SPL[j]+BWeighting[j];

    double dB_C[w];
    dB_C[j]=SPL[j]+CWeighting[j];

            stream << qSetFieldWidth(14)  <<
                      Frequency[j]  << ";" <<
                      Sts[j]  << ";" <<
                      a[j]  << ";" <<
                      A_min[j] << ";" <<
                      A_max[j] << ";" <<
                      Alin_a[j]  << ";" <<
                      SPL[j] << ";" <<
                      AWeighting[j]  << ";" <<
                      dB_A[j] << ";" <<
                      BWeighting[j] << ";" <<
                      dB_B[j] << ";" <<
                      CWeighting[j]  << ";" <<
                      dB_C[j]  << ";" <<
                      endl;

    if(y_aux==w){y_aux=0;}
            }
            stream << endl;
//todo
        }

        //SPL_Alpha ********************************************************


        //teste

        }
                z=z+ldelta;

        //planilha oculta

        //dados de entrada
        int sectionu = 21;
        double cl = 1.26;
        double cd = 0.0126;
        double ac = 0.2;
        int lbrac = 6;

        //definiçoes
        int sectionx = (sectionu-1);
        double wr=bdata->m_pos.value(sectionx)*omega;
        double chordx = bdata->m_c_local.value(sectionx);
        double radiusx = bdata->m_pos.value(sectionx);
        double doubt_calc=(chordx*3/(2*M_PI*radiusx));

//definição de arrays
double a_calc[number_of_segments+1];
double a_corrected_calc[number_of_segments+1];
double a_lin_calc[number_of_segments+1];
double a_variation_calc[number_of_segments+1];
double a_lin_variation_calc[number_of_segments+1];
double phi_calc[number_of_segments+1];
double F_calc[number_of_segments+1];
double K_calc[number_of_segments+1];
double ct_calc[number_of_segments+1];
double one_a_calc[number_of_segments+1];
double one_a_lin_calc[number_of_segments+1];
double theta_calc[number_of_segments+1];
double cl_calc[number_of_segments+1];
double cn_calc[number_of_segments+1];
double f_calc[number_of_segments+1];
double Fma_calc[number_of_segments+1];

       for (int i = -1; i < number_of_segments; ++i) {

//equações
//if (i==-1){

//if (z>=pSimuWidget->m_pctrlLELineEdit->getValue()+ldelta){
//    stream << endl;
//    stream << endl;
//    stream << "With Prandtl's Tip Loss Correction and Glauert Correction" << endl;
//    stream << "BEM Method. Hansen 2008 p. 50" << endl;
//    stream << endl;
//    stream << "Sect Defined:"   <<  (sectionx+1);
//    stream << endl;
//    stream << endl;
//    stream <<   qSetFieldWidth(14)  <<
//    "Sect"  << ";" <<
//    "a"  << ";" <<
//    "a (a>0.2)"  << ";" <<
//    "a'"  << ";" <<
//    "a var [%]"  << ";" <<
//    "a' var [%]"  << ";" <<
//    "(1-a)"  << ";" <<
//    "(1-a')"  << ";" <<
//    "V0"  << ";" <<
//    "Wr"  << ";" <<
//    "phi [deg]"  << ";" <<
//    "["  << ";" <<
//    "theta [deg]"  << ";" <<
//    "cl"  << ";" <<
//    "cd"  << ";" <<
//    "cn"  << ";" <<
//    "ct"  << ";" <<
//    "?"  << ";" <<
//    "f"  << ";" <<
//    "F"  << ";" <<
//    "K"  << ";" <<
//    endl;
//}
//    a_calc[-1]=0;
//    a_corrected_calc[-1]=0;
//    a_lin_calc[-1]=0;
//    a_variation_calc[-1]=0;
//    a_lin_variation_calc[-1]=0;
//    phi_calc[-1]=(qRadiansToDegrees(qAtan((one_a_calc[-1]*(double)approaxing_wind_speed)/(one_a_lin_calc[-1]*(double)wr))));
//    F_calc[-1]=0;
//    K_calc[-1]=0;
//    one_a_calc[-1]=1;
//    one_a_lin_calc[-1]=1;
//    theta_calc[-1]=(phi_calc[-1]-lbrac);
//    cn_calc[-1]=(cl*qCos(qDegreesToRadians(phi_calc[-1]))+cd*(double)qSin(qDegreesToRadians(phi_calc[-1])));
//    ct_calc[-1]=(cl*qSin(qDegreesToRadians(phi_calc[-1]))-cd*(double)qCos(qDegreesToRadians(phi_calc[-1])));
//    f_calc[-1]=(3*(double)(50-radiusx)/(double)(2*radiusx*(double)qSin(qDegreesToRadians(phi_calc[-1]))));
//    Fma_calc[-1]=(2/M_PI*(1/(double)qCos(qExp(-f_calc[-1]))));
//    K_calc[-1]=(4*(double)Fma_calc[-1]*pow(qSin(qDegreesToRadians(phi_calc[-1])),2)/doubt_calc*cn_calc[-1]);
//}
//else{
//        a_lin_calc[i] = (1.f/((4.f*Fma_calc[i-1]*qSin(qDegreesToRadians(phi_calc[i-1]))*qCos(qDegreesToRadians(phi_calc[i-1]))/doubt_calc*ct_calc[i-1])-1.f));

//        a_calc[i] = (1/(double)(K_calc[i-1]+1));
//        a_corrected_calc[i] = (0.5*(2+K_calc[i-1]*(double)(1-2*ac)-qSqrt(qAbs(pow((K_calc[i-1]*(1-2*ac)+2),2)+4*(double)(K_calc[i-1]*pow(ac,2)-1)))));

//        if (i==0){a_variation_calc[i] = ((a_calc[i]-a_calc[i-1])*(double)100);} else {
//        a_variation_calc[i] = (((a_calc[i]-a_calc[i-1])/a_calc[i-1])*(double)100);}

//        if (i==0){a_lin_variation_calc[i] = ((a_lin_calc[i]-a_lin_calc[i-1])*100);} else {
//        a_lin_variation_calc[i] = (((a_lin_calc[i]-a_lin_calc[i-1])/(double)a_lin_calc[i-1])*100);}

//        if (a_calc[i]<ac){one_a_calc[i]=(1-a_calc[i]);} else {one_a_calc[i]=(1-a_corrected_calc[i]);}

//        one_a_lin_calc[i]=(1+a_lin_calc[i]);

//        phi_calc[i]=(qRadiansToDegrees(qAtan((one_a_calc[i]*approaxing_wind_speed)/(double)(one_a_lin_calc[i]*wr))));

//        theta_calc[i]=(phi_calc[i]-lbrac);

//        cn_calc[i] = (cl*qCos(qDegreesToRadians(phi_calc[i]))+cd*qSin(qDegreesToRadians(phi_calc[i])));

//        ct_calc[i] = (cl*qSin(qDegreesToRadians(phi_calc[i]))-cd*qCos(qDegreesToRadians(phi_calc[i])));

//        f_calc[i]=(3*(50-radiusx)/(double)(2*radiusx*qSin(qDegreesToRadians(phi_calc[i]))));

//        Fma_calc[i]=(2/M_PI*(1/(double)qCos(qExp(-f_calc[i]))));

//        K_calc[i]=(4*Fma_calc[i]*pow(qSin(qDegreesToRadians(phi_calc[i])),2)/(double)doubt_calc*cn_calc[i]);
//}

//if (z>=pSimuWidget->m_pctrlLELineEdit->getValue()+ldelta){
//stream << qSetFieldWidth(14)  <<
//                  (i+1) << ";" <<
//                  a_calc[i] << ";" <<
//                  a_corrected_calc[i] << ";" <<
//                  a_lin_calc[i] << ";" <<
//                  a_variation_calc[i] << ";" <<
//                  a_lin_variation_calc[i] << ";" <<
//                  one_a_calc[i] << ";" <<
//                  one_a_lin_calc[i] << ";" <<
//                  approaxing_wind_speed << ";" <<
//                  wr << ";" <<
//                  phi_calc[i] << ";" <<
//                  lbrac << ";" <<
//                  theta_calc[i] << ";" <<
//                  cl << ";" <<
//                  cd << ";" <<
//                  cn_calc[i] << ";" <<
//                  ct_calc[i] << ";" <<
//                  doubt_calc << ";" <<
//                  f_calc[i] << ";" <<
//                  Fma_calc[i] << ";" <<
//                  K_calc[i] << ";" <<
//                endl;
//}
}
}
//       delete[] variables;
       qDeleteAll(noiseOpPoints);
}

//extra to make the 3d graphics
//        for (int j = 1; j < (m_parameter.sects+1); ++j) {
//            if (j > m_parameter.sects*0.85){
//               x = log10(m_parameter.sects*0.85-(j-m_parameter.sects*0.85))/log10(m_parameter.sects*0.85);
//        }
//                     else {x = log10(j)/log10(m_parameter.sects*0.85);}
//                        stream << j <<
//                        x <<
//                    m_parameter.chordBasedReynolds*(log10(j)/log10(m_parameter.sects))/0.3048  << // l*Re(for l=0.3048)/0.3048
//                      0.21/0.3048*log10(j)/log10(m_parameter.sects) << endl; // 0.21 (Mach for l=0.3048)/0.3048*l
//        }

//for (int j = int(noiseOpPoints.size()/2-2); j < int(noiseOpPoints.size()/2+2); ++j) {

//Sara

void NoiseSimulation::setAnalyzedOpPoints(QVector<OpPoint *> newList) {
    removeAllParents();
    for (int i = 0; i < newList.size(); ++i) {
        addParent(newList[i]);
    }

    m_parameter.analyzedOpPoints = newList;
}

QVariant NoiseSimulation::accessParameter(Parameter::NoiseSimulation::Key key, QVariant value) {
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

        //Sara
    case P::sects:
        if(set) {m_parameter.sects = value.toDouble(); if(m_parameter.sects<13) value=13;}
        else {
QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
value = pBEM->dlg_elements; break;
        }

    case P::rot_speed:
        if(set) {m_parameter.sects = value.toDouble();}
        else {
QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
value =0; break;
        }

    case P::u_wind_speed:
        if(set) {m_parameter.sects = value.toDouble();}
        else {
QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;

SimuWidget *pSimuWidget = (SimuWidget *) pBEM->m_pSimuWidget;
value = pSimuWidget->m_pctrlWindspeed->getValue();
 break;
        }

    case P::TSR:
        if(set) {m_parameter.sects = value.toDouble();}
        else {
QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
value = 0; break;
        }

//rot = m_pTData->Lambda0*windspeed*60/2/PI/m_pTData->OuterRadius;

        //Sara
    }

    return (set ? QVariant() : value);
}
