#include "NoiseSimulation.h"
//Sara
#include "NoiseCalculation.h"
#include "NoiseParameter.h"
#include "../XUnsteadyBEM/WindFieldModule.h"
#include "../XLLT/QLLTSimulation.h"
#include "../XLLT/QLLTCreatorDialog.h"
//Sara

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
//Sara begin
        if (m_parameter.qs3DSim==0){//rotor
        const int index = getAvailableVariables_rotor().indexOf(i == 0 ? xAxis : yAxis);
        QVector<double> *vector = (i == 0 ? &xVector : &yVector);
        switch (index) {
            case 0: *vector = NoiseCalculation::CENTRAL_BAND_FREQUENCY; break;
            case 1: *vector = m_calculation.SPLadB()[opPointIndex]; zeroY = true; break;
            case 2: *vector = m_calculation.SPLsdB()[opPointIndex]; zeroY = true; break;
            case 3: *vector = m_calculation.SPLpdB()[opPointIndex]; zeroY = true; break;
            case 4: *vector = m_calculation.SPL_LEdB()[opPointIndex]; zeroY = true; break; //Alexandre MOD Sara
            case 5: *vector = m_calculation.SPLdB()[opPointIndex]; break;
            case 6: *vector = m_calculation.SPLdBAW()[opPointIndex]; break;
            case 7: *vector = m_calculation.SPLdBBW()[opPointIndex]; break;
            case 8: *vector = m_calculation.SPLdBCW()[opPointIndex]; break;
            case 9: *vector = m_calculation.SPLadB3d()[opPointIndex]; zeroY = true; break;
            case 10: *vector = m_calculation.SPLsdB3d()[opPointIndex]; zeroY = true; break;
            case 11: *vector = m_calculation.SPLpdB3d()[opPointIndex]; zeroY = true; break;
            case 12: *vector = m_calculation.SPL_LEdB3d()[opPointIndex]; zeroY = true; break;
            case 13: *vector = m_calculation.SPLdB3d()[opPointIndex]; zeroY = true; break;
            case 14: *vector = m_calculation.SPLdBAW3d()[opPointIndex]; break;
            case 15: *vector = m_calculation.SPLdBBW3d()[opPointIndex]; break;
            case 16: *vector = m_calculation.SPLdBCW3d()[opPointIndex]; break;
            case 17: *vector = m_calculation.SPLadB3d_final()[opPointIndex]; zeroY = true; break;
            case 18: *vector = m_calculation.SPLsdB3d_final()[opPointIndex]; zeroY = true; break;
            case 19: *vector = m_calculation.SPLpdB3d_final()[opPointIndex]; zeroY = true; break;
            case 20: *vector = m_calculation.SPL_LEdB3d_final()[opPointIndex]; zeroY = true; break;
            case 21: *vector = m_calculation.SPLdB3d_final()[opPointIndex]; zeroY = true; break;
            case 22: *vector = m_calculation.SPLdBAW3d_final()[opPointIndex]; break;
            case 23: *vector = m_calculation.SPLdBBW3d_final()[opPointIndex]; break;
            case 24: *vector = m_calculation.SPLdBCW3d_final()[opPointIndex]; break;
            case 25: *vector = m_calculation.SPLadB3d_final_rotor()[opPointIndex]; zeroY = true; break;
            case 26: *vector = m_calculation.SPLsdB3d_final_rotor()[opPointIndex]; zeroY = true; break;
            case 27: *vector = m_calculation.SPLpdB3d_final_rotor()[opPointIndex]; zeroY = true; break;
            case 28: *vector = m_calculation.SPL_LEdB3d_final_rotor()[opPointIndex]; zeroY = true; break;
            case 29: *vector = m_calculation.SPLdB3d_final_rotor()[opPointIndex]; zeroY = true; break;
            case 30: *vector = m_calculation.SPLdBAW3d_final_rotor()[opPointIndex]; break;
            case 31: *vector = m_calculation.SPLdBBW3d_final_rotor()[opPointIndex]; break;
            case 32: *vector = m_calculation.SPLdBCW3d_final_rotor()[opPointIndex]; break;
            case 33: *vector = m_calculation.SPLadB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 34: *vector = m_calculation.SPLsdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 35: *vector = m_calculation.SPLpdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 36: *vector = m_calculation.SPL_LEdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 37: *vector = m_calculation.SPLdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 38: *vector = m_calculation.SPLdBAW3d_final_rotor_loops()[opPointIndex]; break;
            case 39: *vector = m_calculation.SPLdBBW3d_final_rotor_loops()[opPointIndex]; break;
            case 40: *vector = m_calculation.SPLdBCW3d_final_rotor_loops()[opPointIndex]; break;

            default: return nullptr;
            }
        }

        if (m_parameter.qs3DSim==1){//blade
            const int index = getAvailableVariables_blade().indexOf(i == 0 ? xAxis : yAxis);
            QVector<double> *vector = (i == 0 ? &xVector : &yVector);
            switch (index) {
            case 0: *vector = NoiseCalculation::CENTRAL_BAND_FREQUENCY; break;
            case 1: *vector = m_calculation.SPLadB()[opPointIndex]; zeroY = true; break;
            case 2: *vector = m_calculation.SPLsdB()[opPointIndex]; zeroY = true; break;
            case 3: *vector = m_calculation.SPLpdB()[opPointIndex]; zeroY = true; break;
            case 4: *vector = m_calculation.SPL_LEdB()[opPointIndex]; zeroY = true; break; //Alexandre MOD Sara
            case 5: *vector = m_calculation.SPLdB()[opPointIndex]; break;
            case 6: *vector = m_calculation.SPLdBAW()[opPointIndex]; break;
            case 7: *vector = m_calculation.SPLdBBW()[opPointIndex]; break;
            case 8: *vector = m_calculation.SPLdBCW()[opPointIndex]; break;
            case 9: *vector = m_calculation.SPLadB3d()[opPointIndex]; zeroY = true; break;
            case 10: *vector = m_calculation.SPLsdB3d()[opPointIndex]; zeroY = true; break;
            case 11: *vector = m_calculation.SPLpdB3d()[opPointIndex]; zeroY = true; break;
            case 12: *vector = m_calculation.SPL_LEdB3d()[opPointIndex]; zeroY = true; break;
            case 13: *vector = m_calculation.SPLdB3d()[opPointIndex]; zeroY = true; break;
            case 14: *vector = m_calculation.SPLdBAW3d()[opPointIndex]; break;
            case 15: *vector = m_calculation.SPLdBBW3d()[opPointIndex]; break;
            case 16: *vector = m_calculation.SPLdBCW3d()[opPointIndex]; break;
            case 17: *vector = m_calculation.SPLadB3d_final()[opPointIndex]; zeroY = true; break;
            case 18: *vector = m_calculation.SPLsdB3d_final()[opPointIndex]; zeroY = true; break;
            case 19: *vector = m_calculation.SPLpdB3d_final()[opPointIndex]; zeroY = true; break;
            case 20: *vector = m_calculation.SPL_LEdB3d_final()[opPointIndex]; zeroY = true; break;
            case 21: *vector = m_calculation.SPLdB3d_final()[opPointIndex]; zeroY = true; break;
            case 22: *vector = m_calculation.SPLdBAW3d_final()[opPointIndex]; break;
            case 23: *vector = m_calculation.SPLdBBW3d_final()[opPointIndex]; break;
            case 24: *vector = m_calculation.SPLdBCW3d_final()[opPointIndex]; break;

            default: return nullptr;
            }
        }

        if (m_parameter.qs3DSim==2){
            const int index = getAvailableVariables().indexOf(i == 0 ? xAxis : yAxis);
            QVector<double> *vector = (i == 0 ? &xVector : &yVector);
        switch (index) {
            case 0: *vector = NoiseCalculation::CENTRAL_BAND_FREQUENCY; break;
            case 1: *vector = m_calculation.SPLadB()[opPointIndex]; zeroY = true; break;
            case 2: *vector = m_calculation.SPLsdB()[opPointIndex]; zeroY = true; break;
            case 3: *vector = m_calculation.SPLpdB()[opPointIndex]; zeroY = true; break;
            case 4: *vector = m_calculation.SPL_LEdB()[opPointIndex]; zeroY = true; break; //Alexandre MOD Sara
            case 5: *vector = m_calculation.SPLdB()[opPointIndex]; break;
            case 6: *vector = m_calculation.SPLdBAW()[opPointIndex]; break;
            case 7: *vector = m_calculation.SPLdBBW()[opPointIndex]; break;
            case 8: *vector = m_calculation.SPLdBCW()[opPointIndex]; break;
            case 9: *vector = m_calculation.SPLadB3d()[opPointIndex]; zeroY = true; break;
            case 10: *vector = m_calculation.SPLsdB3d()[opPointIndex]; zeroY = true; break;
            case 11: *vector = m_calculation.SPLpdB3d()[opPointIndex]; zeroY = true; break;
            case 12: *vector = m_calculation.SPL_LEdB3d()[opPointIndex]; zeroY = true; break;
            case 13: *vector = m_calculation.SPLdB3d()[opPointIndex]; zeroY = true; break;
            case 14: *vector = m_calculation.SPLdBAW3d()[opPointIndex]; break;
            case 15: *vector = m_calculation.SPLdBBW3d()[opPointIndex]; break;
            case 16: *vector = m_calculation.SPLdBCW3d()[opPointIndex]; break;

            default: return nullptr;
            }
        }
//Sara end
}

    NewCurve *curve = new NewCurve (this);
//    curve->setAllPoints(xVector.data(), yVector.data(), xVector.size());
    for (int i = 0; i < xVector.size(); ++i) {  // zero the y values lower 0 for certain outputs
        curve->addPoint(xVector[i], (zeroY && yVector[i] < 0 ? 0.0 : yVector[i]));
    }
    return curve;
}

//Sara begin
QStringList NoiseSimulation::getAvailableVariables_rotor(NewGraph::GraphType /*graphType*/) {
    QStringList variables;
    // WARNING: when changing any variables list, change newCurve as well!
    variables << "Freq [Hz]" << "SPL_alpha" << "SPL_S" << "SPL_P" << "SPL_LE (dB)" << "SPL (dB)" << "SPL (dB(A))" << "SPL (dB(B))" << "SPL (dB(C))" <<  "SPL_alpha[multi]" << "SPL_S[multi]" << "SPL_P[multi]" << "SPL_LE[multi] (dB)" << "SPL[multi] (dB)" << "SPL[multi] (dB(A))" << "SPL[multi] (dB(B))" << "SPL[multi] (dB(C))" << "SPL_alpha_blade[qs3D]" << "SPL_S_blade[qs3D]" << "SPL_P_blade[qs3D]" << "SPL_LE_blade[qs3D] (dB)" << "SPL_blade[qs3D] (dB)" << "SPL_blade[qs3D] (dB(A))" << "SPL_blade[qs3D] (dB(B))" << "SPL_blade[qs3D] (dB(C))" << "SPL_alpha_rotor[qs3D]" << "SPL_S_rotor[qs3D]" << "SPL_P_rotor[qs3D]" << "SPL_LE_rotor[qs3D] (dB)" << "SPL_rotor[qs3D] (dB)" << "SPL_rotor[qs3D] (dB(A))" << "SPL_rotor[qs3D] (dB(B))" << "SPL_rotor[qs3D] (dB(C))" << "SPL_alpha_rotor_loops[qs3D]" << "SPL_S_rotor_loops[qs3D]" << "SPL_P_rotor_loops[qs3D]" << "SPL_LE_rotor_loops[qs3D] (dB)" << "SPL_rotor_loops[qs3D] (dB)" << "SPL_rotor_loops[qs3D] (dB(A))" << "SPL_rotor_loops[qs3D] (dB(B))" << "SPL_rotor_loops[qs3D] (dB(C))"; //Alexandre MOD Sara

    return variables;
}

QStringList NoiseSimulation::getAvailableVariables_blade(NewGraph::GraphType /*graphType*/) {
    QStringList variables;
    // WARNING: when changing any variables list, change newCurve as well!
    variables << "Freq [Hz]" << "SPL_alpha" << "SPL_S" << "SPL_P" << "SPL_LE (dB)" << "SPL (dB)" << "SPL (dB(A))" << "SPL (dB(B))" << "SPL (dB(C))" <<  "SPL_alpha[multi]" << "SPL_S[multi]" << "SPL_P[multi]" << "SPL_LE[multi] (dB)" << "SPL[multi] (dB)" << "SPL[multi] (dB(A))" << "SPL[multi] (dB(B))" << "SPL[multi] (dB(C))" << "SPL_alpha_blade[qs3D]" << "SPL_S_blade[qs3D]" << "SPL_P_blade[qs3D]" << "SPL_LE_blade[qs3D] (dB)" << "SPL_blade[qs3D] (dB)" << "SPL_blade[qs3D] (dB(A))" << "SPL_blade[qs3D] (dB(B))" << "SPL_blade[qs3D] (dB(C))"; //Alexandre MOD Sara

    return variables;
}
//Sara end

QStringList NoiseSimulation::getAvailableVariables(NewGraph::GraphType /*graphType*/) {
    QStringList variables;
    // WARNING: when changing any variables list, change newCurve as well!
    variables << "Freq [Hz]" << "SPL_alpha" << "SPL_S" << "SPL_P" << "SPL_LE (dB)" << "SPL (dB)" << "SPL (dB(A))" << "SPL (dB(B))" << "SPL (dB(C))"; //Alexandre MOD Sara

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
        //Sara begin
        if (g_360PolarStore.isEmpty() && g_qbem->m_pCur360Polar == NULL) {
            message.prepend(tr("- No 360 Polar in Database"));
            }
        if (g_rotorStore.isEmpty() && g_qbem->m_pBlade == NULL) {
            message.prepend(tr("- No HAWT Blade in Database"));
        }
        if (g_windFieldStore.size() == 0) {
            message.prepend(tr("- No Windfield in Database (Unsteady Case)"));
        }
        message.prepend("- No Noise Simulation in Database");
            //Sara end
            return message;
    } else {
        return QStringList();
    }
}

void NoiseSimulation::simulate() {
    NoiseCreatorDialog *pNoiseCreatorDialog = (NoiseCreatorDialog *) g_mainFrame->m_pBEM;//Sara
    pNoiseCreatorDialog->OnProgressDlg();//Sara
    m_calculation.setNoiseParam(&m_parameter);
    m_calculation.setInitialValues();//Sara

if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculate();}
    if (m_parameter.qs3DSim!=2){
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculateqs3d_graphics_loops();}//Sara
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculateqs3d_blade();}//Sara
    }
    if(m_parameter.qs3DSim==0){
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculateqs3d_rotor();}//Sara
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculateqs3d_rotor_loops();}//Sara
    }
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){pNoiseCreatorDialog->m_progress_dlg->cancel();}
}

void NoiseSimulation::exportCalculation(QTextStream &stream) {
    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Noise prediction file export" << endl;
    stream << endl;
    if(m_parameter.Lowson_type==1){stream << "leading edge model: Von Karman" <<endl;}
    if(m_parameter.Lowson_type==2){stream << "leading edge model: Rapid Distortion" <<endl;}//era m_calculation.m_
//    stream << "constante D: " << m_calculation.d_const << endl;
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();
    for (int i = 0; i < noiseOpPoints.size(); ++i) {
        stream << qSetFieldWidth(0);
        stream << "Alpha: " << noiseOpPoints[i]->getAlphaDegree() <<
                  ", Re = " << noiseOpPoints[i]->getReynolds() << endl;
        stream << "SPL_a: " << m_calculation.SPLALOG()[i] << "" << endl;
        stream << "SPL_s: " << m_calculation.SPLSLOG()[i] << "" << endl;
        stream << "SPL_p: " << m_calculation.SPLPLOG()[i] << "" << endl;
        if(m_parameter.Lowson_type!=0){
        stream << "SPL_LE: " << m_calculation.SPLlogLE()[i] << " dB" << endl;
        }
        stream << "OASPL: " << m_calculation.OASPL()[i] << " dB" << endl;
        stream << "OASPL (A): " << m_calculation.OASPLA()[i] << " dB(A)" << endl;
        stream << "OASPL (B): " << m_calculation.OASPLB()[i] << " dB(B)" << endl;
        stream << "OASPL (C): " << m_calculation.OASPLC()[i] << " dB(C)" << endl;
        stream << endl;
               if(m_parameter.Lowson_type!=0){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" <<
                  "SPLa" <<
                  "SPLs" <<
                  "SPLp" <<
                  "SPL_LE (dB)" <<
                  "SPL (dB)" <<
                  "SPL (dB(A))" <<
                  "SPL (dB(B))" <<
                  "SPL (dB(C))" <<
                  endl; //Alexandre MOD
               }
               else{
                   stream << qSetFieldWidth(14) <<
                             "Freq [Hz]" <<
                             "SPLa" <<
                             "SPLs" <<
                             "SPLp" <<
                             "SPL (dB)" <<
                             "SPL (dB(A))" <<
                             "SPL (dB(B))" <<
                             "SPL (dB(C))" <<endl; //Alexandre MOD
               }

        for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

 if(m_parameter.Lowson_type!=0){
            stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] <<
                      m_calculation.SPLadB()[i][j] <<
                      m_calculation.SPLsdB()[i][j] <<
                      m_calculation.SPLpdB()[i][j] <<
                      m_calculation.SPL_LEdB()[i][j] <<
                      m_calculation.SPLdB()[i][j] <<
                      m_calculation.SPLdBAW()[i][j] <<
                      m_calculation.SPLdBBW()[i][j] <<
                      m_calculation.SPLdBCW()[i][j] <<
                      endl; //Alexandre MOD
        }
else{
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] <<
               m_calculation.SPLadB()[i][j] <<
               m_calculation.SPLsdB()[i][j] <<
               m_calculation.SPLpdB()[i][j] <<
               m_calculation.SPLdB()[i][j] <<
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

void NoiseSimulation::exportCalculationqs3DNoise_blade(QTextStream &stream) {
    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Noise prediction file export quasi 3D for blade" << endl;
    stream << endl;
    if(m_parameter.Lowson_type==1){stream << "leading edge model: Von Karman" <<endl;}
    if(m_parameter.Lowson_type==2){stream << "leading edge model: Rapid Distortion" <<endl;}
    stream << "Tip Speed Ratio: " << m_parameter.TSRtd << endl;
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    int number_of_segments = pBEM->dlg_elements;
    stream << endl;
    stream << endl;
    stream << "*********************************************" << endl;
    stream << "Total Values:" << endl;
    stream << endl;
    stream << "OASPL: " << m_calculation.Final_qs3d << endl;
    stream << "SPL alfa: " << m_calculation.Final_qs3d_alpha << endl;
    stream << "SPL S: " << m_calculation.Final_qs3d_S << endl;
    stream << "SPL P: " << m_calculation.Final_qs3d_P << endl;
    if (m_parameter.Lowson_type!=0){stream << "SPL LE: " << m_calculation.Final_qs3d_LE;}
    stream << endl;
    stream << "*********************************************" << endl;
    stream << endl;
    stream << endl;

    for (int i = 0; i < number_of_segments; ++i){
stream << qSetFieldWidth(0);
//        stream << endl;
stream << "Section: " << (i+1)<<"/"<<number_of_segments << endl;
        stream << endl;
               if(m_parameter.Lowson_type!=0){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" << ";" <<
                  "SPLa" << ";" <<
                  "SPLs" << ";" <<
                  "SPLp" << ";" <<
                  "SPL_LE (dB)" << ";" <<
                  "SPL (dB)" << ";" <<
                  "SPL (dB(A))" << ";" <<
                  "SPL (dB(B))" << ";" <<
                  "SPL (dB(C))" << ";" <<endl; //Alexandre MOD
               }
               else{
                   stream << qSetFieldWidth(14) <<
                             "Freq [Hz]" << ";" <<
                             "SPLa" << ";" <<
                             "SPLs" << ";" <<
                             "SPLp" << ";" <<
                             "SPL (dB)" << ";" <<
                             "SPL (dB(A))" << ";" <<
                             "SPL (dB(B))" << ";" <<
                             "SPL (dB(C))" <<endl; //Alexandre MOD
               }
//int i=0;
        for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

            if(m_parameter.Lowson_type!=0){
                       stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                                 m_calculation.SPLadB3d()[i][j] << ";" <<
                                 m_calculation.SPLsdB3d()[i][j] << ";" <<
                                 m_calculation.SPLpdB3d()[i][j] << ";" <<
                                 m_calculation.SPL_LEdB3d()[i][j] << ";" <<
                                 m_calculation.SPLdB3d()[i][j] << ";" <<
                                 m_calculation.SPLdBAW3d()[i][j] << ";" <<
                                 m_calculation.SPLdBBW3d()[i][j] << ";" <<
                                 m_calculation.SPLdBCW3d()[i][j] << ";" <<
                                 endl; //Alexandre MOD
                   }
           else{
                stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                          m_calculation.SPLadB3d()[i][j] << ";" <<
                          m_calculation.SPLsdB3d()[i][j] << ";" <<
                          m_calculation.SPLpdB3d()[i][j] << ";" <<
                          m_calculation.SPLdB3d()[i][j] << ";" <<
                          m_calculation.SPLdBAW3d()[i][j] << ";" <<
                          m_calculation.SPLdBBW3d()[i][j] << ";" <<
                          m_calculation.SPLdBCW3d()[i][j] << ";" << endl;
            }
        }
        stream << endl;
        stream << qSetFieldWidth(0);
        stream << "SPL_a: " << m_calculation.SPLALOG3d()[i] << "" << endl;
        stream << "SPL_s: " << m_calculation.SPLSLOG3d()[i] << "" << endl;
        stream << "SPL_p: " << m_calculation.SPLPLOG3d()[i] << "" << endl;
        if(m_parameter.Lowson_type!=0){
        stream << "SPL_LE: " << m_calculation.SPLlogLE3d()[i] << " dB" << endl;}
        stream << "OASPL: " << m_calculation.OASPL3d()[i] << " dB" << endl;
        stream << "OASPL (A): " << m_calculation.OASPLA3d()[i] << " dB(A)" << endl;
        stream << "OASPL (B): " << m_calculation.OASPLB3d()[i] << " dB(B)" << endl;
        stream << "OASPL (C): " << m_calculation.OASPLC3d()[i] << " dB(C)" << endl;
        stream << endl;
        stream << endl;

if(i==(number_of_segments-1)){
stream << "********** FINAL **********" << endl;
if(m_parameter.Lowson_type!=0){
stream << qSetFieldWidth(14) <<
   "Freq [Hz]" << ";" <<
   "SPLa" << ";" <<
   "SPLs" << ";" <<
   "SPLp" << ";" <<
   "SPL_LE (dB)" << ";" <<
   "SPL (dB)" << ";" <<
   "SPL (dB(A))" << ";" <<
   "SPL (dB(B))" << ";" <<
   "SPL (dB(C))" << ";" <<endl; //Alexandre MOD
}
else{
    stream << qSetFieldWidth(14) <<
              "Freq [Hz]" << ";" <<
              "SPLa" << ";" <<
              "SPLs" << ";" <<
              "SPLp" << ";" <<
              "SPL (dB)" << ";" <<
              "SPL (dB(A))" << ";" <<
              "SPL (dB(B))" << ";" <<
              "SPL (dB(C))" <<endl; //Alexandre MOD
}
for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

if(m_parameter.Lowson_type!=0){
        stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                  m_calculation.SPLadB3d_final()[i][j] << ";" <<
                  m_calculation.SPLsdB3d_final()[i][j] << ";" <<
                  m_calculation.SPLpdB3d_final()[i][j] << ";" <<
                  m_calculation.SPL_LEdB3d_final()[i][j] << ";" <<
                  m_calculation.SPLdB3d_final()[i][j] << ";" <<
                  m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
                  m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
                  m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
                  endl; //Alexandre MOD
    }
else{
 stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
           m_calculation.SPLadB3d_final()[i][j] << ";" <<
           m_calculation.SPLsdB3d_final()[i][j] << ";" <<
           m_calculation.SPLpdB3d_final()[i][j] << ";" <<
           m_calculation.SPLdB3d_final()[i][j] << ";" <<
           m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
           m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
           m_calculation.SPLdBCW3d_final()[i][j] << ";" << endl;
}
}
stream << endl;
}}}

void NoiseSimulation::exportCalculationqs3DNoise_rotor(QTextStream &stream) {
    double E=m_parameter.initial_azimuth;

    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Noise prediction file export quasi 3D for rotor" << endl;
    stream << endl;
    if(m_parameter.Lowson_type==1){stream << "leading edge model: Von Karman" <<endl;}
    if(m_parameter.Lowson_type==2){stream << "leading edge model: Rapid Distortion" <<endl;}
    stream << "Tip Speed Ratio: " << m_parameter.TSRtd << endl;
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    int number_of_segments = pBEM->dlg_elements;
    stream << endl;
    stream << endl;
    stream << "*********************************************" << endl;
    stream << "Total Values:" << endl;
    stream << endl;
    stream << "OASPL: " << m_calculation.Final_qs3d_rotor << endl;
    stream << "SPL alfa: " << m_calculation.Final_qs3d_alpha_rotor << endl;
    stream << "SPL S: " << m_calculation.Final_qs3d_S_rotor << endl;
    stream << "SPL P: " << m_calculation.Final_qs3d_P_rotor << endl;
    if (m_parameter.Lowson_type!=0){stream << "SPL LE: " << m_calculation.Final_qs3d_LE_rotor;}
    stream << endl;
    stream << "*********************************************" << endl;
    stream << endl;
    stream << endl;

    for (int i = 0; i < number_of_segments; ++i){
stream << qSetFieldWidth(0);
//        stream << endl;
stream << "Section: " << (i+1)<<"/"<<number_of_segments << endl;
        stream << endl;
               if(m_parameter.Lowson_type!=0){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" << ";" <<
                  "SPLa" << ";" <<
                  "SPLs" << ";" <<
                  "SPLp" << ";" <<
                  "SPL_LE (dB)" << ";" <<
                  "SPL (dB)" << ";" <<
                  "SPL (dB(A))" << ";" <<
                  "SPL (dB(B))" << ";" <<
                  "SPL (dB(C))" << ";" <<endl; //Alexandre MOD
               }
               else{
                   stream << qSetFieldWidth(14) <<
                             "Freq [Hz]" << ";" <<
                             "SPLa" << ";" <<
                             "SPLs" << ";" <<
                             "SPLp" << ";" <<
                             "SPL (dB)" << ";" <<
                             "SPL (dB(A))" << ";" <<
                             "SPL (dB(B))" << ";" <<
                             "SPL (dB(C))" <<endl; //Alexandre MOD
               }
//int i=0;
        for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

            if(m_parameter.Lowson_type!=0){
                       stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                                 m_calculation.SPLadB3d_rotor()[i][j] << ";" <<
                                 m_calculation.SPLsdB3d_rotor()[i][j] << ";" <<
                                 m_calculation.SPLpdB3d_rotor()[i][j] << ";" <<
                                 m_calculation.SPL_LEdB3d_rotor()[i][j] << ";" <<
                                 m_calculation.SPLdB3d_rotor()[i][j] << ";" <<
                                 m_calculation.SPLdBAW3d_rotor()[i][j] << ";" <<
                                 m_calculation.SPLdBBW3d_rotor()[i][j] << ";" <<
                                 m_calculation.SPLdBCW3d_rotor()[i][j] << ";" <<
                                 endl; //Alexandre MOD
                   }
           else{
                stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                          m_calculation.SPLadB3d_rotor()[i][j] << ";" <<
                          m_calculation.SPLsdB3d_rotor()[i][j] << ";" <<
                          m_calculation.SPLpdB3d_rotor()[i][j] << ";" <<
                          m_calculation.SPLdB3d_rotor()[i][j] << ";" <<
                          m_calculation.SPLdBAW3d_rotor()[i][j] << ";" <<
                          m_calculation.SPLdBBW3d_rotor()[i][j] << ";" <<
                          m_calculation.SPLdBCW3d_rotor()[i][j] << ";" << endl;
            }
        }
        stream << endl;
        stream << qSetFieldWidth(0);
        stream << "SPL_a: " << m_calculation.SPLALOG3d_rotor()[i] << "" << endl;
        stream << "SPL_s: " << m_calculation.SPLSLOG3d_rotor()[i] << "" << endl;
        stream << "SPL_p: " << m_calculation.SPLPLOG3d_rotor()[i] << "" << endl;
        if(m_parameter.Lowson_type!=0){
        stream << "SPL_LE: " << m_calculation.SPLlogLE3d_rotor()[i] << " dB" << endl;}
        stream << "OASPL: " << m_calculation.OASPL3d_rotor()[i] << " dB" << endl;
        stream << "OASPL (A): " << m_calculation.OASPLA3d_rotor()[i] << " dB(A)" << endl;
        stream << "OASPL (B): " << m_calculation.OASPLB3d_rotor()[i] << " dB(B)" << endl;
        stream << "OASPL (C): " << m_calculation.OASPLC3d_rotor()[i] << " dB(C)" << endl;
        stream << endl;
        stream << endl;

if(i==(number_of_segments-1)){
stream << "********** FINAL **********" << endl;
if(m_parameter.Lowson_type!=0){
stream << qSetFieldWidth(14) <<
   "Freq [Hz]" << ";" <<
   "SPLa" << ";" <<
   "SPLs" << ";" <<
   "SPLp" << ";" <<
   "SPL_LE (dB)" << ";" <<
   "SPL (dB)" << ";" <<
   "SPL (dB(A))" << ";" <<
   "SPL (dB(B))" << ";" <<
   "SPL (dB(C))" << ";" <<endl; //Alexandre MOD
}
else{
    stream << qSetFieldWidth(14) <<
              "Freq [Hz]" << ";" <<
              "SPLa" << ";" <<
              "SPLs" << ";" <<
              "SPLp" << ";" <<
              "SPL (dB)" << ";" <<
              "SPL (dB(A))" << ";" <<
              "SPL (dB(B))" << ";" <<
              "SPL (dB(C))" <<endl; //Alexandre MOD
}
for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

if(m_parameter.Lowson_type!=0){
        stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                  m_calculation.SPLadB3d_final_rotor()[i][j] << ";" <<
                  m_calculation.SPLsdB3d_final_rotor()[i][j] << ";" <<
                  m_calculation.SPLpdB3d_final_rotor()[i][j] << ";" <<
                  m_calculation.SPL_LEdB3d_final_rotor()[i][j] << ";" <<
                  m_calculation.SPLdB3d_final_rotor()[i][j] << ";" <<
                  m_calculation.SPLdBAW3d_final_rotor()[i][j] << ";" <<
                  m_calculation.SPLdBBW3d_final_rotor()[i][j] << ";" <<
                  m_calculation.SPLdBCW3d_final_rotor()[i][j] << ";" <<
                  endl; //Alexandre MOD
    }
else{
 stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
           m_calculation.SPLadB3d_final_rotor()[i][j] << ";" <<
           m_calculation.SPLsdB3d_final_rotor()[i][j] << ";" <<
           m_calculation.SPLpdB3d_final_rotor()[i][j] << ";" <<
           m_calculation.SPLdB3d_final_rotor()[i][j] << ";" <<
           m_calculation.SPLdBAW3d_final_rotor()[i][j] << ";" <<
           m_calculation.SPLdBBW3d_final_rotor()[i][j] << ";" <<
           m_calculation.SPLdBCW3d_final_rotor()[i][j] << ";" << endl;
}
}
stream << endl;
}}}

void NoiseSimulation::exportCalculationqs3DNoise_rotor_loops(QTextStream &stream) {
    double E=m_parameter.initial_azimuth;

    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Noise prediction file export quasi 3D for rotor in rotation movement" << endl;
    stream << endl;
    if(m_parameter.Lowson_type==1){stream << "leading edge model: Von Karman" <<endl;}
    if(m_parameter.Lowson_type==2){stream << "leading edge model: Rapid Distortion" <<endl;}
    stream << "Tip Speed Ratio: " << m_parameter.TSRtd << endl;
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
     int number_of_segments = pBEM->dlg_elements;
     stream << endl;
     stream << endl;
     stream << "*********************************************" << endl;
     stream << "Total Values:" << endl;
     stream << endl;
     stream << "OASPL: " << m_calculation.Final_qs3d_rotor_loops << endl;
     stream << "SPL alfa: " << m_calculation.Final_qs3d_alpha_rotor_loops << endl;
     stream << "SPL S: " << m_calculation.Final_qs3d_S_rotor_loops << endl;
     stream << "SPL P: " << m_calculation.Final_qs3d_P_rotor_loops << endl;
     if (m_parameter.Lowson_type!=0){stream << "SPL LE: " << m_calculation.Final_qs3d_LE_rotor_loops;}
     stream << endl;
     stream << "*********************************************" << endl;
     stream << endl;
     stream << endl;

     for (int i = 0; i < number_of_segments; ++i){
 stream << qSetFieldWidth(0);

 if(i==(number_of_segments-1)){
// stream << "********** FINAL **********" << endl;
 if(m_parameter.Lowson_type!=0){
 stream << qSetFieldWidth(14) <<
    "Freq [Hz]" << ";" <<
    "SPLa" << ";" <<
    "SPLs" << ";" <<
    "SPLp" << ";" <<
    "SPL_LE (dB)" << ";" <<
    "SPL (dB)" << ";" <<
    "SPL (dB(A))" << ";" <<
    "SPL (dB(B))" << ";" <<
    "SPL (dB(C))" << ";" <<endl; //Alexandre MOD
 }
 else{
     stream << qSetFieldWidth(14) <<
               "Freq [Hz]" << ";" <<
               "SPLa" << ";" <<
               "SPLs" << ";" <<
               "SPLp" << ";" <<
               "SPL (dB)" << ";" <<
               "SPL (dB(A))" << ";" <<
               "SPL (dB(B))" << ";" <<
               "SPL (dB(C))" <<endl; //Alexandre MOD
 }
 for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

 if(m_parameter.Lowson_type!=0){
         stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                   m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
                   endl; //Alexandre MOD
     }
 else{
  stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
            m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
            m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
            m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
            m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
            m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
            m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
            m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" << endl;
}
}
stream << endl;
}}}

void NoiseSimulation::exportqs3DLog(QTextStream &stream) {
    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Quasi 3D Noise Log" << endl;
    stream << endl;

QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
double z=lstart;
double approaxing_wind_speed = m_parameter.u_wind_speed;

QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
double blade_pitch=pbem->m_pctrlFixedPitch->getValue();

    foreach(BData * bdata, pbem->m_pBEMData->GetBData()){
            if (z==m_parameter.TSRtd){

    int number_of_segments = bdata->m_pos.size();
    double lambda = pbem->dlg_lambda;
    int mpos_size = bdata->m_pos.size(); //total number of segments
    double finalradius = bdata->m_pos.value(mpos_size-1);
    double nom_tg_speed = approaxing_wind_speed*lambda;
    double omega = nom_tg_speed/finalradius;

    //definitions
    double axial_ind_fact[number_of_segments];
    double axial_ind_fact_n[number_of_segments];
    double axial_velocity[number_of_segments];
    double tangential_speed[number_of_segments];
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
    double r_R[number_of_segments];
    double c_Rx[number_of_segments];
    double D_starred_C_HT[number_of_segments];
    double D_starred_HT[number_of_segments];
    double D_starred_C_N[number_of_segments];
    double D_starred_N[number_of_segments];
    double Mach[number_of_segments];
    double corr_fact[number_of_segments];
    double D_starred_HT_S[number_of_segments];
    double D_starred_HT_P[number_of_segments];
    double D_starred_N_S[number_of_segments];
    double D_starred_N_P[number_of_segments];
    double SwAlpha[number_of_segments];
    double SwAlpha_1[number_of_segments];
    double SwAlpha_2[number_of_segments];
    double gamma[number_of_segments];
    double gamma0[number_of_segments];
    double beta[number_of_segments];
    double beta0[number_of_segments];
    double gamma0_gamma_min[number_of_segments];
    double gamma0_gamma_plus[number_of_segments];
    double resultant_local_speed[number_of_segments];
    double phi_BEM[number_of_segments];
    double theta_BEM[number_of_segments];
    double cl_cd[number_of_segments];
    double K1[number_of_segments];
    double K2[number_of_segments];
    double EddyMach_calc[number_of_segments];
    double dist_obs[number_of_segments];
    double D_starred_S[number_of_segments];
    double D_starred_P[number_of_segments];
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
    double theta[number_of_segments];
    double phi[number_of_segments];
    double calc_int_a[number_of_segments];
    double phi_e[number_of_segments];
    double theta_e[number_of_segments];
    double L[number_of_segments];
    double Dh[number_of_segments];
    double Dl[number_of_segments];
    double r_rt[number_of_segments];
    double r_e[number_of_segments];
    double r_1[number_of_segments];
    double c_1[number_of_segments];
    double r_0[number_of_segments];
    double c_0[number_of_segments];
    double local_twist[number_of_segments];

    double DStarXFoilS[number_of_segments];
    double DStarXFoilP[number_of_segments];

    double sp_OASPL_alpha=0;
    double splog_OASPL_alpha=0;
    double sp_OASPL_S=0;
    double splog_OASPL_S=0;
    double sp_OASPL_P=0;
    double splog_OASPL_P=0;
    double sp_OASPL=0;
    double splog_OASPL=0;
    double sp_dBA=0;
    double splog_dBA=0;
    double sp_dBB=0;
    double splog_dBB=0;
    double sp_dBC=0;
    double splog_dBC=0;

    double r_R0  =  0.05; double c_R0 = 0.05500;
    double r_R1  =  0.25; double c_R1 = 0.07500;
    double r_R2  =  1.00; double c_R2 = 0.02000;

    double z=m_parameter.TSRtd;
    QString str= QString::number(z, 'f', 1);

        stream << "Tip Speed Ratio: " << str << endl;
        stream << endl;

        stream << qSetFieldWidth(14)  <<
                  "Section"  << ";" <<
                  "Radius [m]"  << ";" <<
                  "r/R"  << ";" <<
                  "Chord [m]" << ";" <<
                  "Wind Speed [m/s]" << ";" <<
                  "Rotational Speed [rpm]" << ";" <<
                  "Tip Speed Ratio" << ";" <<
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
            if (i==(number_of_segments-1)){axial_ind_fact_n[i] = bdata->m_a_axial.value(i);}else {axial_ind_fact_n[i] = bdata->m_a_axial.value(i+1);}

            if (i<number_of_segments/2) {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact[i]);}
            else {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact_n[i]);}

            tangential_speed[i] = omega*bdata->m_pos.value(i)*(1.+bdata->m_a_tangential.value(i));
            resultant_local_speed[i] = qSqrt(pow(axial_velocity[i],2)+pow(tangential_speed[i],2));
            chord[i] = bdata->m_c_local.value(i);
            Reynolds[i] = bdata->m_Reynolds.value(i);

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
            alpha_error[i]=qFabs(alpha_polar[i]-alpha_BEM[i])/alpha_BEM[i]*100.;

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
if (i==(number_of_segments-1)){L[i]=0;}else{
L[i]=bdata->m_pos.value(i+1)-bdata->m_pos.value(i);}

//Calculate the Switching Angle
SwAlpha_1[i]=23.43*Mach[i]+4.651;
SwAlpha_2[i]=12.5;

if (SwAlpha_1[i]<SwAlpha_2[i]){SwAlpha[i]=SwAlpha_1[i];}
else {SwAlpha[i]=SwAlpha_2[i];}

double EddyMach = m_parameter.eddyConvectionMach;

Dh[i]=(2.*pow(sin(qDegreesToRadians(theta_e[i]/2.)),2)*pow(sin(qDegreesToRadians(phi_e[i])),2))/pow(1+Mach[i]*cos(qDegreesToRadians(theta_e[i]))*(1.+(Mach[i]-Mach[i]*EddyMach)*cos(qDegreesToRadians(phi_e[i]))),2);

Dl[i]=(2.*pow(sin(qDegreesToRadians(theta_e[i])),2)*pow(sin(qDegreesToRadians(phi_e[i])),2))/pow((1+(Mach[i]*cos(qDegreesToRadians(theta_e[i])))),4);

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
FoilPolarDlg *pFoilPolarDlg = (FoilPolarDlg *) g_mainFrame->m_pctrlXDirectWidget;

    double TopTrip=pFoilPolarDlg->m_XTopTr;
    double BotTrip=pFoilPolarDlg->m_XBotTr;

    if(((TopTrip<=0) & (TopTrip>=0)) & ((BotTrip<=0) & (BotTrip>=0))) {
//        natural transition
    D_starred_S[i]=D_starred_HT_S[i];
    D_starred_P[i]=D_starred_HT_P[i];
}
else {
        //heavy tripping
            D_starred_S[i]=DStarXFoilS[i];
            D_starred_P[i]=DStarXFoilP[i];
        }}
else if (m_parameter.dstar_type==1){
    D_starred_S[i]=D_starred_N_S[i];
    D_starred_P[i]=D_starred_N_P[i];
}
else if (m_parameter.dstar_type==2){
    //user
    NoiseParameter *pNoiseParameter = (NoiseParameter *) g_mainFrame->m_pSimuWidget;
    D_starred_S[i]=pNoiseParameter->D_starred_S_user[i];
    D_starred_P[i]=pNoiseParameter->D_starred_P_user[i];
        }

double XB;
double YB;
double ZB;

//    re phi_e and theta_e calculation p 77 C_Project_Log_Text_Jan_16.pdf
//    Input X e , Y e , Z e
//    Attribute their respective values to X B , Y B , Z B
    XB=m_parameter.obs_x_pos;
    YB=m_parameter.obs_y_pos;
    ZB=m_parameter.obs_z_pos;

    c_0[i]=bdata->m_c_local.value(i);
    r_0[i]=bdata->m_pos.value(i);

    double hub_radius=pbem->m_pBlade->m_HubRadius;
    double outer_radius=pbem->m_pTData->OuterRadius;
    double blade_radius=(outer_radius-hub_radius);

if(i==(number_of_segments-1)){
c_1[i]=bdata->m_c_local.value(i);
r_1[i]=blade_radius;
}
else
{
c_1[i]=bdata->m_c_local.value(i+1);
r_1[i]=bdata->m_pos.value(i+1);
}

local_twist[i]=theta_BEM[i];

    b[i]=qRadiansToDegrees(qAtan((c_1[i]-c_0[i])/(r_1[i]-r_0[i])));

//    the angle a is the total angle between the Y B Z B blade reference system plane and the local midsection chord line p 75 handout
        a[i]=local_twist[i]+blade_pitch;

    XRS[i]=XB*cos(qDegreesToRadians(a[i]))+YB*sin(qDegreesToRadians(a[i]));
    YRS[i]=-XB*sin(qDegreesToRadians(a[i]))+YB*cos(qDegreesToRadians(a[i]));
    ZRS[i]=ZB-(r_1[i]-r_0[i])/2.;

//    Calculate Y RS − 0.75 ∗ (C r i+1 − C r i )/2
    calc_int_a[i]=(YRS[i]-0.75*(c_1[i]-c_0[i])/2.);

    XRT[i]=XRS[i];
    YRT[i]=cos(qDegreesToRadians(b[i]))*calc_int_a[i]+sin(qDegreesToRadians(b[i]))*ZRS[i];
    ZRT[i]=-sin(qDegreesToRadians(b[i]))*calc_int_a[i]+cos(qDegreesToRadians(b[i]))*ZRS[i];

    r_e[i]=sqrt(pow(XRT[i],2)+pow(YRT[i],2)+pow(ZRT[i],2));
    r_rt[i]=r_e[i];
    theta_e[i]=qRadiansToDegrees(qAtan(ZRT[i]/YRT[i]));
    phi_e[i]=qRadiansToDegrees(qAtan(XRT[i]/ZRT[i]));

    if(i==number_of_segments){
        theta_e[i]=0;
        phi_e[i]=0;
    }
    dist_obs[i]=r_e[i];

    //phi type and theta type fixed 90º or calculated
   phi_rad[i]=qDegreesToRadians(phi[i]);
   theta_rad[i]=qDegreesToRadians(theta[i]);

   alpha_error[i]=qFabs(alpha_polar[i]-alpha_BEM[i])/alpha_BEM[i]*100.;

QString observations_x("");
//new validation for greater Reynolds
if(!(((((alpha[i]<=19.8) & (Mach[i]<0.21)) & (Reynolds[i]>0)) & (Mach[i]>0)))){observations_x.append("1 ");}
if (!(((((m_parameter.Lowson_type!=0) & (Mach[i]<=0.18)) & (Mach[i]>0))) & (Reynolds[i]>0))){observations_x.append("2");}

//uncomment to input data
// if((z<=m_parameter.TSRtd) & (z>=m_parameter.TSRtd)){
        stream << qSetFieldWidth(14)  <<
                  (i+1) << ";" <<
                  bdata->m_pos.value(i) << ";" <<
                                        r_R[i] << ";" <<
                                        chord[i] << ";" <<
                                        bdata->m_Windspeed.value(i) << ";" <<
                                        m_parameter.rot_speed << ";" <<
                                        m_parameter.TSRtd << ";" <<
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
}
    }
         z=z+ldelta;
}
     stream << endl;
     stream << "[1]1 - Out of range for BPM method. 2 - Out of range for LE method."<< endl;
}

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
    float simulation_time;
    float number_time_steps;

    if(g_windFieldStore.size() == 0){
        simulation_time = 0;
        number_time_steps = 0;
    }else{
        simulation_time = g_windFieldModule->getShownWindField()->getSimulationTime();
        number_time_steps = g_windFieldModule->getShownWindField()->getNumberOfTimesteps();}
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
    case P::TSR_check:
        if(set) m_parameter.TSR_check = value.toBool();
        else {value=m_parameter.TSR_check;}
break;

    case P::shear_check:
        if(set) m_parameter.shear_check = value.toBool();
        else {value=m_parameter.shear_check;}
break;

    case P::u_wind_speed_check:
        if(set) m_parameter.u_wind_speed_check = value.toBool();
        else {value=m_parameter.u_wind_speed_check;}
break;

    case P::rot_speed_check:
        if(set) {m_parameter.rot_speed_check = value.toBool();}
else {value=m_parameter.rot_speed_check;}
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

    case P::dstar_type:
        if(set) m_parameter.dstar_type = value.toInt();
        else {value = m_parameter.dstar_type;}break;

    case P::phi_type:
        if(set) m_parameter.phi_type = value.toInt();
        else {value = m_parameter.phi_type;}break;

    case P::theta_type:
        if(set) m_parameter.theta_type = value.toInt();
        else {value = m_parameter.theta_type;}break;

    case P::obs_x_pos:
        if(set) m_parameter.obs_x_pos = value.toDouble();
        else {value = m_parameter.obs_x_pos;}break;

    case P::obs_y_pos:
        if(set) m_parameter.obs_y_pos = value.toDouble();
        else {value = m_parameter.obs_y_pos;}break;

    case P::obs_z_pos:
        if(set) m_parameter.obs_z_pos = value.toDouble();
        else {value = m_parameter.obs_z_pos;}break;

    case P::obs_x_pos_rotor:
         if(set) m_parameter.obs_x_pos_rotor = value.toDouble();
        else {value = m_parameter.obs_x_pos_rotor;}break;

    case P::obs_y_pos_rotor:
        if(set) m_parameter.obs_y_pos_rotor = value.toDouble();
        else {value = m_parameter.obs_y_pos_rotor;}break;

    case P::obs_z_pos_rotor:
        if(set) m_parameter.obs_z_pos_rotor = value.toDouble();
        else {value = m_parameter.obs_z_pos_rotor;}break;

    case P::shear_roughness:
        if(set) m_parameter.shear_roughness = value.toDouble();
        else {value=m_parameter.shear_roughness;}
break;

    case P::shear_height:
        if(set) m_parameter.shear_height = value.toDouble();
        else {value=m_parameter.shear_height;}
break;

    case P::shear_speed:
        if(set) m_parameter.shear_speed = value.toDouble();
        else {value=m_parameter.shear_speed;}
break;

    case P::tower_height:
        if(set) m_parameter.tower_height = value.toDouble();
        else {value = m_parameter.tower_height;}
        break;

    case P::yaw_angle:
        if(set) m_parameter.yaw_angle = value.toDouble();
        else {value = m_parameter.yaw_angle;}
        break;

    case P::initial_azimuth:
        if(set) m_parameter.initial_azimuth = value.toDouble();
        else {value = m_parameter.initial_azimuth;}
        break;

    case P::tower_to_hub_distance:
        if(set) m_parameter.tower_to_hub_distance = value.toDouble();
        else {value = m_parameter.tower_to_hub_distance;}
        break;

    case P::Lowson_type:
        if(set) {m_parameter.Lowson_type = value.toInt();}
else {value=m_parameter.Lowson_type;}
        break;

    case P::qs3DSim:
        if(set) {m_parameter.qs3DSim = value.toInt();}
else {value=m_parameter.qs3DSim;}
        break;

    case P::state_ss_us:
        if(set) {m_parameter.state_ss_us = value.toInt();}
else {value=m_parameter.state_ss_us;}
        break;

    case P::rotation_type:
       if(set) m_parameter.rotation_type = value.toInt();
       else {value = m_parameter.rotation_type;}
       break;

    case P::number_loops:
       if(set) m_parameter.number_loops = value.toInt();
       else {value = m_parameter.number_loops;}
       break;

    case P::time:
       if(set) m_parameter.time = value.toDouble();
       else {value = m_parameter.time;}
       break;

    case P::timesteps:
       if(set) m_parameter.timesteps = value.toInt();
       else {value = m_parameter.timesteps;}
       break;

    case P::anglesteps:
        if(set) m_parameter.anglesteps = value.toInt();
        else {value = m_parameter.anglesteps;}
        break;
    }
// Sara

    return (set ? QVariant() : value);
}
