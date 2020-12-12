#include "NoiseSimulation.h"

//Sara
#include "NoiseCalculation.h"
#include "NoiseParameter.h"
#include "../XUnsteadyBEM/WindFieldModule.h"
#include "../XLLT/QLLTSimulation.h"
#include "../XLLT/QLLTCreatorDialog.h"
#include "../XDirect/FoilPolarDlg.h"
#include "../XBEM/BData.h"
#include "../XBEM/TData.h"
#include "../XBEM/TBEMData.h"
#include "../XDMS/DData.h"
#include "../XBEM/Blade.h"
#include "../XDirect/XDirect.h"
#include "NoiseException.h"
#include "NoiseOpPoint.h"
#include "NoiseCreatorDialog.h"
#include <cmath>
#include <sstream>
#include <string>
#include <QtMath>
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

//Sara changed all endl; for Qt:: endl;

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
        if (m_parameter.qs3DSim==0){//2d
            const int index = getAvailableVariables().indexOf(i == 0 ? xAxis : yAxis);
            QVector<double> *vector = (i == 0 ? &xVector : &yVector);
        switch (index) {
            case 0: *vector = NoiseCalculation::CENTRAL_BAND_FREQUENCY; break;
            case 1: *vector = m_calculation.SPLadB()[opPointIndex]; zeroY = true; break;
            case 2: *vector = m_calculation.SPLsdB()[opPointIndex]; zeroY = true; break;
            case 3: *vector = m_calculation.SPLpdB()[opPointIndex]; zeroY = true; break;
            case 4: *vector = m_calculation.SPL_LEdB()[opPointIndex]; zeroY = true; break; //Alexandre MOD Sara
            case 5: *vector = m_calculation.SPL_LEdBAW()[opPointIndex]; zeroY = true; break;
            case 6: *vector = m_calculation.SPL_LEdBBW()[opPointIndex]; zeroY = true; break;
            case 7: *vector = m_calculation.SPL_LEdBCW()[opPointIndex]; zeroY = true; break;
            case 8: *vector = m_calculation.SPL_LBLVSdB()[opPointIndex]; zeroY = true; break;
            case 9: *vector = m_calculation.SPL_bluntdB()[opPointIndex]; zeroY = true; break;
            case 10: *vector = m_calculation.SPLdB()[opPointIndex]; zeroY = true; break;
            case 11: *vector = m_calculation.SPLdBAW()[opPointIndex]; zeroY = true; break;
            case 12: *vector = m_calculation.SPLdBBW()[opPointIndex]; zeroY = true; break;
            case 13: *vector = m_calculation.SPLdBCW()[opPointIndex]; zeroY = true; break;
            case 14: *vector = m_calculation.SPLadB3d()[opPointIndex]; zeroY = true; break;
            case 15: *vector = m_calculation.SPLsdB3d()[opPointIndex]; zeroY = true; break;
            case 16: *vector = m_calculation.SPLpdB3d()[opPointIndex]; zeroY = true; break;
            case 17: *vector = m_calculation.SPL_LEdB3d()[opPointIndex]; zeroY = true; break;
            case 18: *vector = m_calculation.SPL_LEdBAW3d()[opPointIndex]; zeroY = true; break;
            case 19: *vector = m_calculation.SPL_LEdBBW3d()[opPointIndex]; zeroY = true; break;
            case 20: *vector = m_calculation.SPL_LEdBCW3d()[opPointIndex]; zeroY = true; break;
            case 21: *vector = m_calculation.SPL_LBLVSdB3d()[opPointIndex]; zeroY = true; break;
            case 22: *vector = m_calculation.SPL_bluntdB3d()[opPointIndex]; zeroY = true; break;
            case 23: *vector = m_calculation.SPL_tipvortexdB3d()[opPointIndex]; zeroY = true; break;
            case 24: *vector = m_calculation.SPLdB3d()[opPointIndex]; zeroY = true; break;
            case 25: *vector = m_calculation.SPLdBAW3d()[opPointIndex]; zeroY = true; break;
            case 26: *vector = m_calculation.SPLdBBW3d()[opPointIndex]; zeroY = true; break;
            case 27: *vector = m_calculation.SPLdBCW3d()[opPointIndex]; zeroY = true; break;

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
            case 5: *vector = m_calculation.SPL_LEdBAW()[opPointIndex]; zeroY = true; break;
            case 6: *vector = m_calculation.SPL_LEdBBW()[opPointIndex]; zeroY = true; break;
            case 7: *vector = m_calculation.SPL_LEdBCW()[opPointIndex]; zeroY = true; break;
            case 8: *vector = m_calculation.SPL_LBLVSdB()[opPointIndex]; zeroY = true; break; //Sara
            case 9: *vector = m_calculation.SPL_bluntdB()[opPointIndex]; zeroY = true; break; //Sara
            case 10: *vector = m_calculation.SPL_tipvortexdB()[opPointIndex]; zeroY = true; break; //Sara
            case 11: *vector = m_calculation.SPLdB()[opPointIndex]; zeroY = true; break;
            case 12: *vector = m_calculation.SPLdBAW()[opPointIndex]; zeroY = true; break;
            case 13: *vector = m_calculation.SPLdBBW()[opPointIndex]; zeroY = true; break;
            case 14: *vector = m_calculation.SPLdBCW()[opPointIndex]; zeroY = true; break;
            case 15: *vector = m_calculation.SPLadB3d_final()[opPointIndex]; zeroY = true; break;
            case 16: *vector = m_calculation.SPLsdB3d_final()[opPointIndex]; zeroY = true; break;
            case 17: *vector = m_calculation.SPLpdB3d_final()[opPointIndex]; zeroY = true; break;
            case 18: *vector = m_calculation.SPL_LEdB3d_final()[opPointIndex]; zeroY = true; break;
            case 19: *vector = m_calculation.SPL_LEdBAW3d_final()[opPointIndex]; zeroY = true; break;
            case 20: *vector = m_calculation.SPL_LEdBBW3d_final()[opPointIndex]; zeroY = true; break;
            case 21: *vector = m_calculation.SPL_LEdBCW3d_final()[opPointIndex]; zeroY = true; break;
            case 22: *vector = m_calculation.SPL_LBLVSdB3d_final()[opPointIndex]; zeroY = true; break;
            case 23: *vector = m_calculation.SPL_bluntdB3d_final()[opPointIndex]; zeroY = true; break;
            case 24: *vector = m_calculation.SPL_tipvortexdB3d_final()[opPointIndex]; zeroY = true; break;
            case 25: *vector = m_calculation.SPLdB3d_final()[opPointIndex]; zeroY = true; break;
            case 26: *vector = m_calculation.SPLdBAW3d_final()[opPointIndex]; zeroY = true; break;
            case 27: *vector = m_calculation.SPLdBBW3d_final()[opPointIndex]; zeroY = true; break;
            case 28: *vector = m_calculation.SPLdBCW3d_final()[opPointIndex]; zeroY = true; break;
            case 29: *vector = m_calculation.SPLadB3d()[opPointIndex];zeroY = true; break;
            case 30: *vector = m_calculation.SPLsdB3d()[opPointIndex]; zeroY = true; break;
            case 31: *vector = m_calculation.SPLpdB3d()[opPointIndex]; zeroY = true; break;
            case 32: *vector = m_calculation.SPL_LEdB3d()[opPointIndex]; zeroY = true; break;
            case 33: *vector = m_calculation.SPL_LEdBAW3d()[opPointIndex]; zeroY = true; break;
            case 34: *vector = m_calculation.SPL_LEdBBW3d()[opPointIndex]; zeroY = true; break;
            case 35: *vector = m_calculation.SPL_LEdBCW3d()[opPointIndex]; zeroY = true; break;
            case 36: *vector = m_calculation.SPL_LBLVSdB3d()[opPointIndex]; zeroY = true; break;
            case 37: *vector = m_calculation.SPL_bluntdB3d()[opPointIndex]; zeroY = true; break;
            case 38: *vector = m_calculation.SPL_tipvortexdB3d()[opPointIndex]; zeroY = true; break;
            case 39: *vector = m_calculation.SPLdB3d()[opPointIndex]; zeroY = true; break;

            default: return nullptr;
            }
        }

        if (m_parameter.qs3DSim==2){//rotor
        const int index = getAvailableVariables_rotor().indexOf(i == 0 ? xAxis : yAxis);
        QVector<double> *vector = (i == 0 ? &xVector : &yVector);
        switch (index) {
            case 0: *vector = NoiseCalculation::CENTRAL_BAND_FREQUENCY; break;
            case 1: *vector = m_calculation.SPLadB()[opPointIndex]; zeroY = true; break;
            case 2: *vector = m_calculation.SPLsdB()[opPointIndex]; zeroY = true; break;
            case 3: *vector = m_calculation.SPLpdB()[opPointIndex]; zeroY = true; break;
            case 4: *vector = m_calculation.SPL_LEdB()[opPointIndex]; zeroY = true; break; //Alexandre MOD Sara
            case 5: *vector = m_calculation.SPL_LEdBAW()[opPointIndex]; zeroY = true; break;
            case 6: *vector = m_calculation.SPL_LEdBBW()[opPointIndex]; zeroY = true; break;
            case 7: *vector = m_calculation.SPL_LEdBCW()[opPointIndex]; zeroY = true; break;
            case 8: *vector = m_calculation.SPL_LBLVSdB()[opPointIndex]; zeroY = true; break; //Sara
            case 9: *vector = m_calculation.SPL_bluntdB()[opPointIndex]; zeroY = true; break; //Sara
            case 10: *vector = m_calculation.SPL_tipvortexdB()[opPointIndex]; zeroY = true; break; //Sara
            case 11: *vector = m_calculation.SPLdB()[opPointIndex];  zeroY = true; break;//Sara
            case 12: *vector = m_calculation.SPLdBAW()[opPointIndex]; zeroY = true; break;
            case 13: *vector = m_calculation.SPLdBBW()[opPointIndex]; zeroY = true; break;
            case 14: *vector = m_calculation.SPLdBCW()[opPointIndex]; zeroY = true; break;
            case 15: *vector = m_calculation.SPLadB3d_final()[opPointIndex]; zeroY = true; break;
            case 16: *vector = m_calculation.SPLsdB3d_final()[opPointIndex]; zeroY = true; break;
            case 17: *vector = m_calculation.SPLpdB3d_final()[opPointIndex]; zeroY = true; break;
            case 18: *vector = m_calculation.SPL_LEdB3d_final()[opPointIndex]; zeroY = true; break;
            case 19: *vector = m_calculation.SPL_LEdBAW3d_final()[opPointIndex]; zeroY = true; break;
            case 20: *vector = m_calculation.SPL_LEdBBW3d_final()[opPointIndex]; zeroY = true; break;
            case 21: *vector = m_calculation.SPL_LEdBCW3d_final()[opPointIndex]; zeroY = true; break;
            case 22: *vector = m_calculation.SPL_LBLVSdB3d_final()[opPointIndex]; zeroY = true; break;
            case 23: *vector = m_calculation.SPL_bluntdB3d_final()[opPointIndex]; zeroY = true; break;
            case 24: *vector = m_calculation.SPL_tipvortexdB3d_final()[opPointIndex]; zeroY = true; break;
            case 25: *vector = m_calculation.SPLdB3d_final()[opPointIndex]; zeroY = true; break;
            case 26: *vector = m_calculation.SPLdBAW3d_final()[opPointIndex]; zeroY = true; break;
            case 27: *vector = m_calculation.SPLdBBW3d_final()[opPointIndex]; zeroY = true; break;
            case 28: *vector = m_calculation.SPLdBCW3d_final()[opPointIndex]; zeroY = true; break;
            case 29: *vector = m_calculation.SPLadB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 30: *vector = m_calculation.SPLsdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 31: *vector = m_calculation.SPLpdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 32: *vector = m_calculation.SPL_LEdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 33: *vector = m_calculation.SPL_LEdBAW3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 34: *vector = m_calculation.SPL_LEdBBW3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 35: *vector = m_calculation.SPL_LEdBCW3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 36: *vector = m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 37: *vector = m_calculation.SPL_bluntdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 38: *vector = m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 39: *vector = m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 40: *vector = m_calculation.SPLdB3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 41: *vector = m_calculation.SPLdBAW3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 42: *vector = m_calculation.SPLdBBW3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 43: *vector = m_calculation.SPLdBCW3d_final_rotor_loops()[opPointIndex]; zeroY = true; break;
            case 44: *vector = m_calculation.SPLadB3d()[opPointIndex];zeroY = true; break;
            case 45: *vector = m_calculation.SPLsdB3d()[opPointIndex];zeroY = true; break;
            case 46: *vector = m_calculation.SPLpdB3d()[opPointIndex];zeroY = true; break;
            case 47: *vector = m_calculation.SPL_LEdB3d()[opPointIndex]; zeroY = true; break;
            case 48: *vector = m_calculation.SPL_LEdBAW3d()[opPointIndex];zeroY = true; break;
            case 49: *vector = m_calculation.SPL_LEdBBW3d()[opPointIndex]; zeroY = true; break;
            case 50: *vector = m_calculation.SPL_LEdBCW3d()[opPointIndex]; zeroY = true; break;
            case 51: *vector = m_calculation.SPL_LBLVSdB3d()[opPointIndex]; zeroY = true; break;
            case 52: *vector = m_calculation.SPL_bluntdB3d()[opPointIndex]; zeroY = true; break;
            case 53: *vector = m_calculation.SPL_tipvortexdB3d()[opPointIndex]; zeroY = true; break;
            case 54: *vector = m_calculation.SPLdB3d()[opPointIndex]; zeroY = true; break;

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
    variables << "Freq [Hz]" << "SPL_alpha" << "SPL_S" << "SPL_P" << "SPL_LE (dB)" << "SPL_LE (dB(A))" << "SPL_LE (dB(B))" << "SPL_LE (dB(C))" << "SPL_LBL_VS (dB)" <<  "SPL_blunt (dB)" <<  "SPL_tipvortex (dB)" << "SPL (dB)" << "SPL (dB(A))" << "SPL (dB(B))" << "SPL (dB(C))" << "SPL_alpha_blade[qs3D]" << "SPL_S_blade[qs3D]" << "SPL_P_blade[qs3D]" << "SPL_LE_blade[qs3D] (dB)" << "SPL_LE_blade[qs3D] (dB(A))" << "SPL_LE_blade[qs3D] (dB(B))" << "SPL_LE_blade[qs3D] (dB(C))" << "SPL_LBL_VS_blade[qs3D] (dB)" << "SPL_blunt_blade[qs3D] (dB)" << "SPL_tipvortex_blade[qs3D] (dB)" << "SPL_blade[qs3D] (dB)" << "SPL_blade[qs3D] (dB(A))" << "SPL_blade[qs3D] (dB(B))" << "SPL_blade[qs3D] (dB(C))" << "SPL_alpha_rotor[qs3D]" << "SPL_S_rotor[qs3D]" << "SPL_P_rotor[qs3D]" << "SPL_LE_rotor[qs3D] (dB)" << "SPL_LE_rotor[qs3D] (dB(A))" << "SPL_LE_rotor[qs3D] (dB(B))" << "SPL_LE_rotor[qs3D] (dB(C))" << "SPL_LBL_VS_rotor[qs3D] (dB)" << "SPL_blunt_rotor[qs3D] (dB)" << "SPL_tipvortex_rotor[qs3D] (dB)" << "SPL_rotor[qs3D] (dB)" << "SPL_rotor[qs3D] (dB(A))" << "SPL_rotor[qs3D] (dB(B))" << "SPL_rotor[qs3D] (dB(C))" <<  "SPL_alpha_multi[qs3D]" << "SPL_S_multi[qs3D]" << "SPL_P_multi[qs3D]" << "SPL_LE_multi[qs3D] (dB)" << "SPL_LE_multi[qs3D] (dB(A))" << "SPL_LE_multi[qs3D] (dB(B))" << "SPL_LE_multi[qs3D] (dB(C))" << "SPL_LBL_VS_multi[qs3D] (dB)" << "SPL_blunt_multi[qs3D] (dB)" << "SPL_tipvortex_multi[qs3D] (dB)" << "SPL_multi[qs3D] (dB)"; //Alexandre MOD Sara

    return variables;
}

QStringList NoiseSimulation::getAvailableVariables_blade(NewGraph::GraphType /*graphType*/) {
    QStringList variables;
    // WARNING: when changing any variables list, change newCurve as well!
    variables << "Freq [Hz]" << "SPL_alpha" << "SPL_S" << "SPL_P" << "SPL_LE (dB)" << "SPL_LE (dB(A))" << "SPL_LE (dB(B))" << "SPL_LE (dB(C))" << "SPL_LBL_VS (dB)" <<  "SPL_blunt (dB)" <<  "SPL_tipvortex (dB)" << "SPL (dB)" << "SPL (dB(A))" << "SPL (dB(B))" << "SPL (dB(C))" << "SPL_alpha_blade[qs3D]" << "SPL_S_blade[qs3D]" << "SPL_P_blade[qs3D]" << "SPL_LE_blade[qs3D] (dB)" << "SPL_LE_blade[qs3D] (dB(A))" << "SPL_LE_blade[qs3D] (dB(B))" << "SPL_LE_blade[qs3D] (dB(C))" << "SPL_LBL_VS_blade[qs3D] (dB)" << "SPL_blunt_blade[qs3D] (dB)" << "SPL_tipvortex_blade[qs3D] (dB)" << "SPL_blade[qs3D] (dB)" << "SPL_blade[qs3D] (dB(A))" << "SPL_blade[qs3D] (dB(B))" << "SPL_blade[qs3D] (dB(C))" <<  "SPL_alpha_multi[qs3D]" << "SPL_S_multi[qs3D]" << "SPL_P_multi[qs3D]" << "SPL_LE_multi[qs3D] (dB)" << "SPL_LE_multi[qs3D] (dB(A))" << "SPL_LE_multi[qs3D] (dB(B))" << "SPL_LE_multi[qs3D] (dB(C))" << "SPL_LBL_VS_multi[qs3D] (dB)" << "SPL_blunt_multi[qs3D] (dB)" << "SPL_tipvortex_multi[qs3D] (dB)" << "SPL_multi[qs3D] (dB)"; //Alexandre MOD Sara

    return variables;
}
//Sara end

QStringList NoiseSimulation::getAvailableVariables(NewGraph::GraphType /*graphType*/) {
    QStringList variables;
    // WARNING: when changing any variables list, change newCurve as well!
    variables << "Freq [Hz]" << "SPL_alpha" << "SPL_S" << "SPL_P" << "SPL_LE (dB)" << "SPL_LE (dB(A))" << "SPL_LE (dB(B))" << "SPL_LE (dB(C))" << "SPL_LBL_VS (dB)" <<  "SPL_blunt (dB)" <<  "SPL_tipvortex (dB)" << "SPL (dB)" << "SPL (dB(A))" << "SPL (dB(B))" << "SPL (dB(C))"; //Alexandre MOD Sara

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
        if (g_tbemdataStore.isEmpty()){
            message.prepend("\n - No Turbine BEM Simulation in Database");
        }
        if (g_rotorStore.isEmpty() && g_qbem->m_pBlade == NULL) {
            message.prepend(tr("- No HAWT Blade in Database (for quasi 3D)"));
        }
        if (g_windFieldStore.size() == 0) {
            message.prepend(tr("- No Windfield in Database (for quasi 3D unsteady)"));
        }
        message.prepend("- No Noise Simulation in Database");
            //Sara end
            return message;
    } else {
        return QStringList();
    }
}

//Sara
void NoiseSimulation::pre_simulate() {
    m_calculation.setNoiseParam(&m_parameter);

//Sara
//verify and generate new polars when quasi 3d multiple polars options are active
if ((m_parameter.qs3DSim!=0) & (m_parameter.opPointSource == NoiseParameter::MultiplePolars)){
    m_calculation.onVerifyDeltaandValFor3D();
    if(m_parameter.autopolars_check){loopsReMaalpha();}
}
}
//Sara

void NoiseSimulation::simulate() {
    NoiseCreatorDialog *pNoiseCreatorDialog = (NoiseCreatorDialog *) g_mainFrame->m_pBEM;//Sara
    m_calculation.setNoiseParam(&m_parameter);

//Sara
//verify and generate new polars when quasi 3d multiple polars options are active
if (m_parameter.qs3DSim!=0){
    m_calculation.onVerifyDeltaandValFor3D();
    m_calculation.onVerifyDeltaandValFor3DAlerts();
}

if(m_calculation.simulate_yes_no){
    pNoiseCreatorDialog->OnProgressDlg();
    m_calculation.setInitialValues();
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculate();}

if (m_parameter.qs3DSim!=0){
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculateqs3d_graphics_loops();}
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculateqs3d_blade();}
    }
if(m_parameter.qs3DSim==2){
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){m_calculation.calculateqs3d_rotor();}
    }
if (!pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){pNoiseCreatorDialog->m_progress_dlg->cancel();}
}
//Sara
}

void NoiseSimulation::exportCalculation(QTextStream &stream) {
    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Noise prediction file export" << Qt::endl;
    stream << Qt::endl;
    if(m_parameter.Lowson_type==1){stream << "leading edge model: Von Karman" << Qt::endl;}
    if(m_parameter.Lowson_type==2){stream << "leading edge model: Rapid Distortion" << Qt::endl;}//era m_calculation.m_
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();
    for (int i = 0; i < noiseOpPoints.size(); ++i) {
        stream << qSetFieldWidth(0);
        stream << "Alpha: " << noiseOpPoints[i]->getAlphaDegree() <<
                  ", Re = " << noiseOpPoints[i]->getReynolds() << Qt::endl;
        stream << "SPL_a: " << m_calculation.SPLALOG()[i] << "" << Qt::endl;
        stream << "SPL_s: " << m_calculation.SPLSLOG()[i] << "" << Qt::endl;
        stream << "SPL_p: " << m_calculation.SPLPLOG()[i] << "" << Qt::endl;
        if(m_parameter.Lowson_type!=0){
        stream << "SPL_LE: " << m_calculation.SPLlogLE()[i] << " dB" << Qt::endl;
        }
        if(m_parameter.LBLVS!=0){
        stream << "SPL_LBLVS: " << m_calculation.SPLlogLBLVS()[i] << " dB" << Qt::endl;
        }
        if(m_parameter.blunt_check!=0){
        stream << "SPL_blunt: " << m_calculation.SPLlogblunt()[i] << " dB" << Qt::endl;
        }
        if(m_parameter.tipvortex_check!=0){
        stream << "SPL_tipvortex: " << m_calculation.SPLlogtipvortex()[i] << " dB" << Qt::endl;
        }
        stream << "OASPL: " << m_calculation.OASPL()[i] << " dB" << Qt::endl;
        stream << "OASPL (A): " << m_calculation.OASPLA()[i] << " dB(A)" << Qt::endl;
        stream << "OASPL (B): " << m_calculation.OASPLB()[i] << " dB(B)" << Qt::endl;
        stream << "OASPL (C): " << m_calculation.OASPLC()[i] << " dB(C)" << Qt::endl;
        stream << Qt::endl;
               if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" <<
                  "SPLa" <<
                  "SPLs" <<
                  "SPLp" <<
                  "SPL_LE (dB)" <<
                  "SPL_LBL_VS (dB)" <<
                  "SPL_blunt (dB)" <<
                  "SPL_tipvortex (dB)" <<
                  "SPL (dB)" <<
                  "SPL (dB(A))" <<
                  "SPL (dB(B))" <<
                  "SPL (dB(C))" <<
                  Qt::endl; //Alexandre MOD
               }
               else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
                   stream << qSetFieldWidth(14) <<
                             "Freq [Hz]" <<
                             "SPLa" <<
                             "SPLs" <<
                             "SPLp" <<
                             "SPL_LE (dB)" <<
                             "SPL_LBL_VS (dB)" <<
                             "SPL_blunt (dB)" <<
                             "SPL (dB)" <<
                             "SPL (dB(A))" <<
                             "SPL (dB(B))" <<
                             "SPL (dB(C))" <<
                             Qt::endl; //Alexandre MOD
                          }

else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
                   stream << qSetFieldWidth(14) <<
                             "Freq [Hz]" <<
                             "SPLa" <<
                             "SPLs" <<
                             "SPLp" <<
                             "SPL_LE (dB)" <<
                             "SPL_LBL_VS (dB)" <<
                             "SPL_tipvortex (dB)" <<
                             "SPL (dB)" <<
                             "SPL (dB(A))" <<
                             "SPL (dB(B))" <<
                             "SPL (dB(C))" <<
                             Qt::endl; //Alexandre MOD
                          }

               else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_LE (dB)" <<
                                            "SPL_LBL_VS (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }

else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
                   stream << qSetFieldWidth(14) <<
                             "Freq [Hz]" <<
                             "SPLa" <<
                             "SPLs" <<
                             "SPLp" <<
                             "SPL_LE (dB)" <<
                             "SPL_blunt (dB)" <<
                             "SPL_tipvortex (dB)" <<
                             "SPL (dB)" <<
                             "SPL (dB(A))" <<
                             "SPL (dB(B))" <<
                             "SPL (dB(C))" <<
                             Qt::endl; //Alexandre MOD
                          }

               else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_LE (dB)" <<
                                            "SPL_blunt (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }

               else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_LE (dB)" <<
                                            "SPL_tipvortex (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }

               else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
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
                                            Qt::endl; //Alexandre MOD
                                         }

               else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_blunt (dB)" <<
                                            "SPL_tipvortex (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }
               else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_blunt (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }
               else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_tipvortex (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }
               else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }
               else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_LBL_VS (dB)" <<
                                            "SPL_blunt (dB)" <<
                                            "SPL_tipvortex (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }
               else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_LBL_VS (dB)" <<
                                            "SPL_blunt (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }
               else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_LBL_VS (dB)" <<
                                            "SPL_tipvortex (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }
               else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
                                  stream << qSetFieldWidth(14) <<
                                            "Freq [Hz]" <<
                                            "SPLa" <<
                                            "SPLs" <<
                                            "SPLp" <<
                                            "SPL_LBL_VS (dB)" <<
                                            "SPL (dB)" <<
                                            "SPL (dB(A))" <<
                                            "SPL (dB(B))" <<
                                            "SPL (dB(C))" <<
                                            Qt::endl; //Alexandre MOD
                                         }

        for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

 if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
            stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                      m_calculation.SPLadB()[i][j] <<
                      m_calculation.SPLsdB()[i][j] <<
                      m_calculation.SPLpdB()[i][j] <<
                      m_calculation.SPL_LEdB()[i][j] <<
                      m_calculation.SPL_LBLVSdB()[i][j] <<
                      m_calculation.SPL_bluntdB()[i][j] <<
                      m_calculation.SPL_tipvortexdB()[i][j] <<
                      m_calculation.SPLdB()[i][j] <<
                      m_calculation.SPLdBAW()[i][j] <<
                      m_calculation.SPLdBBW()[i][j] <<
                      m_calculation.SPLdBCW()[i][j] <<
                      Qt::endl; //Alexandre MOD
        }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
            stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                      m_calculation.SPLadB()[i][j] <<
                      m_calculation.SPLsdB()[i][j] <<
                      m_calculation.SPLpdB()[i][j] <<
                      m_calculation.SPL_LEdB()[i][j] <<
                      m_calculation.SPL_LBLVSdB()[i][j] <<
                      m_calculation.SPL_tipvortexdB()[i][j] <<
                      m_calculation.SPLdB()[i][j] <<
                      m_calculation.SPLdBAW()[i][j] <<
                      m_calculation.SPLdBBW()[i][j] <<
                      m_calculation.SPLdBCW()[i][j] <<
                      Qt::endl; //Alexandre MOD
        }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
            stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                      m_calculation.SPLadB()[i][j] <<
                      m_calculation.SPLsdB()[i][j] <<
                      m_calculation.SPLpdB()[i][j] <<
                      m_calculation.SPL_LEdB()[i][j] <<
                      m_calculation.SPL_LBLVSdB()[i][j] <<
                      m_calculation.SPL_tipvortexdB()[i][j] <<
                      m_calculation.SPLdB()[i][j] <<
                      m_calculation.SPLdBAW()[i][j] <<
                      m_calculation.SPLdBBW()[i][j] <<
                      m_calculation.SPLdBCW()[i][j] <<
                      Qt::endl; //Alexandre MOD
        }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
            stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                      m_calculation.SPLadB()[i][j] <<
                      m_calculation.SPLsdB()[i][j] <<
                      m_calculation.SPLpdB()[i][j] <<
                      m_calculation.SPL_LEdB()[i][j] <<
                      m_calculation.SPL_LBLVSdB()[i][j] <<
                      m_calculation.SPLdB()[i][j] <<
                      m_calculation.SPLdBAW()[i][j] <<
                      m_calculation.SPLdBBW()[i][j] <<
                      m_calculation.SPLdBCW()[i][j] <<
                      Qt::endl; //Alexandre MOD
        }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_LEdB()[i][j] <<
                       m_calculation.SPL_LBLVSdB()[i][j] <<
                       m_calculation.SPL_bluntdB()[i][j] <<
                       m_calculation.SPL_tipvortexdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_LEdB()[i][j] <<
                       m_calculation.SPL_LBLVSdB()[i][j] <<
                       m_calculation.SPL_bluntdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_LEdB()[i][j] <<
                       m_calculation.SPL_tipvortexdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_LEdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_LBLVSdB()[i][j] <<
                       m_calculation.SPL_bluntdB()[i][j] <<
                       m_calculation.SPL_tipvortexdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_LBLVSdB()[i][j] <<
                       m_calculation.SPL_bluntdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_LBLVSdB()[i][j] <<
                       m_calculation.SPL_tipvortexdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_LBLVSdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_bluntdB()[i][j] <<
                       m_calculation.SPL_tipvortexdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_bluntdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPL_tipvortexdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
             stream << QString::number(NoiseCalculation::CENTRAL_BAND_FREQUENCY[j], 'f', 0) << //Sara
                       m_calculation.SPLadB()[i][j] <<
                       m_calculation.SPLsdB()[i][j] <<
                       m_calculation.SPLpdB()[i][j] <<
                       m_calculation.SPLdB()[i][j] <<
                       m_calculation.SPLdBAW()[i][j] <<
                       m_calculation.SPLdBBW()[i][j] <<
                       m_calculation.SPLdBCW()[i][j] <<
                       Qt::endl; //Alexandre MOD
         }
        }
        stream << Qt::endl;
        stream << Qt::endl;
    }
    qDeleteAll(noiseOpPoints);
}
//Sara

void NoiseSimulation::exportCalculationqs3DNoise_blade(QTextStream &stream) {
    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Noise prediction file export quasi 3D for blade" << Qt::endl;
    stream << Qt::endl;
    if(m_parameter.Lowson_type==1){stream << "leading edge model: Von Karman" << Qt::endl;}
    if(m_parameter.Lowson_type==2){stream << "leading edge model: Rapid Distortion" << Qt::endl;}
    stream << "Tip Speed Ratio: " << m_parameter.TSRtd << Qt::endl;
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    int number_of_segments = pBEM->m_pBData->m_pos.size();

    stream << Qt::endl;
    stream << Qt::endl;
    stream << "*********************************************" << Qt::endl;
    stream << "Total Values:" << Qt::endl;
    stream << Qt::endl;
    stream << "OASPL: " << m_calculation.Final_qs3d << Qt::endl;
    stream << "SPL alfa: " << m_calculation.Final_qs3d_alpha << Qt::endl;
    stream << "SPL S: " << m_calculation.Final_qs3d_S << Qt::endl;
    stream << "SPL P: " << m_calculation.Final_qs3d_P << Qt::endl;
    if (m_parameter.Lowson_type!=0){stream << "SPL LE: " << m_calculation.Final_qs3d_LE;}
    if (m_parameter.LBLVS!=0){stream << "SPL LBL-VS: " << m_calculation.Final_qs3d_LBLVS;}
    if (m_parameter.blunt_check!=0){stream << "SPL blunt: " << m_calculation.Final_qs3d_blunt;}
    if (m_parameter.tipvortex_check!=0){stream << "SPL tip vortex: " << m_calculation.Final_qs3d_tipvortex;}
    stream << Qt::endl;
    stream << "*********************************************" << Qt::endl;
    stream << Qt::endl;
    stream << Qt::endl;

    for (int i = 0; i < number_of_segments; ++i){
stream << qSetFieldWidth(0);
//        stream << Qt::endl;
stream << "Section: " << (i+1)<<"/"<<number_of_segments << Qt::endl;
        stream << Qt::endl;
if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0)){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" << ";" <<
                  "SPLa" << ";" <<
                  "SPLs" << ";" <<
                  "SPLp" << ";" <<
                  "SPL_LE (dB)" << ";" <<
                  "SPL_LBL_VS (dB)" << ";" <<
                  "SPL_blunt (dB)" << ";" <<
                  "SPL (dB)" << ";" <<
                  "SPL (dB(A))" << ";" <<
                  "SPL (dB(B))" << ";" <<
                  "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
               }
else if ((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0)){
    stream << qSetFieldWidth(14) <<
              "Freq [Hz]" << ";" <<
              "SPLa" << ";" <<
              "SPLs" << ";" <<
              "SPLp" << ";" <<
              "SPL_LE (dB)" << ";" <<
              "SPL_LBL_VS (dB)" << ";" <<
              "SPL (dB)" << ";" <<
              "SPL (dB(A))" << ";" <<
              "SPL (dB(B))" << ";" <<
              "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
           }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0)){
    stream << qSetFieldWidth(14) <<
              "Freq [Hz]" << ";" <<
              "SPLa" << ";" <<
              "SPLs" << ";" <<
              "SPLp" << ";" <<
              "SPL_LE (dB)" << ";" <<
              "SPL_blunt (dB)" << ";" <<
              "SPL (dB)" << ";" <<
              "SPL (dB(A))" << ";" <<
              "SPL (dB(B))" << ";" <<
              "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
           }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
    stream << qSetFieldWidth(14) <<
              "Freq [Hz]" << ";" <<
              "SPLa" << ";" <<
              "SPLs" << ";" <<
              "SPLp" << ";" <<
              "SPL_LE (dB)" << ";" <<
              "SPL (dB)" << ";" <<
              "SPL (dB(A))" << ";" <<
              "SPL (dB(B))" << ";" <<
              "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
           }
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0)){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" << ";" <<
                  "SPLa" << ";" <<
                  "SPLs" << ";" <<
                  "SPLp" << ";" <<
                  "SPL_LBL_VS (dB)" << ";" <<
                  "SPL_blunt (dB)" << ";" <<
                  "SPL (dB)" << ";" <<
                  "SPL (dB(A))" << ";" <<
                  "SPL (dB(B))" << ";" <<
                  "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
               }
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0)){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" << ";" <<
                  "SPLa" << ";" <<
                  "SPLs" << ";" <<
                  "SPLp" << ";" <<
                  "SPL_LBL_VS (dB)" << ";" <<
                  "SPL (dB)" << ";" <<
                  "SPL (dB(A))" << ";" <<
                  "SPL (dB(B))" << ";" <<
                  "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
               }
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0)){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" << ";" <<
                  "SPLa" << ";" <<
                  "SPLs" << ";" <<
                  "SPLp" << ";" <<
                  "SPL_blunt (dB)" << ";" <<
                  "SPL (dB)" << ";" <<
                  "SPL (dB(A))" << ";" <<
                  "SPL (dB(B))" << ";" <<
                  "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
               }
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
        stream << qSetFieldWidth(14) <<
                  "Freq [Hz]" << ";" <<
                  "SPLa" << ";" <<
                  "SPLs" << ";" <<
                  "SPLp" << ";" <<
                  "SPL (dB)" << ";" <<
                  "SPL (dB(A))" << ";" <<
                  "SPL (dB(B))" << ";" <<
                  "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
               }
//int i=0;
        for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

            if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0)){
                       stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                                 m_calculation.SPLadB3d()[i][j] << ";" <<
                                 m_calculation.SPLsdB3d()[i][j] << ";" <<
                                 m_calculation.SPLpdB3d()[i][j] << ";" <<
                                 m_calculation.SPL_LEdB3d()[i][j] << ";" <<
                                 m_calculation.SPL_LBLVSdB3d()[i][j] << ";" <<
                                 m_calculation.SPL_bluntdB3d()[i][j] << ";" <<
                                 m_calculation.SPLdB3d()[i][j] << ";" <<
                                 m_calculation.SPLdBAW3d()[i][j] << ";" <<
                                 m_calculation.SPLdBBW3d()[i][j] << ";" <<
                                 m_calculation.SPLdBCW3d()[i][j] << ";" <<
                                 Qt::endl; //Alexandre MOD
                   }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0)){
                stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                          m_calculation.SPLadB3d()[i][j] << ";" <<
                          m_calculation.SPLsdB3d()[i][j] << ";" <<
                          m_calculation.SPLpdB3d()[i][j] << ";" <<
                          m_calculation.SPL_LEdB3d()[i][j] << ";" <<
                          m_calculation.SPL_LBLVSdB3d()[i][j] << ";" <<
                          m_calculation.SPLdB3d()[i][j] << ";" <<
                          m_calculation.SPLdBAW3d()[i][j] << ";" <<
                          m_calculation.SPLdBBW3d()[i][j] << ";" <<
                          m_calculation.SPLdBCW3d()[i][j] << ";" <<
                          Qt::endl; //Alexandre MOD
            }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0)){
                    stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                              m_calculation.SPLadB3d()[i][j] << ";" <<
                              m_calculation.SPLsdB3d()[i][j] << ";" <<
                              m_calculation.SPLpdB3d()[i][j] << ";" <<
                              m_calculation.SPL_LEdB3d()[i][j] << ";" <<
                              m_calculation.SPL_bluntdB3d()[i][j] << ";" <<
                              m_calculation.SPLdB3d()[i][j] << ";" <<
                              m_calculation.SPLdBAW3d()[i][j] << ";" <<
                              m_calculation.SPLdBBW3d()[i][j] << ";" <<
                              m_calculation.SPLdBCW3d()[i][j] << ";" <<
                              Qt::endl; //Alexandre MOD
                }
            else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
                                stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                                          m_calculation.SPLadB3d()[i][j] << ";" <<
                                          m_calculation.SPLsdB3d()[i][j] << ";" <<
                                          m_calculation.SPLpdB3d()[i][j] << ";" <<
                                          m_calculation.SPL_LEdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdBAW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBBW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBCW3d()[i][j] << ";" <<
                                          Qt::endl; //Alexandre MOD
                            }
            else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0)){
                                stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                                          m_calculation.SPLadB3d()[i][j] << ";" <<
                                          m_calculation.SPLsdB3d()[i][j] << ";" <<
                                          m_calculation.SPLpdB3d()[i][j] << ";" <<
                                          m_calculation.SPL_LBLVSdB3d()[i][j] << ";" <<
                                          m_calculation.SPL_bluntdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdBAW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBBW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBCW3d()[i][j] << ";" <<
                                          Qt::endl; //Alexandre MOD
                            }
            else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0)){
                                stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                                          m_calculation.SPLadB3d()[i][j] << ";" <<
                                          m_calculation.SPLsdB3d()[i][j] << ";" <<
                                          m_calculation.SPLpdB3d()[i][j] << ";" <<
                                          m_calculation.SPL_LBLVSdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdBAW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBBW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBCW3d()[i][j] << ";" <<
                                          Qt::endl; //Alexandre MOD
                            }
            else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0)){
                                stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                                          m_calculation.SPLadB3d()[i][j] << ";" <<
                                          m_calculation.SPLsdB3d()[i][j] << ";" <<
                                          m_calculation.SPLpdB3d()[i][j] << ";" <<
                                          m_calculation.SPL_bluntdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdBAW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBBW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBCW3d()[i][j] << ";" <<
                                          Qt::endl; //Alexandre MOD
                            }
            else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
                                stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                                          m_calculation.SPLadB3d()[i][j] << ";" <<
                                          m_calculation.SPLsdB3d()[i][j] << ";" <<
                                          m_calculation.SPLpdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdB3d()[i][j] << ";" <<
                                          m_calculation.SPLdBAW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBBW3d()[i][j] << ";" <<
                                          m_calculation.SPLdBCW3d()[i][j] << ";" <<
                                          Qt::endl; //Alexandre MOD
                            }
        }
        stream << Qt::endl;
        stream << qSetFieldWidth(0);
        stream << "SPL_a: " << m_calculation.SPLALOG3d()[i] << "" << Qt::endl;
        stream << "SPL_s: " << m_calculation.SPLSLOG3d()[i] << "" << Qt::endl;
        stream << "SPL_p: " << m_calculation.SPLPLOG3d()[i] << "" << Qt::endl;
        if(m_parameter.Lowson_type!=0){
        stream << "SPL_LE: " << m_calculation.SPLlogLE3d()[i] << " dB" << Qt::endl;}
        if(m_parameter.LBLVS!=0){
        stream << "SPL_LBL_VS: " << m_calculation.SPLlogLBLVS3d()[i] << " dB" << Qt::endl;}
        if(m_parameter.blunt_check!=0){
        stream << "SPL_blunt: " << m_calculation.SPLlogblunt3d()[i] << " dB" << Qt::endl;}
        stream << "OASPL: " << m_calculation.OASPL3d()[i] << " dB" << Qt::endl;
        stream << "OASPL (A): " << m_calculation.OASPLA3d()[i] << " dB(A)" << Qt::endl;
        stream << "OASPL (B): " << m_calculation.OASPLB3d()[i] << " dB(B)" << Qt::endl;
        stream << "OASPL (C): " << m_calculation.OASPLC3d()[i] << " dB(C)" << Qt::endl;
        stream << Qt::endl;
        stream << Qt::endl;

if(i==(number_of_segments-1)){
stream << "********** FINAL **********" << Qt::endl;
if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0)){
stream << qSetFieldWidth(14) <<
   "Freq [Hz]" << ";" <<
   "SPLa" << ";" <<
   "SPLs" << ";" <<
   "SPLp" << ";" <<
   "SPL_LE (dB)" << ";" <<
   "SPL_LBL_VS (dB)" << ";" <<
   "SPL_blunt (dB)" << ";" <<
   "SPL (dB)" << ";" <<
   "SPL (dB(A))" << ";" <<
   "SPL (dB(B))" << ";" <<
   "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
}
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0)){
    stream << qSetFieldWidth(14) <<
       "Freq [Hz]" << ";" <<
       "SPLa" << ";" <<
       "SPLs" << ";" <<
       "SPLp" << ";" <<
       "SPL_LE (dB)" << ";" <<
       "SPL_LBL_VS (dB)" << ";" <<
       "SPL (dB)" << ";" <<
       "SPL (dB(A))" << ";" <<
       "SPL (dB(B))" << ";" <<
       "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
    }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0)){
    stream << qSetFieldWidth(14) <<
       "Freq [Hz]" << ";" <<
       "SPLa" << ";" <<
       "SPLs" << ";" <<
       "SPLp" << ";" <<
       "SPL_LE (dB)" << ";" <<
       "SPL_blunt (dB)" << ";" <<
       "SPL (dB)" << ";" <<
       "SPL (dB(A))" << ";" <<
       "SPL (dB(B))" << ";" <<
       "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
    }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
    stream << qSetFieldWidth(14) <<
       "Freq [Hz]" << ";" <<
       "SPLa" << ";" <<
       "SPLs" << ";" <<
       "SPLp" << ";" <<
       "SPL_LE (dB)" << ";" <<
       "SPL (dB)" << ";" <<
       "SPL (dB(A))" << ";" <<
       "SPL (dB(B))" << ";" <<
       "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
    }
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0)){
    stream << qSetFieldWidth(14) <<
       "Freq [Hz]" << ";" <<
       "SPLa" << ";" <<
       "SPLs" << ";" <<
       "SPLp" << ";" <<
       "SPL_LBL_VS (dB)" << ";" <<
       "SPL_blunt (dB)" << ";" <<
       "SPL (dB)" << ";" <<
       "SPL (dB(A))" << ";" <<
       "SPL (dB(B))" << ";" <<
       "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
    }
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0)){
    stream << qSetFieldWidth(14) <<
       "Freq [Hz]" << ";" <<
       "SPLa" << ";" <<
       "SPLs" << ";" <<
       "SPLp" << ";" <<
       "SPL_LBL_VS (dB)" << ";" <<
       "SPL (dB)" << ";" <<
       "SPL (dB(A))" << ";" <<
       "SPL (dB(B))" << ";" <<
       "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
    }
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0)){
stream << qSetFieldWidth(14) <<
   "Freq [Hz]" << ";" <<
   "SPLa" << ";" <<
   "SPLs" << ";" <<
   "SPLp" << ";" <<
   "SPL_blunt (dB)" << ";" <<
   "SPL (dB)" << ";" <<
   "SPL (dB(A))" << ";" <<
   "SPL (dB(B))" << ";" <<
   "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
}
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
stream << qSetFieldWidth(14) <<
   "Freq [Hz]" << ";" <<
   "SPLa" << ";" <<
   "SPLs" << ";" <<
   "SPLp" << ";" <<
   "SPL (dB)" << ";" <<
   "SPL (dB(A))" << ";" <<
   "SPL (dB(B))" << ";" <<
   "SPL (dB(C))" << ";" << Qt::endl; //Alexandre MOD
}
for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0)){
        stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                  m_calculation.SPLadB3d_final()[i][j] << ";" <<
                  m_calculation.SPLsdB3d_final()[i][j] << ";" <<
                  m_calculation.SPLpdB3d_final()[i][j] << ";" <<
                  m_calculation.SPL_LEdB3d_final()[i][j] << ";" <<
                  m_calculation.SPL_LBLVSdB3d_final()[i][j] << ";" <<
                  m_calculation.SPL_bluntdB3d_final()[i][j] << ";" <<
                  m_calculation.SPLdB3d_final()[i][j] << ";" <<
                  m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
                  m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
                  m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
                  Qt::endl; //Alexandre MOD
    }
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0)){
    stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
              m_calculation.SPLadB3d_final()[i][j] << ";" <<
              m_calculation.SPLsdB3d_final()[i][j] << ";" <<
              m_calculation.SPLpdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_LEdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_LBLVSdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
              Qt::endl; //Alexandre MOD
}
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0)){
    stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
              m_calculation.SPLadB3d_final()[i][j] << ";" <<
              m_calculation.SPLsdB3d_final()[i][j] << ";" <<
              m_calculation.SPLpdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_LEdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_bluntdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
              Qt::endl; //Alexandre MOD
}
else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
    stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
              m_calculation.SPLadB3d_final()[i][j] << ";" <<
              m_calculation.SPLsdB3d_final()[i][j] << ";" <<
              m_calculation.SPLpdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_LEdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
              Qt::endl; //Alexandre MOD
}
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0)){
    stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
              m_calculation.SPLadB3d_final()[i][j] << ";" <<
              m_calculation.SPLsdB3d_final()[i][j] << ";" <<
              m_calculation.SPLpdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_LBLVSdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_bluntdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
              Qt::endl; //Alexandre MOD
}
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0)){
    stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
              m_calculation.SPLadB3d_final()[i][j] << ";" <<
              m_calculation.SPLsdB3d_final()[i][j] << ";" <<
              m_calculation.SPLpdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_LBLVSdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
              Qt::endl; //Alexandre MOD
}
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0)){
    stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
              m_calculation.SPLadB3d_final()[i][j] << ";" <<
              m_calculation.SPLsdB3d_final()[i][j] << ";" <<
              m_calculation.SPLpdB3d_final()[i][j] << ";" <<
              m_calculation.SPL_bluntdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
              Qt::endl; //Alexandre MOD
}
else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0)){
    stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
              m_calculation.SPLadB3d_final()[i][j] << ";" <<
              m_calculation.SPLsdB3d_final()[i][j] << ";" <<
              m_calculation.SPLpdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdB3d_final()[i][j] << ";" <<
              m_calculation.SPLdBAW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBBW3d_final()[i][j] << ";" <<
              m_calculation.SPLdBCW3d_final()[i][j] << ";" <<
              Qt::endl; //Alexandre MOD
}
}
stream << Qt::endl;
}}}

void NoiseSimulation::exportCalculationqs3DNoise_rotor(QTextStream &stream) {
    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Noise prediction file export quasi 3D for a moving rotor" << Qt::endl;
    stream << Qt::endl;
    if(m_parameter.Lowson_type==1){stream << "leading edge model: Von Karman" <<Qt::endl;}
    if(m_parameter.Lowson_type==2){stream << "leading edge model: Rapid Distortion" <<Qt::endl;}
    stream << "Tip Speed Ratio: " << m_parameter.TSRtd << Qt::endl;
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    int number_of_segments = pBEM->m_pBData->m_pos.size();

     stream << Qt::endl;
     stream << Qt::endl;
     stream << "*********************************************" << Qt::endl;
     stream << "Total Values:" << Qt::endl;
     stream << Qt::endl;
     stream << "OASPL: " << m_calculation.Final_qs3d_rotor_loops << Qt::endl;
     stream << "SPL alfa: " << m_calculation.Final_qs3d_alpha_rotor_loops << Qt::endl;
     stream << "SPL S: " << m_calculation.Final_qs3d_S_rotor_loops << Qt::endl;
     stream << "SPL P: " << m_calculation.Final_qs3d_P_rotor_loops << Qt::endl;
     if (m_parameter.Lowson_type!=0){stream << "SPL LE: " << m_calculation.Final_qs3d_LE_rotor_loops;}
     if (m_parameter.LBLVS!=0){stream << "SPL LBL-VS: " << m_calculation.Final_qs3d_LBLVS_rotor_loops;}
     if (m_parameter.blunt_check!=0){stream << "SPL blunt: " << m_calculation.Final_qs3d_blunt_rotor_loops;}
     if (m_parameter.tipvortex_check!=0){stream << "SPL tip vortex: " << m_calculation.Final_qs3d_tipvortex_rotor_loops;}
     stream << Qt::endl;
     stream << "*********************************************" << Qt::endl;
     stream << Qt::endl;
     stream << Qt::endl;

     for (int i = 0; i < number_of_segments; ++i){
 stream << qSetFieldWidth(0);

 if(i==(number_of_segments-1)){
// stream << "********** FINAL **********" << Qt::endl;
 if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
 stream << qSetFieldWidth(14) <<
    "Freq [Hz]" << ";" <<
    "SPLa" << ";" <<
    "SPLs" << ";" <<
    "SPLp" << ";" <<
    "SPL_LE (dB)" << ";" <<
    "SPL_LBL_VS (dB)" << ";" <<
    "SPL_blunt (dB)" << ";" <<
    "SPL_tipvortex (dB)" << ";" <<
    "SPL (dB)" << ";" <<
    "SPL (dB(A))" << ";" <<
    "SPL (dB(B))" << ";" <<
    "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
 }
 else  if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LE (dB)" << ";" <<
        "SPL_LBL_VS (dB)" << ";" <<
        "SPL_tipvortex (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LE (dB)" << ";" <<
        "SPL_LBL_VS (dB)" << ";" <<
        "SPL_tipvortex (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LE (dB)" << ";" <<
        "SPL_LBL_VS (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LE (dB)" << ";" <<
        "SPL_blunt (dB)" << ";" <<
        "SPL_tipvortex (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LE (dB)" << ";" <<
        "SPL_blunt (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LE (dB)" << ";" <<
        "SPL_tipvortex (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LE (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LBL_VS (dB)" << ";" <<
        "SPL_blunt (dB)" << ";" <<
        "SPL_tipvortex (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LBL_VS (dB)" << ";" <<
        "SPL_blunt (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LBL_VS (dB)" << ";" <<
        "SPL_tipvortex (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_LBL_VS (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_blunt (dB)" << ";" <<
        "SPL_tipvortex (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_blunt (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL_tipvortex (dB)" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 else  if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
     stream << qSetFieldWidth(14) <<
        "Freq [Hz]" << ";" <<
        "SPLa" << ";" <<
        "SPLs" << ";" <<
        "SPLp" << ";" <<
        "SPL (dB)" << ";" <<
        "SPL (dB(A))" << ";" <<
        "SPL (dB(B))" << ";" <<
        "SPL (dB(C))" << ";" <<Qt::endl; //Alexandre MOD
     }
 for (int j = 0; j < NoiseCalculation::FREQUENCY_TABLE_SIZE; ++j) {

 if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
         stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                   m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
                   Qt::endl; //Alexandre MOD
     }
 if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
         stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
                   m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
                   m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
                   Qt::endl; //Alexandre MOD
     }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS!=0) &  (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type!=0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LEdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS!=0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_LBLVSdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check!=0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check!=0) & (m_parameter.tipvortex_check==0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_bluntdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check!=0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPL_tipvortexdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
 else if((m_parameter.Lowson_type==0) & (m_parameter.LBLVS==0) & (m_parameter.blunt_check==0) & (m_parameter.tipvortex_check==0)){
     stream << NoiseCalculation::CENTRAL_BAND_FREQUENCY[j] << ";" <<
               m_calculation.SPLadB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLsdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLpdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdB3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBAW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBBW3d_final_rotor_loops()[i][j] << ";" <<
               m_calculation.SPLdBCW3d_final_rotor_loops()[i][j] << ";" <<
               Qt::endl; //Alexandre MOD
 }
}
stream << Qt::endl;
}}}

void NoiseSimulation::exportqs3DLog(QTextStream &stream) {
    stream.setRealNumberNotation(QTextStream::FixedNotation);
    stream.setRealNumberPrecision(5);

    stream << "Quasi 3D Noise Log" << Qt::endl;
    stream << Qt::endl;

QList<NoiseOpPoint*> noiseOpPoints = m_parameter.prepareNoiseOpPointList();
    double z=m_parameter.TSRtd;
    QString str= QString::number(z, 'f', 1);

        stream << "Tip Speed Ratio: " << str << Qt::endl;
        stream << Qt::endl;

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
                  Qt::endl;

        m_calculation.qs3D_log(stream);

     stream << Qt::endl;
     stream << "[1]1 - Out of TE validation range. 2 - Out of LE validation range. 3 - Reynolds outside of error margin defined. 4 - Mach outside of error margin defined. 5 - Alpha outside of error margin defined."<< Qt::endl;

             stream.setRealNumberNotation(QTextStream::FixedNotation);
             stream.setRealNumberPrecision(5);
             bool blade_values = m_calculation.qs3D_val_blade!=NULL;
             bool rotor_values = m_calculation.qs3D_val_rotor!=NULL;

             if((m_parameter.qs3DSim!=0) & blade_values){
                 stream << Qt::endl;
                 stream << Qt::endl;
                 stream << "Quasi 3D Noise Validation Error Log for Blade Simulation" << Qt::endl;
                 stream << Qt::endl;
                 stream << "Tip Speed Ratio: " << str << Qt::endl;
                 stream << Qt::endl;

                 stream << qSetFieldWidth(14)  <<
                           "Section"  << ";" <<
                           "Reynolds"  << ";" <<
                           "Mach" << ";" <<
                           "AOA" << ";" <<
                           "Error number[1]" <<  ";" <<
                           Qt::endl;

                 stream << qSetFieldWidth(14)  <<
                 m_calculation.qs3D_val_blade <<
                 Qt::endl;

                 stream << Qt::endl;
                 stream << "[1]1 - Out of TE validation range. 2 - Out of LE validation range."<< Qt::endl;
                 stream << Qt::endl;
         }

             if((m_parameter.qs3DSim==2) & rotor_values){
                 stream << Qt::endl;
                 stream << Qt::endl;
                 stream << "Quasi 3D Noise Validation Error Log for Rotor Simulation" << Qt::endl;
                 stream << Qt::endl;
                 stream << "Tip Speed Ratio: " << str << Qt::endl;
                 stream << Qt::endl;

                 stream << qSetFieldWidth(14)  <<
                           "Blade"  << ";" <<
                           "Azimuthal Angle [deg]"  << ";" <<
                           "Section"  << ";" <<
                           "Reynolds"  << ";" <<
                           "Mach" << ";" <<
                           "AOA" << ";" <<
                           "Error number[1]" <<  ";" <<
                           Qt::endl;

                 stream << qSetFieldWidth(14)  <<
                 m_calculation.qs3D_val_rotor <<
                           Qt::endl;

                 stream << Qt::endl;
                 stream << "[1]1 - Out of TE validation range. 2 - Out of LE validation range."<< Qt::endl;
                 stream << Qt::endl;
         }
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
    case P::LBLVS:
        if(set) m_parameter.LBLVS = value.toBool();
        else value = m_parameter.LBLVS; break;

    case P::blunt_check:
        if(set) m_parameter.blunt_check = value.toBool();
        else value = m_parameter.blunt_check; break;

    case P::hblunt_check:
        if(set) m_parameter.hblunt_check = value.toBool();
        else value = m_parameter.hblunt_check; break;

    case P::hblunt:
        if(set) m_parameter.hblunt = value.toDouble();
        else value = m_parameter.hblunt; break;

    case P::tipvortex_check:
        if(set) m_parameter.tipvortex_check = value.toBool();
        else value = m_parameter.tipvortex_check; break;

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

    case P::valRel_TE_check:
        if(set) {m_parameter.valRel_TE_check = value.toBool();}
else {value=m_parameter.valRel_TE_check;}
        break;

    case P::valReu_TE_check:
        if(set) {m_parameter.valReu_TE_check = value.toBool();}
else {value=m_parameter.valReu_TE_check;}
        break;

    case P::valRel_TE:
        if(set) {m_parameter.valRel_TE = value.toDouble();}
else {value=m_parameter.valRel_TE;}
        break;

    case P::valReu_TE:
        if(set) {m_parameter.valReu_TE = value.toDouble();}
else {value=m_parameter.valReu_TE;}
        break;

    case P::valMal_TE_check:
        if(set) {m_parameter.valMal_TE_check = value.toBool();}
else {value=m_parameter.valMal_TE_check;}
        break;

    case P::valMau_TE_check:
        if(set) {m_parameter.valMau_TE_check = value.toBool();}
else {value=m_parameter.valMau_TE_check;}
        break;

    case P::valMal_TE:
        if(set) {m_parameter.valMal_TE = value.toDouble();}
else {value=m_parameter.valMal_TE;}
        break;

    case P::valMau_TE:
        if(set) {m_parameter.valMau_TE = value.toDouble();}
else {value=m_parameter.valMau_TE;}
        break;

    case P::valAOAl_TE_check:
        if(set) {m_parameter.valAOAl_TE_check = value.toBool();}
else {value=m_parameter.valAOAl_TE_check;}
        break;

    case P::valAOAu_TE_check:
        if(set) {m_parameter.valAOAu_TE_check = value.toBool();}
else {value=m_parameter.valAOAu_TE_check;}
        break;

    case P::valAOAl_TE:
        if(set) {m_parameter.valAOAl_TE = value.toDouble();}
else {value=m_parameter.valAOAl_TE;}
        break;

    case P::valAOAu_TE:
        if(set) {m_parameter.valAOAu_TE = value.toDouble();}
else {value=m_parameter.valAOAu_TE;}
        break;

    case P::valRel_LE_check:
        if(set) {m_parameter.valRel_LE_check = value.toBool();}
else {value=m_parameter.valRel_LE_check;}
        break;

    case P::valReu_LE_check:
        if(set) {m_parameter.valReu_LE_check = value.toBool();}
else {value=m_parameter.valReu_LE_check;}
        break;

    case P::valRel_LE:
        if(set) {m_parameter.valRel_LE = value.toDouble();}
else {value=m_parameter.valRel_LE;}
        break;

    case P::valReu_LE:
        if(set) {m_parameter.valReu_LE = value.toDouble();}
else {value=m_parameter.valReu_LE;}
        break;

    case P::valMal_LE_check:
        if(set) {m_parameter.valMal_LE_check = value.toBool();}
else {value=m_parameter.valMal_LE_check;}
        break;

    case P::valMau_LE_check:
        if(set) {m_parameter.valMau_LE_check = value.toBool();}
else {value=m_parameter.valMau_LE_check;}
        break;

    case P::autopolars_check:
        if(set) {m_parameter.autopolars_check = value.toBool();}
else {value=m_parameter.autopolars_check;}
        break;

    case P::valMal_LE:
        if(set) {m_parameter.valMal_LE = value.toDouble();}
else {value=m_parameter.valMal_LE;}
        break;

    case P::valMau_LE:
        if(set) {m_parameter.valMau_LE = value.toDouble();}
else {value=m_parameter.valMau_LE;}
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

//Sara
void NoiseSimulation::createPolars(int size){
    QXDirect *pXDirect = (QXDirect*) g_mainFrame->m_pXDirect;

    double Reynolds[size];
    double Mach[size];
    double alpha[size];
    double acrit[size];
    double xbot[size];
    double xtop[size];
    double aspec[size];
    double polar_type[size];
    double alpha_max[size];
    enumPolarType ptype;

//    g_mainFrame->OnXDirect();

    for (int i=0;i<size;++i){
    Reynolds[i]=m_calculation.Reynolds_error_value().at(i);
    Mach[i]=m_calculation.Mach_error_value().at(i);
    alpha[i]=m_calculation.alpha_error_value().at(i);
    acrit[i]=m_calculation.acrit_error().at(i);
    xbot[i]=m_calculation.xbot_error().at(i);
    xtop[i]=m_calculation.xtop_error().at(i);
    aspec[i]=m_calculation.aspec_error().at(i);
    polar_type[i]=m_calculation.polar_type_error().at(i);
    alpha_max[i]=m_calculation.alpha_error_value_max().at(i);

    if(polar_type[i]==1){ptype = FIXEDSPEEDPOLAR;}
    else if(polar_type[i]==2){ptype = FIXEDLIFTPOLAR;}
    else if(polar_type[i]==3){ptype = RUBBERCHORDPOLAR;}
    else if(polar_type[i]==4){ptype = FIXEDAOAPOLAR;}

    double delta = 0.01;

    pXDirect->OnNoiseLoop(acrit[i], xbot[i], xtop[i], Mach[i], Reynolds[i], ptype, aspec[i],alpha[i],alpha_max[i],delta);
    }

//    g_noiseModule->onActivationActionTriggered();
}

void NoiseSimulation::create360Polars(){
    QBEM *pQBEM = (QBEM *) g_mainFrame->m_pBEM;
//    pQBEM->On360View();
    pQBEM->m_pctrlNew360All->click();
}

void NoiseSimulation::loopsReMaalpha(){
if(m_parameter.qs3DSim!=0){
int size = m_calculation.Reynolds_error_value().size();

if(size>0){
createPolars(size);}
}}
//Sara
