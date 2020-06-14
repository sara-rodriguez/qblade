#include "NoiseCalculation.h"

#include <QDebug>

#include "NoiseParameter.h"
#include "../Serializer.h"
#include "NoiseOpPoint.h"
#include "NoiseException.h"
#include "../XBEM/BEM.h" //Sara
#include "../Objects/Polar.h"//Sara

//Sara
#include "../Noise/NoiseParameter.h"
#include "../StorableObject.h"
#include "../Graph/ShowAsGraphInterface.h"
#include "../ParameterObject.h"
#include "../XDirect/FoilPolarDlg.h"
#include "../XDirect/XDirect.h"
#include "../XBEM/BData.h"
#include "../XBEM/Blade.h"
#include "../MainFrame.h"
#include "../XUnsteadyBEM/WindFieldModule.h"
#include "../XUnsteadyBEM/WindField.h"
#include "../XLLT/QLLTSimulation.h"
#include "../XLLT/QLLTCreatorDialog.h"
#include "../XUnsteadyBEM/FASTSimulation.h"
#include "../Noise/NoiseCreatorDialog.h"
#include "../Noise/NoiseSimulation.h"
#include "../MainFrame.h"
#include "../XWidgets.h"
#include <QtMath>
#include <cmath>
//Sara

const double NoiseCalculation::AWeighting[] = {-70.4, -63.4, -56.7,  -50.5, -44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1, -13.4, -10.9,  -8.6,  -6.6,  -4.8,  -3.2,  -1.9,  -0.8, 0.0,   0.6,   1.0,   1.2,   1.3,   1.2,   1.0,   0.5,-0.1,-1.1,  -2.5,  -4.3,  -6.6,-  9.3};
const double NoiseCalculation::BWeighting[] = {38.2, -33.3, -28.3, -24.2, -20.4, -17.1, -14.2, -11.6,  -9.3,  -7.4,  -5.6,  -4.2, -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.1,   0.0, 0.0,   0.0,   0.0,  -0.1,  -0.2,  -0.4,  -0.7,  -1.2,-1.9,  -2.9,  -4.3,  -6.1,  -8.4, -11.1};
const double NoiseCalculation::CWeighting[] = {-14.3, -11.2, -8.5, -6.2, -4.4,  -3.0,  -2.0,  -1.3,  -0.8,  -0.5,  -0.3,  -0.2, -0.1,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 0.0,   0.0,  -0.1,  -0.2,  -0.3,  -0.5,  -0.8,  -1.3,-2.0,-3.0,  -4.4,  -6.2,  -8.5, -11.2};
const QVector<double> NoiseCalculation::CENTRAL_BAND_FREQUENCY ({10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160,200, 250, 315, 400, 500, 630, 800, 1000,1250, 1600, 2000, 2500, 3150, 4000, 5000,6300, 8000, 10000, 12500, 16000, 20000});

NoiseCalculation::NoiseCalculation() {
    m_parameter = nullptr; //Alexandre MOD - changed NULL for nullptr
    m_CalcSeparatedFlow = false;
    m_CalcPressureSide = false;
    m_CalcSuctionSide = false;
}

void NoiseCalculation::serialize() {
    g_serializer.readOrWriteDoubleVector2D(&m_SPLadB);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLsdB);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLpdB);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdB);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBAW);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBBW);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBCW);

    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdB); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBAW); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBBW); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBCW); //Alexandre MOD

    g_serializer.readOrWriteDoubleVector1D(&m_OASPL);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLA);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLB);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLC);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLALOG);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLSLOG);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLPLOG);

//Sara

//    multi graphs
    g_serializer.readOrWriteDoubleVector2D(&m_SPLadB3d);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLsdB3d);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLpdB3d);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdB3d);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBAW3d);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBBW3d);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBCW3d);

    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdB3d); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBAW3d); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBBW3d); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBCW3d); //Alexandre MOD

//    blade
    g_serializer.readOrWriteDoubleVector2D(&m_SPLadB3d_final);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLsdB3d_final);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLpdB3d_final);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdB3d_final);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBAW3d_final);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBBW3d_final);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBCW3d_final);

    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdB3d_final); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBAW3d_final); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBBW3d_final); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBCW3d_final); //Alexandre MOD

//    rotor
    g_serializer.readOrWriteDoubleVector2D(&m_SPLadB3d_rotor);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLsdB3d_rotor);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLpdB3d_rotor);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdB3d_rotor);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBAW3d_rotor);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBBW3d_rotor);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBCW3d_rotor);

    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdB3d_rotor); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBAW3d_rotor); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBBW3d_rotor); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBCW3d_rotor); //Alexandre MOD

    //    rotor loops
    g_serializer.readOrWriteDoubleVector2D(&m_SPLadB3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLsdB3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLpdB3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdB3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBAW3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBBW3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBCW3d_rotor_loops);

    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdB3d_rotor_loops); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBAW3d_rotor_loops); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBBW3d_rotor_loops); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBCW3d_rotor_loops); //Alexandre MOD

    //    final rotor loops
    g_serializer.readOrWriteDoubleVector2D(&m_SPLadB3d_final_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLsdB3d_final_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLpdB3d_final_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdB3d_final_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBAW3d_final_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBBW3d_final_rotor_loops);
    g_serializer.readOrWriteDoubleVector2D(&m_SPLdBCW3d_final_rotor_loops);

    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdB3d_final_rotor_loops); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBAW3d_final_rotor_loops); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBBW3d_final_rotor_loops); //Alexandre MOD
    g_serializer.readOrWriteDoubleVector2D(&m_SPL_LEdBCW3d_final_rotor_loops); //Alexandre MOD

//4d
    g_serializer.readOrWriteDoubleVector4D(&m_SPLadB3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLsdB3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLpdB3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLdB3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLdBAW3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLdBBW3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLdBCW3d_4d);

    g_serializer.readOrWriteDoubleVector4D(&m_SPL_LEdB3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPL_LEdBAW3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPL_LEdBBW3d_4d);
    g_serializer.readOrWriteDoubleVector4D(&m_SPL_LEdBCW3d_4d);

    g_serializer.readOrWriteDoubleVector4D(&m_SPLadB3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLsdB3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLpdB3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLdB3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLdBAW3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLdBBW3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPLdBCW3d_4d_blade);

    g_serializer.readOrWriteDoubleVector4D(&m_SPL_LEdB3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPL_LEdBAW3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPL_LEdBBW3d_4d_blade);
    g_serializer.readOrWriteDoubleVector4D(&m_SPL_LEdBCW3d_4d_blade);

//OASPLs
    g_serializer.readOrWriteDoubleVector1D(&m_OASPL3d);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLA3d);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLB3d);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLC3d);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLALOG3d);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLSLOG3d);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLPLOG3d);

    g_serializer.readOrWriteDoubleVector1D(&m_OASPL3d_rotor);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLA3d_rotor);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLB3d_rotor);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLC3d_rotor);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLALOG3d_rotor);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLSLOG3d_rotor);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLPLOG3d_rotor);

    g_serializer.readOrWriteDoubleVector1D(&m_OASPL3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLA3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLB3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector1D(&m_OASPLC3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLALOG3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLSLOG3d_rotor_loops);
    g_serializer.readOrWriteDoubleVector1D(&m_SPLPLOG3d_rotor_loops);

//Sara
}

double NoiseCalculation::getDStarInterpolated(bool top,NoiseOpPoint * nop) {
    bool upDownFind = false;
    double chordUpStream = 0;
    double chordDownStream = 0;
    double dStarUpStream = 0;
    double dStarDownStream = 0;

    //For positive alpha use TopSide else BottomSide
    //int side = top ? 1 : 2;
    int side = top ? 2 : 1;
    int nside = top ? nop->getNSide2() : nop->getNSide1();

    double currentChord = 0;
    double currentDStar = 0;
    double previousChord = 0;
    double previousDStar = 0;

    //Find closest station assuming crescent order on chordStation
    for (int i = 2; i <= nside; ++i) {
        currentChord = nop->getXValue(i, side);
        currentDStar = nop->getDstrAt(i, side);
        previousChord = i == 0 ? currentChord : nop->getXValue(i-1, side);
        previousDStar = i == 0 ? currentDStar: nop->getDstrAt(i-1, side);

        //qDebug() << "i: " << i << " - " << ccur;
//        qDebug() << "currentChord" << nop->getXValue(i, side);
//        qDebug() << "currentDStar" << nop->getDstrAt(i, side);
//        qDebug() << nop->getXValue(i, side);
//        qDebug() << nop->getDstrAt(i, side);
//        qDebug() << "";

        if (currentChord > m_parameter->dStarChordStation) {
            chordUpStream = previousChord;
            chordDownStream = currentChord;
            dStarUpStream = previousDStar;
            dStarDownStream = currentDStar;

//			qDebug() << "Chord UpStream: " << chordUpStream;
//			qDebug() << "Chord DownStream: " << chordDownStream;
//			qDebug() << "D* UpStream: " << dStarUpStream;
//			qDebug() << "D* DownStream: " << dStarDownStream;

            upDownFind = true;
            break;
        }
    }

    if (!upDownFind) {
        qWarning() << "Can not find upstream and downstream. D* Interpolated will be zero ! D* ChordStation target: "
                   << m_parameter->dStarChordStation << " - Last found X ( "<<currentChord<<" ) D* ("
                   << currentDStar <<")";
        throw NoiseException(NoiseException::EXPT_DSTAR_NOT_FOUND, "There is no data to interpolate D* from, at the "
                             "specified chord station");
    }

//    qDebug() << ((dStarUpStream-dStarDownStream) * (m_parameter->dStarChordStation-chordDownStream) /                                       (chordUpStream-chordDownStream)) + dStarDownStream;

    return ((dStarUpStream-dStarDownStream) * (m_parameter->dStarChordStation-chordDownStream)/(chordUpStream-chordDownStream)) + dStarDownStream;
}

//Sara
double NoiseCalculation::getDStarInterpolated3d(bool top,double chord,NoiseOpPoint * nop) {
    bool upDownFind = false;
    double chordUpStream = 0;
    double chordDownStream = 0;
    double dStarUpStream = 0;
    double dStarDownStream = 0;

    //For positive alpha use TopSide else BottomSide
    //int side = top ? 1 : 2;
    int side = top ? 2 : 1;
    int nside = top ? nop->getNSide2() : nop->getNSide1();

    double currentChord = 0;
    double currentDStar = 0;
    double previousChord = 0;
    double previousDStar = 0;

    //Find closest station assuming crescent order on chordStation
    for (int i = 2; i <= nside; ++i) {
        currentChord = nop->getXValue(i, side);
        currentDStar = nop->getDstrAt(i, side);
        previousChord = i == 0 ? currentChord : nop->getXValue(i-1, side);
        previousDStar = i == 0 ? currentDStar: nop->getDstrAt(i-1, side);

        //qDebug() << "i: " << i << " - " << ccur;
//        qDebug() << "currentChord" << nop->getXValue(i, side);
//        qDebug() << "currentDStar" << nop->getDstrAt(i, side);
//        qDebug() << nop->getXValue(i, side);
//        qDebug() << nop->getDstrAt(i, side);
//        qDebug() << "";

        if (currentChord > chord) {
            chordUpStream = previousChord;
            chordDownStream = currentChord;
            dStarUpStream = previousDStar;
            dStarDownStream = currentDStar;

//			qDebug() << "Chord UpStream: " << chordUpStream;
//			qDebug() << "Chord DownStream: " << chordDownStream;
//			qDebug() << "D* UpStream: " << dStarUpStream;
//			qDebug() << "D* DownStream: " << dStarDownStream;

            upDownFind = true;
            break;
        }
    }

    if (!upDownFind) {
        qWarning() << "Can not find upstream and downstream. D* Interpolated will be zero ! D* ChordStation target: "
                   << chord << " - Last found X ( "<<currentChord<<" ) D* ("
                   << currentDStar <<")";
        throw NoiseException(NoiseException::EXPT_DSTAR_NOT_FOUND, "There is no data to interpolate D* from, at the "
                             "specified chord station");
    }

//    qDebug() << ((dStarUpStream-dStarDownStream) * (m_parameter->dStarChordStation-chordDownStream) /                                       (chordUpStream-chordDownStream)) + dStarDownStream;

    return ((dStarUpStream-dStarDownStream) * (chord-chordDownStream) /
            (chordUpStream-chordDownStream)) + dStarDownStream;
}

double NoiseCalculation::getDL() {
    //Dh, Low Freq. Directivity Factor
    if (m_parameter->lowFreq) {
        return (2 * (pow(sin(qDegreesToRadians(m_parameter->directivityGreek)),2)) * (pow(sin(qDegreesToRadians(m_parameter->directivityPhi)),2))) /
                pow((1+m_parameter->originalMach * cos(qDegreesToRadians(m_parameter->directivityGreek))),4); //Sara
    } else {
        return 1;
    }
}

double NoiseCalculation::getDH() {
    //Dh, High Freq. Directivity Factor
    if (m_parameter->lowFreq)
    {return (2 * pow((sin((qDegreesToRadians(m_parameter->directivityGreek)/2))),2) * pow((sin(qDegreesToRadians(m_parameter->directivityPhi))),2)) /
                    ((1+m_parameter->originalMach * cos(qDegreesToRadians(m_parameter->directivityGreek))) *
                     (pow((1+(m_parameter->originalMach-m_EddyMachNumber) * cos(qDegreesToRadians(m_parameter->directivityPhi))),2))); //Sara
     } else {
        return 1;

    }
}

double NoiseCalculation::getSt1() {
    return 0.02 * pow(m_parameter->originalMach, -0.6);
}

double NoiseCalculation::getSt2(NoiseOpPoint* nop) {
    const double st1 = getSt1();
    double st2;

    if (nop->getAlphaDegreeAbsolute() < 1.33) {
        st2 = st1;
    } else if (nop->getAlphaDegreeAbsolute() > 12.5) {
        st2 = st1*4.72;
    } else {
        st2 = st1 * pow(10,0.0054* pow((nop->getAlphaDegreeAbsolute()-1.33),2));
    }

    return st2;
}

double NoiseCalculation::getBPMThickness(NoiseOpPoint *nop, AirfoilSide as) {
    double dStarCF = 0;
    double dStarCT = 0;
    double dStar = 0;
    double bpm = 0;

    if (m_parameter->transition == NoiseParameter::FullyTurbulent) {
        if(nop->getReynolds() > 300000) {
            dStarCF = pow(10,(3.411-1.5397 * log10(nop->getReynolds())+0.1059 * pow(log10(nop->getReynolds()),2)));
        } else {
            dStarCF = 0.0601 * (pow(nop->getReynolds(),-0.114));
        }
        dStar = dStarCF * m_parameter->originalChordLength;
//		qDebug() << "BPM FullyTurbulent dStarCF: " << dStarCF;
    } else {
        dStarCT = pow(10,(3.0187-1.5397*log10(nop->getReynolds())+0.1059* pow(log10(nop->getReynolds()),2)));
//		qDebug() << "BPM TransitionFlow dStarCT: " << dStarCT;
        dStar = dStarCT * m_parameter->originalChordLength;
    }

    if (nop->getAlphaDegreeAbsolute() == 0.) { //Sara
        bpm = dStar;
    } else {
        double corFactor = 0;

        if (as == PressureSide) {
            corFactor = pow(10,(-0.0432*nop->getAlphaDegreeAbsolute()+0.00113*pow(nop->getAlphaDegreeAbsolute(),2)));
            bpm = corFactor * dStar;
        } else {
            if (m_parameter->transition == NoiseParameter::FullyTurbulent) {
                if (nop->getAlphaDegreeAbsolute() < 5) {
                    corFactor = pow(10,(0.0679 * nop->getAlphaDegreeAbsolute()));
                } else if (nop->getAlphaDegreeAbsolute() >= 5 && nop->getAlphaDegreeAbsolute() <= 12.5) {
                    corFactor = 0.381*(pow(10,(0.1516*nop->getAlphaDegreeAbsolute())));
                } else {
                    corFactor = 14.296*(pow(10,(0.0258*nop->getAlphaDegreeAbsolute())));
                }
            } else {
                if (nop->getAlphaDegreeAbsolute() < 7.5) {
                    corFactor = pow(10,(0.0679*nop->getAlphaDegreeAbsolute()));
                } else if (nop->getAlphaDegreeAbsolute() >= 7.5 && nop->getAlphaDegreeAbsolute() <= 12.5) {
                    corFactor = 0.0162*(pow(10,(0.3066*nop->getAlphaDegreeAbsolute())));
                } else {
                    corFactor = 54.42*(pow(10,(0.0258*nop->getAlphaDegreeAbsolute())));
                }
            }
            bpm = corFactor * dStar;
        }
    }
//	qDebug() << "BPM D*: " << bpm;

    return bpm;
}

double NoiseCalculation::getK1(NoiseOpPoint* nop) {
    double k1;

    if (nop->getReynolds() < 247000) {
        k1 = -4.31 * log10(nop->getReynolds()) + 156.3;
    } else if (nop->getReynolds() > 800000) {
        k1 = 128.5;
    } else {
        k1 = -9 * log10(nop->getReynolds()) + 181.6;
    }

    return k1;
}

void NoiseCalculation::preCalcA1(NoiseOpPoint* nop) {
    double m_A1ChordBasedReynolds;
    double m_A1Ao;
    double m_A1AMin;
    double m_A1AMax;
//    double m_A1FirstTerm;  // NM was not used...

//	m_A1FirstTerm = 10 * log10(pow(m_parameter->m_OriginalMach, 5) * m_parameter->m_WettedLength * getDL() *
//							   m_DStarFinalS / pow(m_parameter->m_DistanceObsever, 2));
    m_A1ChordBasedReynolds = nop->getReynolds() * 3;

    if(m_A1ChordBasedReynolds < 95200){
        m_A1Ao = 0.57;
    }else if(m_A1ChordBasedReynolds > 857000){
        m_A1Ao = 1.13;
    }else{
        m_A1Ao = (-0.000000000000957*(pow((m_A1ChordBasedReynolds-857000),2))+1.13);
    }

    if(m_A1Ao < 0.204){
        m_A1AMin = (sqrt(67.552-886.788*pow(m_A1Ao,2))-8.219);
    }else if(m_A1Ao > 0.244){
        m_A1AMin = (-142.795*pow(m_A1Ao,3)+103.656*pow(m_A1Ao,2)-57.757*m_A1Ao+6.006);
    }else{
        m_A1AMin = (-32.665*m_A1Ao+3.981);
    }

    if(m_A1Ao < 0.13){
        m_A1AMax = (sqrt(67.552-886.788*pow(m_A1Ao,2))-8.219);
    }else if(m_A1Ao > 0.321){
        m_A1AMax = (-4.669*pow(m_A1Ao,3)+3.491*pow(m_A1Ao,2)-16.699*m_A1Ao+1.149);
    }else{
        m_A1AMax = (-15.901*m_A1Ao+1.098);
    }

    m_A1Ar = (-20-m_A1AMin)/(m_A1AMax-m_A1AMin);

//	qDebug() << "A1 ChordBasedReynolds " << m_A1ChordBasedReynolds;
//	qDebug() << "A1 Ao " << m_A1Ao;
//	qDebug() << "A1 aMin " << m_A1AMin;
//	qDebug() << "A1 aMax " << m_A1AMax;
//	qDebug() << "A1 Ar " << m_A1Ar;
}

void NoiseCalculation::preCalcSPLa(NoiseOpPoint* nop) {
//    qDebug() << "---> SPLa CALCULATION";

    m_SplaFirstTerm=0;

    //If angle is smaller than the switching Angle
    //use dH and B else use dL and A1
    if (nop->getAlphaDegree() <= m_SwAlpha) {
//		qDebug() << "SPLa dH: " << getDH();
        m_SplaFirstTerm = 10 * log10(pow(m_parameter->originalMach, 5) * m_parameter->wettedLength * getDH() *
                                     m_DStarFinalS / pow(m_parameter->distanceObsever, 2));

        if(nop->getReynolds() < 95200){
            m_SplaBo = 0.3;
        }else if(nop->getReynolds() > 857000){
            m_SplaBo = 0.56;
        }else{
            m_SplaBo = (-0.000000000000448*( pow((nop->getReynolds()-857000),2)+0.56));
        }

        if(m_SplaBo < 0.13){
            m_SplaBMin = (sqrt(16.888-886.788*pow(m_SplaBo,2))-4.109);
        }else if(m_SplaBo > 0.145){
            m_SplaBMin = (-817.81*pow(m_SplaBo,3)+335.21*pow(m_SplaBo,2)-135.024*m_SplaBo+10.619);
        }else{
            m_SplaBMin = (-83.607*m_SplaBo+8.138);
        }

        if(m_SplaBo < 0.1){
            m_SplaBMax = (sqrt(16.888-886.788*pow(m_SplaBo,2))-4.109);
        }else if(m_SplaBo > 0.187){
            m_SplaBMax = (-80.541*pow(m_SplaBo,3)+44.174*pow(m_SplaBo,2)-39.381*m_SplaBo+2.344);
        }else{
            m_SplaBMax = (-31.33*m_SplaBo+1.854);
        }

        m_SplaBr = (-20-m_SplaBMin)/(m_SplaBMax-m_SplaBMin);

//		qDebug() << "SPLA Bo " << m_SplaBo;
//		qDebug() << "SPLA bMin " << m_SplaBMin;
//		qDebug() << "SPLA bMax " << m_SplaBMax;
//		qDebug() << "SPLA Br " << m_SplaBr;
    } else {
//		qDebug() << "SPLa dL: " << getDL();
        m_SplaFirstTerm = 10 * log10(pow(m_parameter->originalMach, 5) * m_parameter->wettedLength *
                                     getDL() * m_DStarFinalS / pow(m_parameter->distanceObsever, 2));
        m_ChordBasedReynolds = nop->getReynolds() * 3;

        if(m_ChordBasedReynolds < 95200){
            m_SplaAo = 0.57;
        }else if(m_ChordBasedReynolds > 857000){
            m_SplaAo = 1.13;
        }else{
            m_SplaAo = (-0.000000000000957*(pow((m_ChordBasedReynolds-857000),2))+1.13);
        }

        if(m_SplaAo < 0.204){
            m_SplaAMin = (sqrt(67.552-886.788*pow(m_SplaAo,2))-8.219);
        }else if(m_SplaAo > 0.244){
            m_SplaAMin = (-142.795*pow(m_SplaAo,3)+103.656*pow(m_SplaAo,2)-57.757*m_SplaAo+6.006);
        }else{
            m_SplaAMin = (-32.665*m_SplaAo+3.981);
        }

        if(m_SplaAo < 0.13){
            m_SplaAMax = (sqrt(67.552-886.788*pow(m_SplaAo,2))-8.219);
        }else if(m_SplaAo > 0.321){
            m_SplaAMax = (-4.669*pow(m_SplaAo,3)+3.491*pow(m_SplaAo,2)-16.699*m_SplaAo+1.149);
        }else{
            m_SplaAMax = (-15.901*m_SplaAo+1.098);
        }

        m_SplaAr = (-20-m_SplaAMin)/(m_SplaAMax-m_SplaAMin);

//		qDebug() << "SPLA ChordBasedReynolds " << m_ChordBasedReynolds;
//		qDebug() << "SPLA Ao " << m_SplaAo;
//		qDebug() << "SPLA aMin " << m_SplaAMin;
//		qDebug() << "SPLA aMax " << m_SplaAMax;
//		qDebug() << "SPLA Ar " << m_SplaAr;
    }

    m_SplaSt1 = getSt1();
    m_SplaSt2 = getSt2(nop);

    m_SplaGamma = 27.094 * m_parameter->originalMach + 3.32;
    m_SplaBeta = 72.65*m_parameter->originalMach+10.74;
    m_SplaBetaZero = -34.19*m_parameter->originalMach-13.82;

    m_SplaK1 = getK1(nop);

    m_SplaK2 = 0;
    if(nop->getAlphaDegreeAbsolute() < (m_SwAlpha1-m_SplaGamma)){
        m_SplaK2 = m_SplaK1 - 1000;
    }else if(nop->getAlphaDegreeAbsolute() > (m_SwAlpha1+m_SplaGamma)){
        m_SplaK2 = m_SplaK1 - 12;
    }else{
        m_SplaK2 = m_SplaK1 + (sqrt( pow(m_SplaBeta,2)-pow((m_SplaBeta/m_SplaGamma),2) *
                                     pow((nop->getAlphaDegreeAbsolute()-m_SwAlpha1),2) ) + m_SplaBetaZero);
    }

//	qDebug() << "SPLa firstTerm: " << m_SplaFirstTerm;
//	qDebug() << "SPLa st1: " << m_SplaSt1;
//	qDebug() << "SPLa st2: " << m_SplaSt2;
//	qDebug() << "SPLa gamma: " << m_SplaGamma;
//	qDebug() << "SPLa gamma_zero: " << m_SwAlpha1;
//	qDebug() << "SPLa beta: " << m_SplaBeta;
//	qDebug() << "SPLa betaZero: " << m_SplaBetaZero;
//	qDebug() << "SPLa reynolds: " << nop->Reynolds();
//	qDebug() << "SPLa k1: " << m_SplaK1;
//	qDebug() << "SPLa k2: " << m_SplaK2;
}

void NoiseCalculation::preCalcSPLs(NoiseOpPoint *nop) {
//	qDebug() << "---> SPLs CALCULATION";
    m_SplsFirstTerm = 0;
    m_SplsFirstTerm = 10 * log10(pow(m_parameter->originalMach,5) * m_parameter->wettedLength * getDH() *
                                 m_DStarFinalS / pow(m_parameter->distanceObsever,2));

//	qDebug() << "SPLs dH: " << getDH();
//	qDebug() << "SPLs firstTerm: " << m_SplsFirstTerm;

    m_SplsSt1 = getSt1();
    m_SplsSt2 = getSt2(nop);
    m_SplsK1 = getK1(nop);
    m_SplsSt1Bar = (m_SplsSt1+m_SplsSt2)/2;
    m_SplsK13 = m_SplsK1-3;

//	qDebug() << "SPLs st1: " << m_SplsSt1;
//	qDebug() << "SPLs st2: " << m_SplsSt2;
//	qDebug() << "SPLs st1Bar: " << m_SplsSt1Bar;
//	qDebug() << "SPLs k1: " << m_SplsK1;
//	qDebug() << "SPLs k1-3: " << m_SplsK13;
}

void NoiseCalculation::preCalcSPLp(NoiseOpPoint *nop) {
//    qDebug() << "---> SPLp CALCULATION";
    m_SplpFirstTerm=0;
    m_SplpFirstTerm = 10 * log10(pow(m_parameter->originalMach,5) * m_parameter->wettedLength * getDH() *
                                 m_DStarFinalP / pow(m_parameter->distanceObsever,2));

//    qDebug() << "SPLp dH: " << getDH();
//    qDebug() << "SPLp firstTerm: " << m_SplpFirstTerm;

    m_SplpSt1 = getSt1();
    m_SplpK1 = getK1(nop);
    m_SplpK13 = m_SplpK1-3;
    m_SplpDeltaK1 = 0;

    //Sara
    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    double rho=pBEM->dlg_rho;

    m_ReynoldsBasedDisplacement = rho * m_parameter->originalVelocity * m_DStarFinalP / 0.0000178;
    // Sara

    if (m_ReynoldsBasedDisplacement > 5000) {
        m_SplpDeltaK1 = 0;
    } else {
        m_SplpDeltaK1 = nop->getAlphaDegreeAbsolute() * (1.43*log10(m_ReynoldsBasedDisplacement)-5.29);
    }

//	qDebug() << "Reynolds Based Displacement: " << m_ReynoldsBasedDisplacement;
//	qDebug() << "SPLp st1: " << m_SplpSt1;
//	qDebug() << "SPLp k1: " << m_SplpK1;
//	qDebug() << "SPLp k1-3: " << m_SplpK13;
//	qDebug() << "SPLp DeltaK1: " << m_SplpDeltaK1;
}

void NoiseCalculation::calcSPLa(double alpha, int posOpPoint, int posFreq) {
    double splDb = 0;

    //If angle is smaller than the switching Angle
    //use dH and B else use dL and A1
    if (alpha <= m_SwAlpha ) {
        double sts = CENTRAL_BAND_FREQUENCY[posFreq]*m_DStarFinalS/m_parameter->originalVelocity;
        double b = fabs(log10(sts/m_SplaSt2));
        double bMin = 0;
        double bMax = 0;
        double db = 0;

        if(b < 0.13){
            bMin = (sqrt(16.888-886.788*pow(b,2))-4.109);
        }else if(b > 0.145){
            bMin = (-817.81*pow(b,3)+335.21*pow(b,2)-135.024*b+10.619);
        }else{
            bMin = (-83.607*b+8.138);
        }

        if(b < 0.1){
            bMax = (sqrt(16.888-886.788*pow(b,2))-4.109);
        }else if(b > 0.187){

            bMax = (-80.541*pow(b,3)+44.174*pow(b,2)-39.381*b+2.344);
        }else{
            bMax = (-83.607*b+8.138);
        }
        db =bMin+ m_SplaBr *(bMax-bMin);
        splDb = m_SplaFirstTerm +db+m_SplaK2;

        m_SPLadB[posOpPoint][posFreq] = splDb;
        m_SPLadBAW[posOpPoint][posFreq] = splDb + AWeighting[posFreq];
        m_SPLadBBW[posOpPoint][posFreq] = splDb + BWeighting[posFreq];
        m_SPLadBCW[posOpPoint][posFreq] = splDb + CWeighting[posFreq];

//		qDebug() << "SPLa -> sts("<< sts <<")\t"<< "b("<< b <<")\t"<< "bMin("<< bMin <<")\t"<< "bMax("<< bMax
//				 << ")\t"<< "db("<< db <<")\t"<< "splDb("<< splDb <<")\t"<< "splDb-AW("
//				 << m_SPLadBAW[posOpPoint][posFreq] <<")\t"<< "splDb-BW("<< m_SPLadBBW[posOpPoint][posFreq] <<")\t"
//				 << "splDb-CW("<< m_SPLadBCW[posOpPoint][posFreq] <<")\t";
    } else {
        double sts = CENTRAL_BAND_FREQUENCY[posFreq]*m_DStarFinalS/m_parameter->originalVelocity;
        double a = fabs(log10(sts/m_SplaSt2));
        double aMin = 0;
        double aMax = 0;
        double a1 = 0;

        if(a<0.204){
            aMin = sqrt(67.552-886.788*pow(a,2));
        }else if(a > 0.244){
            aMin = -142.795*pow(a,3)+103.656*pow(a,2)-57.757*a+6.006;
        }else{
            aMin = -32.665*a+3.981;
        }

        if(a<0.13){
            aMax = (sqrt(67.552-886.788*pow(a,2))-8.219);
        }else if(a > 0.321){
            aMax = (-4.669*pow(a,3)+3.491*pow(a,2)-16.699*a+1.149);
        }else{
            aMax = (-15.901*a+1.098);
        }

        a1 =aMin+ m_SplaAr *(aMax-aMin);
        splDb = m_SplaFirstTerm +a1+m_SplaK2;

        m_SPLadB[posOpPoint][posFreq] = splDb;
        m_SPLadBAW[posOpPoint][posFreq] = splDb + AWeighting[posFreq];
        m_SPLadBBW[posOpPoint][posFreq] = splDb + BWeighting[posFreq];
        m_SPLadBCW[posOpPoint][posFreq] = splDb + CWeighting[posFreq];

//		qDebug() << "SPLa -> sts("<< sts <<")\t"<< "a("<< a <<")\t"<< "aMin("<< aMin <<")\t"<< "aMax("<< aMax
//				 << ")\t"<< "a1("<< a1 <<")\t"<< "splDb("<< splDb <<")\t"<< "splDb-AW("
//				 << m_SPLadBAW[posOpPoint][posFreq] <<")\t"<< "splDb-BW("<< m_SPLadBBW[posOpPoint][posFreq] <<")\t"
//				 << "splDb-CW("<< m_SPLadBCW[posOpPoint][posFreq] <<")\t";
    }
}

void NoiseCalculation::calcSPLs(int posOpPoint,int posFreq) {
    double sts = CENTRAL_BAND_FREQUENCY[posFreq]*m_DStarFinalS/m_parameter->originalVelocity;
    double a = fabs(log10(sts/m_SplsSt1Bar));
    double aMin = 0;
    double aMax = 0;
    double a1 = 0;
    double splDb = 0;

    //If angle is bigger than the switching Angle
    //or suction side is mandatory
    if(m_CalcSuctionSide){

        if(a<0.204){
            aMin = sqrt(67.552-886.788*pow(a,2)-8.219);
        }else if(a > 0.244){
            aMin = -142.795*pow(a,3)+103.656*pow(a,2)-57.757*a+6.006;
        }else{
            aMin = -32.665*a+3.981;
        }

        if(a<0.13){
            aMax = (sqrt(67.552-886.788*pow(a,2))-8.219);
        }else if(a > 0.321){
            aMax = (-4.669*pow(a,3)+3.491*pow(a,2)-16.699*a+1.149);
        }else{
            aMax = (-15.901*a+1.098);
        }

        a1 =aMin+ m_A1Ar *(aMax-aMin);
        splDb = m_SplsFirstTerm +a1+m_SplsK13;

    } else {
        splDb = std::numeric_limits<long>::lowest();
    }

    m_SPLsdB[posOpPoint][posFreq] = splDb;
    m_SPLsdBAW[posOpPoint][posFreq] = splDb + AWeighting[posFreq];
    m_SPLsdBBW[posOpPoint][posFreq] = splDb + BWeighting[posFreq];
    m_SPLsdBCW[posOpPoint][posFreq] = splDb + CWeighting[posFreq];

//	qDebug() << "SPLs -> sts("<< sts <<")\t"<< "a("<< a <<")\t"<< "aMin("<< aMin <<")\t"<< "aMax("<< aMax <<")\t"
//			 << "a1("<< a1 <<")\t" << "splDb("<< splDb <<")\t"<< "splDb-AW("<< m_SPLsdBAW[posOpPoint][posFreq]
//			 << ")\t"<< "splDb-BW("<< m_SPLsdBBW[posOpPoint][posFreq] <<")\t"<< "splDb-CW("
//			 << m_SPLsdBCW[posOpPoint][posFreq] <<")\t";
}

void NoiseCalculation::calcSPLp(int posOpPoint,int posFreq) {
    double stp = CENTRAL_BAND_FREQUENCY[posFreq]*m_DStarFinalP /m_parameter->originalVelocity;
    double a = fabs(log10(stp/m_SplpSt1));
    double aMin = 0;
    double aMax = 0;
    double a1 = 0;
    double splDb = 0;

    //If angle is bigger than the switching Angle
    //or pressure side is mandatory
    if(m_CalcPressureSide){

        if(a<0.204){
            aMin = sqrt(67.552-886.788*pow(a,2)-8.219);
        }else if(a > 0.244){
            aMin = -142.795*pow(a,3)+103.656*pow(a,2)-57.757*a+6.006;
        }else{
            aMin = -32.665*a+3.981;
        }

        if(a<0.13){
            aMax = (sqrt(67.552-886.788*pow(a,2))-8.219);
        }else if(a > 0.321){
            aMax = (-4.669*pow(a,3)+3.491*pow(a,2)-16.699*a+1.149);
        }else{
            aMax = (-15.901*a+1.098);
        }

        a1 =aMin+ m_A1Ar *(aMax-aMin);
        splDb = m_SplpFirstTerm +a1+m_SplpK13+m_SplpDeltaK1;

    } else {
        splDb = std::numeric_limits<long>::lowest();
    }

    m_SPLpdB[posOpPoint][posFreq] = splDb;
    m_SPLpdBAW[posOpPoint][posFreq] = splDb + AWeighting[posFreq];
    m_SPLpdBBW[posOpPoint][posFreq] = splDb + BWeighting[posFreq];
    m_SPLpdBCW[posOpPoint][posFreq] = splDb + CWeighting[posFreq];

//	qDebug() << "SPLp -> stp("<< stp <<")\t"<< "a("<< a <<")\t"<< "aMin("<< aMin <<")\t"<< "aMax("<< aMax <<")\t"
//			 << "a1("<< a1 <<")\t" << "splDb("<< splDb <<")\t"<< "splDb-AW("<< m_SPLpdBAW[posOpPoint][posFreq]
//			 << ")\t"<< "splDb-BW("<< m_SPLpdBBW[posOpPoint][posFreq] <<")\t"<< "splDb-CW("
//			 << m_SPLpdBCW[posOpPoint][posFreq] <<")\t";
}

//Alexandre MOD

void NoiseCalculation::LECalc(int posOpPoint,int posFreq) {
    double aux_SPL_LEdB=0;
    double Aux1=0;
    double Aux4=0;
    double Aux5=0;
    c_const=0;
    d_const = 0;

    m_SPL_LEdB[posOpPoint][posFreq] = 0;
    m_SPL_LEdBAW[posOpPoint][posFreq] = 0;
    m_SPL_LEdBBW[posOpPoint][posFreq] = 0;
    m_SPL_LEdBCW[posOpPoint][posFreq] = 0;

    //Sara
    QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
    //Sara

    const double rho = pbem->dlg_rho; //Sara
    const double c_0 = 34000;
    const double u = (m_parameter->originalVelocity)*100.;
    const double c = 100.*(m_parameter->originalChordLength);
    const double I = (m_parameter->TurbulenceIntensity);
    const double Lambda = 100.*(m_parameter->IntegralLengthScale);
    const double r_e = 100.*(m_parameter->distanceObsever);
    const double Mach = m_parameter->originalMach;
    const double beta = sqrt(1-pow(Mach, 2));
    const double D_L = 0.5*getDL();
    const double L = 100.*(m_parameter->wettedLength);
    double Aux = 0.5*(Lambda*L*pow(rho/1000., 2)*pow(c_0, 2)*pow(u, 2)*pow(Mach, 3)*pow(I, 2)*D_L)/(pow(r_e, 2));
    double K = M_PI*CENTRAL_BAND_FREQUENCY[posFreq]*c/u;
    double S = sqrt(pow((2.*M_PI*K/(pow(beta, 2)))+(pow((1+(2.4*K/pow(beta,2))), -1)), -1));
    double LFC = 10.*Mach*pow(S*K/beta, 2);

if(m_parameter->Lowson_type==2){
    c_const=19./6.;
    d_const = 65.95;
 }
if(m_parameter->Lowson_type==1){
    c_const=7./3.;
    d_const = 58.4;
 }

if (m_parameter->Lowson_type!=0){
    Aux1 = 10.*log10(pow(LFC/(1+LFC), 2))+d_const; //Lowson's standard is pow = 1 and const is 58.4
    Aux4 = pow(K, 3)/pow(1+(pow(K, 2)),c_const); //Lowson's standard is 7/3
    Aux5 = 10.*log10(Aux*Aux4);
    aux_SPL_LEdB = 10.*log10(pow(10,(Aux1+Aux5)/10.));

    m_SPL_LEdB[posOpPoint][posFreq] = aux_SPL_LEdB;
    m_SPL_LEdBAW[posOpPoint][posFreq] = aux_SPL_LEdB + AWeighting[posFreq];
    m_SPL_LEdBBW[posOpPoint][posFreq] = aux_SPL_LEdB + BWeighting[posFreq];
    m_SPL_LEdBCW[posOpPoint][posFreq] = aux_SPL_LEdB + CWeighting[posFreq];
}
else{
    m_SPL_LEdB[posOpPoint][posFreq] = 0;
    m_SPL_LEdBAW[posOpPoint][posFreq] = 0;
    m_SPL_LEdBBW[posOpPoint][posFreq] = 0;
    m_SPL_LEdBCW[posOpPoint][posFreq] = 0;
}
}
//end Alexandre MOD

//calculation for 2D noise
void NoiseCalculation::calculate() {
    setupVectors();

    QList<NoiseOpPoint*> noiseOpPoints = m_parameter->prepareNoiseOpPointList();
    for (int posOpPoint = 0; posOpPoint < noiseOpPoints.size(); ++posOpPoint) {
        NoiseOpPoint *nop = noiseOpPoints[posOpPoint];

        m_OASPL[posOpPoint] = 0;
        m_OASPLA[posOpPoint] = 0;
        m_OASPLB[posOpPoint] = 0;
        m_OASPLC[posOpPoint] = 0;
        m_SPLALOG[posOpPoint] = 0;
        m_SPLSLOG[posOpPoint] = 0;
        m_SPLPLOG[posOpPoint] = 0;

//Sara
        m_SPLLEdB[posOpPoint] = 0;
        m_SPLLEdBAW[posOpPoint] = 0;
        m_SPLLEdBBW[posOpPoint] = 0;
        m_SPLLEdBCW[posOpPoint] = 0;
        m_SPLlogLE[posOpPoint] = 0;
//Sara

//        qDebug() << "======================== OpPoint ========================";
//        qDebug() << "Alpha deg: " << nop->AlphaDeg();
//        qDebug() << "Reynolds: " << nop->Reynolds();

        bool dStarOrder = false;

        //When angle is negative D* search must be inverted
        if(nop->getAlphaDegree() < 0){
            dStarOrder = true;
        }


        if (m_parameter->opPointSource == NoiseParameter::OnePolar ||
            m_parameter->opPointSource == NoiseParameter::MultiplePolars)
        {

            //For XFoil model

            m_DStarInterpolatedS = getDStarInterpolated(dStarOrder,nop);
            m_DStarInterpolatedP = getDStarInterpolated(!dStarOrder,nop);

//qDebug() << "D* S original:" << m_DStarInterpolatedS;

            m_DStarFinalS = m_DStarInterpolatedS * m_parameter->originalChordLength * m_parameter->dStarScalingFactor;
            m_DStarFinalP = m_DStarInterpolatedP * m_parameter->originalChordLength * m_parameter->dStarScalingFactor;
        } else if (m_parameter->opPointSource == NoiseParameter::OriginalBpm) {

            m_DStarInterpolatedS = 0;
            m_DStarInterpolatedP = 0;

            //For BPM model

            m_DStarFinalS = getBPMThickness(nop, SuctionSide) * m_parameter->dStarScalingFactor;
            m_DStarFinalP = getBPMThickness(nop, PressureSide) * m_parameter->dStarScalingFactor;

        }

        m_EddyMachNumber = m_parameter->originalMach * m_parameter->eddyConvectionMach;

//        qDebug() << "Linear DStar interpolated Top/Suction: " << m_DStarInterpolatedS;
//        qDebug() << "Final DStar Top/Suction: " << m_DStarFinalS;
//        qDebug() << "Linear DStar interpolated Bottom/Pressure: " << m_DStarInterpolatedP;
//        qDebug() << "Final DStar Bottom/Pressure: " << m_DStarFinalP;
//        qDebug() << "Mach Number: " << m_NoiseParameter->m_OriginalMach;
//        qDebug() << "Velocity: " << m_NoiseParameter->m_OriginalVelocity;
//        qDebug() << "Eddy Mach Number: " << m_EddyMachNumber;

        m_SwAlpha1 = 23.43 * m_parameter->originalMach + 4.651;
        m_SwAlpha = fmin(m_SwAlpha1, SWITCHING_ANGLE2);

        if( nop->getAlphaDegree() <= m_SwAlpha){
            m_AlphaBigSw = false;
        }else{
            m_AlphaBigSw = true;
        }

//        qDebug() << "SwAngle1: " << m_SwAlpha1;
//        qDebug() << "SwAngle calculated: " << m_SwAlpha;


        m_CalcSeparatedFlow = false;
        m_CalcPressureSide = false;
        m_CalcSuctionSide = false;

        if(m_AlphaBigSw && m_parameter->separatedFlow){
//            qDebug() << "Only separated flow source will be calculated";
            m_CalcSeparatedFlow = true;
            m_CalcPressureSide = false;
            m_CalcSuctionSide = false;
        }else{
            if(m_parameter->separatedFlow){
//                qDebug() << "Separated flow source will be calculated";
                m_CalcSeparatedFlow = true;
            }

            if(m_parameter->pressureSide){
//                qDebug() << "Pressure side source will be calculated";
                m_CalcPressureSide = true;
            }

            if(m_parameter->suctionSide){
//                qDebug() << "Suction side source will be calculated";
                m_CalcSuctionSide = true;
            }
        }

        preCalcA1(nop);

        if(m_CalcSeparatedFlow){
            preCalcSPLa(nop);
        }

        //If angle is bigger than the switching Angle
        //or suction side is mandatory
        if(m_CalcSuctionSide){
            preCalcSPLs(nop);
        }

        //If angle is bigger than the switching Angle
        //or pressure side is mandatory
        if(m_CalcPressureSide){
            preCalcSPLp(nop);
        }

        //For each frequency
        //for (unsigned int posFreq = 0; posFreq < 3; ++posFreq) {
        for (unsigned int posFreq = 0; posFreq < FREQUENCY_TABLE_SIZE; ++posFreq) {
//			qDebug() << "==== Band Frequency ====";
//			qDebug() << "Freq: [" << (posFreq+1) << "] " << Noise::CENTRAL_BAND_FREQUENCY[posFreq] ;

            LECalc(posOpPoint, posFreq); //Alexandre MOD

            if (m_CalcSeparatedFlow) {
                calcSPLa(nop->getAlphaDegree(),posOpPoint,posFreq);
            }

            //If angle is bigger than the switching Angle
            //or suction side is mandatory
//			if(m_CalcSuctionSide){
                calcSPLs(posOpPoint,posFreq);
//			}

            //If angle is bigger than the switching Angle
            //or pressure side is mandatory
//            if(m_CalcPressureSide){
                calcSPLp(posOpPoint,posFreq);
//            }

            double splDbConsolidated = 0.0;

            if(m_CalcSeparatedFlow)
                splDbConsolidated += pow(10,(m_SPLadB[posOpPoint][posFreq]/10));

            if(m_CalcPressureSide)
                splDbConsolidated += pow(10,(m_SPLsdB[posOpPoint][posFreq]/10));

            if(m_CalcSuctionSide)
                splDbConsolidated += pow(10,(m_SPLpdB[posOpPoint][posFreq]/10));

//Sara
            if(m_parameter->Lowson_type!=0){
                splDbConsolidated += pow(10,(m_SPL_LEdB[posOpPoint][posFreq]/10));
            }
//Sara

            m_SPLdB[posOpPoint][posFreq] = 10*log10( splDbConsolidated );
            //m_SPLdB[posOpPoint][posFreq] = 10*log10(pow(10,(m_SPLadB[posOpPoint][posFreq]/10))+pow(10,(m_SPLsdB[posOpPoint][posFreq]/10))+pow(10,(m_SPLpdB[posOpPoint][posFreq]/10)));
            m_SPLdBAW[posOpPoint][posFreq] =m_SPLdB[posOpPoint][posFreq] + AWeighting[posFreq];
            m_SPLdBBW[posOpPoint][posFreq] =m_SPLdB[posOpPoint][posFreq] + BWeighting[posFreq];
            m_SPLdBCW[posOpPoint][posFreq] =m_SPLdB[posOpPoint][posFreq] + CWeighting[posFreq];

//            qDebug() << "SPLdb(" << m_SPLdB[posOpPoint][posFreq] << ") " << "SPLdbAW(" << m_SPLdBAW[posOpPoint][posFreq] << ") " << "SPLdbBW(" << m_SPLdBBW[posOpPoint][posFreq] << ") " << "SPLdbCW(" << m_SPLdBCW[posOpPoint][posFreq] << ") ";


            m_OASPL[posOpPoint] += pow(10,(m_SPLdB[posOpPoint][posFreq]/10));
            m_OASPLA[posOpPoint] += pow(10,(m_SPLdBAW[posOpPoint][posFreq]/10));
            m_OASPLB[posOpPoint] += pow(10,(m_SPLdBBW[posOpPoint][posFreq]/10));
            m_OASPLC[posOpPoint] += pow(10,(m_SPLdBCW[posOpPoint][posFreq]/10));

            m_SPLALOG[posOpPoint] += pow(10,(m_SPLadB[posOpPoint][posFreq]/10));
            m_SPLSLOG[posOpPoint] += pow(10,(m_SPLsdB[posOpPoint][posFreq]/10));
            m_SPLPLOG[posOpPoint] += pow(10,(m_SPLpdB[posOpPoint][posFreq]/10));

//Sara
            if (m_parameter->Lowson_type!=0){
            m_SPLlogLE[posOpPoint] += pow(10,(m_SPL_LEdB[posOpPoint][posFreq]/10));
            m_SPLLEdBAW[posOpPoint] += pow(10,(m_SPL_LEdBAW[posOpPoint][posFreq]/10));
            m_SPLLEdBBW[posOpPoint] +=pow(10,(m_SPL_LEdBBW[posOpPoint][posFreq]/10));
            m_SPLLEdBCW[posOpPoint] += pow(10,(m_SPL_LEdBCW[posOpPoint][posFreq]/10));
            }
            else{
            m_SPLlogLE[posOpPoint] = 0;
            m_SPLLEdBAW[posOpPoint] = 0;
            m_SPLLEdBBW[posOpPoint] = 0;
            m_SPLLEdBCW[posOpPoint] = 0;
            }
//Sara
        }

        m_OASPL[posOpPoint] = 10*log10(m_OASPL[posOpPoint]);
        m_OASPLA[posOpPoint] = 10*log10(m_OASPLA[posOpPoint]);
        m_OASPLB[posOpPoint] = 10*log10(m_OASPLB[posOpPoint]);
        m_OASPLC[posOpPoint] = 10*log10(m_OASPLC[posOpPoint]);

        m_SPLALOG[posOpPoint] = 10*log10(m_SPLALOG[posOpPoint]);
        m_SPLSLOG[posOpPoint] = 10*log10(m_SPLSLOG[posOpPoint]);
        m_SPLPLOG[posOpPoint] = 10*log10(m_SPLPLOG[posOpPoint]);

//Sara
m_SPLlogLE[posOpPoint] = 10*log10(m_SPLlogLE[posOpPoint]);
m_SPLLEdBAW[posOpPoint] = 10*log10(m_SPLLEdBAW[posOpPoint]);
m_SPLLEdBBW[posOpPoint] = 10*log10(m_SPLLEdBBW[posOpPoint]);
m_SPLLEdBCW[posOpPoint] = 10*log10(m_SPLLEdBCW[posOpPoint]);
//Sara
    }
    qDeleteAll(noiseOpPoints);

//	// transform negative values to zero  // NM this apparently wasn't intended...
//	for (TwoDVector *vector2d : {&m_SPLadB, &m_SPLsdB, &m_SPLpdB}) {
//		for (int i = 0; i < vector2d->size(); ++i) {
//			for (int j = 0; j < (*vector2d)[i].size(); ++j) {
//				(*vector2d)[i][j] = std::max(0.0, (*vector2d)[i][j]);
//			}
//		}
//	}

ProgressBar(1);//Sara
}

//setup 2D vectors
void NoiseCalculation::setupVectors() {
    // Clear all results is case of simulation edition
    m_SPLdB.clear();
    m_SPLdBAW.clear();
    m_SPLdBBW.clear();
    m_SPLdBCW.clear();
    m_SPLpdB.clear();
    m_SPLpdBAW.clear();
    m_SPLpdBBW.clear();
    m_SPLpdBCW.clear();
    m_SPLsdB.clear();
    m_SPLsdBAW.clear();
    m_SPLsdBBW.clear();
    m_SPLsdBCW.clear();
    m_SPLadB.clear();
    m_SPLadBAW.clear();
    m_SPLadBBW.clear();
    m_SPLadBCW.clear();
    //Alexandre MOD
    m_SPL_LEdB.clear();
    m_SPL_LEdBAW.clear();
    m_SPL_LEdBBW.clear();
    m_SPL_LEdBCW.clear();
//Sara
    m_SPLLEdB.clear();
    m_SPLLEdBAW.clear();
    m_SPLLEdBBW.clear();
    m_SPLLEdBCW.clear();
    m_SPLlogLE.clear();
//Sara

    m_OASPL.clear();
    m_OASPLA.clear();
    m_OASPLB.clear();
    m_OASPLC.clear();
    m_SPLALOG.clear();
    m_SPLSLOG.clear();
    m_SPLPLOG.clear();

    // Resize vectors acording to OpPoints total
    unsigned int sizea = 0;//Sara size
    if (m_parameter->opPointSource == NoiseParameter::OnePolar ||
        m_parameter->opPointSource == NoiseParameter::MultiplePolars)
    {
        sizea = m_parameter->analyzedOpPoints.size();//Sara size
    } else if (m_parameter->opPointSource == NoiseParameter::OriginalBpm) {
        sizea = 1;//Sara size
    }

    //Sara
  QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
  unsigned int number_of_segments = pBEM->dlg_elements;

  unsigned int size;

  if (number_of_segments>sizea){size = number_of_segments;} else {size = sizea;}
    //Sara

    m_SPLdB.resize   (size);
    m_SPLdBAW.resize (size);
    m_SPLdBBW.resize (size);
    m_SPLdBCW.resize (size);
    m_SPLpdB.resize  (size);
    m_SPLpdBAW.resize(size);
    m_SPLpdBBW.resize(size);
    m_SPLpdBCW.resize(size);
    m_SPLsdB.resize  (size);
    m_SPLsdBAW.resize(size);
    m_SPLsdBBW.resize(size);
    m_SPLsdBCW.resize(size);
    m_SPLadB.resize  (size);
    m_SPLadBAW.resize(size);
    m_SPLadBBW.resize(size);
    m_SPLadBCW.resize(size);
    //Alexandre MOD
    m_SPL_LEdB.resize(size);
    m_SPL_LEdBAW.resize(size);
    m_SPL_LEdBBW.resize(size);
    m_SPL_LEdBCW.resize(size);
    m_SPLLEdB.resize(size);
//Sara
    m_SPLLEdBAW.resize(size);
    m_SPLLEdBBW.resize(size);
    m_SPLLEdBCW.resize(size);
    m_SPLlogLE.resize(size);
//Sara
    m_OASPL.resize  (size);
    m_OASPLA.resize (size);
    m_OASPLB.resize (size);
    m_OASPLC.resize (size);
    m_SPLALOG.resize(size);
    m_SPLSLOG.resize(size);
    m_SPLPLOG.resize(size);

    // Resize for each OpPoint the frequency table
    for (unsigned int i = 0; i < size; ++i) {
        m_SPLdB[i].resize   (FREQUENCY_TABLE_SIZE);
        m_SPLdBAW[i].resize (FREQUENCY_TABLE_SIZE);
        m_SPLdBBW[i].resize (FREQUENCY_TABLE_SIZE);
        m_SPLdBCW[i].resize (FREQUENCY_TABLE_SIZE);
        m_SPLpdB[i].resize  (FREQUENCY_TABLE_SIZE);
        m_SPLpdBAW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdBBW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdBCW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB[i].resize  (FREQUENCY_TABLE_SIZE);
        m_SPLsdBAW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdBBW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdBCW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPLadB[i].resize  (FREQUENCY_TABLE_SIZE);
        m_SPLadBAW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPLadBBW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPLadBCW[i].resize(FREQUENCY_TABLE_SIZE);
        //Alexandre MOD
        m_SPL_LEdB[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW[i].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW[i].resize(FREQUENCY_TABLE_SIZE);
    }

//Sara
    if (sizea<size){
        for (unsigned int i=sizea;i<size;++i){
            for (unsigned int w = 0; w < FREQUENCY_TABLE_SIZE; ++w){
                m_SPLdB[i][w]=0;
                m_SPLdBAW[i][w]=0;
                m_SPLdBBW[i][w]=0;
                m_SPLdBCW[i][w]=0;
                m_SPLpdB[i][w]=0;
                m_SPLpdBAW[i][w]=0;
                m_SPLpdBBW[i][w]=0;
                m_SPLpdBCW[i][w]=0;
                m_SPLsdB[i][w]=0;
                m_SPLsdBAW[i][w]=0;
                m_SPLsdBBW[i][w]=0;
                m_SPLsdBCW[i][w]=0;
                m_SPLadB[i][w]=0;
                m_SPLadBAW[i][w]=0;
                m_SPLadBBW[i][w]=0;
                m_SPLadBCW[i][w]=0;
                m_SPL_LEdB[i][w]=0;
                m_SPL_LEdBAW[i][w]=0;
                m_SPL_LEdBBW[i][w]=0;
                m_SPL_LEdBCW[i][w]=0;
            }
        }
    }
//Sara
}

//setup qs3D vectors
void NoiseCalculation::setupVectorsqs3d() {
    m_SPLdB3d.clear();
    m_SPLdBAW3d.clear();
    m_SPLdBBW3d.clear();
    m_SPLdBCW3d.clear();
    m_SPLpdB3d.clear();
    m_SPLpdBAW3d.clear();
    m_SPLpdBBW3d.clear();
    m_SPLpdBCW3d.clear();
    m_SPLsdB3d.clear();
    m_SPLsdBAW3d.clear();
    m_SPLsdBBW3d.clear();
    m_SPLsdBCW3d.clear();
    m_SPLadB3d.clear();
    m_SPLadBAW3d.clear();
    m_SPLadBBW3d.clear();
    m_SPLadBCW3d.clear();
    m_SPL_LEdB3d.clear();
    m_SPL_LEdBAW3d.clear();
    m_SPL_LEdBBW3d.clear();
    m_SPL_LEdBCW3d.clear();
    m_SPLLEdB3d.clear();
    m_SPLLEdBAW3d.clear();
    m_SPLLEdBBW3d.clear();
    m_SPLLEdBCW3d.clear();
    m_SPLlogLE3d.clear();
    m_OASPL3d.clear();
    m_OASPLA3d.clear();
    m_OASPLB3d.clear();
    m_OASPLC3d.clear();
    m_SPLALOG3d.clear();
    m_SPLSLOG3d.clear();
    m_SPLPLOG3d.clear();

    m_SPLdB3d_final.clear();
    m_SPLdBAW3d_final.clear();
    m_SPLdBBW3d_final.clear();
    m_SPLdBCW3d_final.clear();
    m_SPLpdB3d_final.clear();
    m_SPLsdB3d_final.clear();
    m_SPLadB3d_final.clear();
    m_SPL_LEdB3d_final.clear();
    m_SPL_LEdBAW3d_final.clear();
    m_SPL_LEdBBW3d_final.clear();
    m_SPL_LEdBCW3d_final.clear();

    m_SPLdB3d_rotor.clear();
    m_SPLdBAW3d_rotor.clear();
    m_SPLdBBW3d_rotor.clear();
    m_SPLdBCW3d_rotor.clear();
    m_SPLpdB3d_rotor.clear();
    m_SPLpdBAW3d_rotor.clear();
    m_SPLpdBBW3d_rotor.clear();
    m_SPLpdBCW3d_rotor.clear();
    m_SPLsdB3d_rotor.clear();
    m_SPLsdBAW3d_rotor.clear();
    m_SPLsdBBW3d_rotor.clear();
    m_SPLsdBCW3d_rotor.clear();
    m_SPLadB3d_rotor.clear();
    m_SPLadBAW3d_rotor.clear();
    m_SPLadBBW3d_rotor.clear();
    m_SPLadBCW3d_rotor.clear();
    m_SPL_LEdB3d_rotor.clear();
    m_SPL_LEdBAW3d_rotor.clear();
    m_SPL_LEdBBW3d_rotor.clear();
    m_SPL_LEdBCW3d_rotor.clear();
    m_SPLLEdB3d_rotor.clear();
    m_SPLLEdBAW3d_rotor.clear();
    m_SPLLEdBBW3d_rotor.clear();
    m_SPLLEdBCW3d_rotor.clear();
    m_SPLlogLE3d_rotor.clear();
    m_OASPL3d_rotor.clear();
    m_OASPLA3d_rotor.clear();
    m_OASPLB3d_rotor.clear();
    m_OASPLC3d_rotor.clear();
    m_SPLALOG3d_rotor.clear();
    m_SPLSLOG3d_rotor.clear();
    m_SPLPLOG3d_rotor.clear();

    m_SPLdB3d_final_rotor.clear();
    m_SPLdBAW3d_final_rotor.clear();
    m_SPLdBBW3d_final_rotor.clear();
    m_SPLdBCW3d_final_rotor.clear();
    m_SPLpdB3d_final_rotor.clear();
    m_SPLsdB3d_final_rotor.clear();
    m_SPLadB3d_final_rotor.clear();
    m_SPL_LEdB3d_final_rotor.clear();
    m_SPL_LEdBAW3d_final_rotor.clear();
    m_SPL_LEdBBW3d_final_rotor.clear();
    m_SPL_LEdBCW3d_final_rotor.clear();

    m_SPLdB3d_rotor_loops.clear();
    m_SPLdBAW3d_rotor_loops.clear();
    m_SPLdBBW3d_rotor_loops.clear();
    m_SPLdBCW3d_rotor_loops.clear();
    m_SPLpdB3d_rotor_loops.clear();
    m_SPLpdBAW3d_rotor_loops.clear();
    m_SPLpdBBW3d_rotor_loops.clear();
    m_SPLpdBCW3d_rotor_loops.clear();
    m_SPLsdB3d_rotor_loops.clear();
    m_SPLsdBAW3d_rotor_loops.clear();
    m_SPLsdBBW3d_rotor_loops.clear();
    m_SPLsdBCW3d_rotor_loops.clear();
    m_SPLadB3d_rotor_loops.clear();
    m_SPLadBAW3d_rotor_loops.clear();
    m_SPLadBBW3d_rotor_loops.clear();
    m_SPLadBCW3d_rotor_loops.clear();
    m_SPL_LEdB3d_rotor_loops.clear();
    m_SPL_LEdBAW3d_rotor_loops.clear();
    m_SPL_LEdBBW3d_rotor_loops.clear();
    m_SPL_LEdBCW3d_rotor_loops.clear();
    m_SPLLEdB3d_rotor_loops.clear();
    m_SPLLEdBAW3d_rotor_loops.clear();
    m_SPLLEdBBW3d_rotor_loops.clear();
    m_SPLLEdBCW3d_rotor_loops.clear();
    m_SPLlogLE3d_rotor_loops.clear();
    m_OASPL3d_rotor_loops.clear();
    m_OASPLA3d_rotor_loops.clear();
    m_OASPLB3d_rotor_loops.clear();
    m_OASPLC3d_rotor_loops.clear();
    m_SPLALOG3d_rotor_loops.clear();
    m_SPLSLOG3d_rotor_loops.clear();
    m_SPLPLOG3d_rotor_loops.clear();

    m_SPLdB3d_final_rotor_loops.clear();
    m_SPLdBAW3d_final_rotor_loops.clear();
    m_SPLdBBW3d_final_rotor_loops.clear();
    m_SPLdBCW3d_final_rotor_loops.clear();
    m_SPLpdB3d_final_rotor_loops.clear();
    m_SPLsdB3d_final_rotor_loops.clear();
    m_SPLadB3d_final_rotor_loops.clear();
    m_SPL_LEdB3d_final_rotor_loops.clear();
    m_SPL_LEdBAW3d_final_rotor_loops.clear();
    m_SPL_LEdBBW3d_final_rotor_loops.clear();
    m_SPL_LEdBCW3d_final_rotor_loops.clear();

    m_SPLdB3d_4d.clear();
    m_SPLdBAW3d_4d.clear();
    m_SPLdBBW3d_4d.clear();
    m_SPLdBCW3d_4d.clear();
    m_SPLpdB3d_4d.clear();
    m_SPLpdBAW3d_4d.clear();
    m_SPLpdBBW3d_4d.clear();
    m_SPLpdBCW3d_4d.clear();
    m_SPLsdB3d_4d.clear();
    m_SPLsdBAW3d_4d.clear();
    m_SPLsdBBW3d_4d.clear();
    m_SPLsdBCW3d_4d.clear();
    m_SPLadB3d_4d.clear();
    m_SPLadBAW3d_4d.clear();
    m_SPLadBBW3d_4d.clear();
    m_SPLadBCW3d_4d.clear();
    m_SPL_LEdB3d_4d.clear();
    m_SPL_LEdBAW3d_4d.clear();
    m_SPL_LEdBBW3d_4d.clear();
    m_SPL_LEdBCW3d_4d.clear();

    m_SPLdB3d_4d_blade.clear();
    m_SPLdBAW3d_4d_blade.clear();
    m_SPLdBBW3d_4d_blade.clear();
    m_SPLdBCW3d_4d_blade.clear();
    m_SPLpdB3d_4d_blade.clear();
    m_SPLpdBAW3d_4d_blade.clear();
    m_SPLpdBBW3d_4d_blade.clear();
    m_SPLpdBCW3d_4d_blade.clear();
    m_SPLsdB3d_4d_blade.clear();
    m_SPLsdBAW3d_4d_blade.clear();
    m_SPLsdBBW3d_4d_blade.clear();
    m_SPLsdBCW3d_4d_blade.clear();
    m_SPLadB3d_4d_blade.clear();
    m_SPLadBAW3d_4d_blade.clear();
    m_SPLadBBW3d_4d_blade.clear();
    m_SPLadBCW3d_4d_blade.clear();
    m_SPL_LEdB3d_4d_blade.clear();
    m_SPL_LEdBAW3d_4d_blade.clear();
    m_SPL_LEdBBW3d_4d_blade.clear();
    m_SPL_LEdBCW3d_4d_blade.clear();

    // Resize vectors acording to OpPoints total
    unsigned int sizea = 0;
    if (m_parameter->opPointSource == NoiseParameter::OnePolar ||
        m_parameter->opPointSource == NoiseParameter::MultiplePolars)
    {
        sizea = m_parameter->analyzedOpPoints.size();
     } else if (m_parameter->opPointSource == NoiseParameter::OriginalBpm) {
        sizea = 1;
    }

QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
unsigned int number_of_segments = pBEM->dlg_elements;

unsigned int size;

if (number_of_segments>sizea){size = number_of_segments;} else {size = sizea;}
    m_SPLadB3d.resize(size);
    m_SPLsdB3d.resize(size);
    m_SPLpdB3d.resize(size);
    m_SPLdB3d.resize(size);
    m_SPLdBAW3d.resize(size);
    m_SPLdBBW3d.resize(size);
    m_SPLdBCW3d.resize(size);
    m_SPL_LEdB3d.resize(size);
    m_SPL_LEdBAW3d.resize(size);
    m_SPL_LEdBBW3d.resize(size);
    m_SPL_LEdBCW3d.resize(size);
    m_SPLdB3d_final.resize(size);
    m_SPLdBAW3d_final.resize(size);
    m_SPLdBBW3d_final.resize(size);
    m_SPLdBCW3d_final.resize(size);
    m_SPLpdB3d_final.resize(size);
    m_SPLsdB3d_final.resize(size);
    m_SPLadB3d_final.resize(size);
    m_SPL_LEdB3d_final.resize(size);
    m_SPL_LEdBAW3d_final.resize(size);
    m_SPL_LEdBBW3d_final.resize(size);
    m_SPL_LEdBCW3d_final.resize(size);
    m_SPLLEdB3d.resize(size);
    m_SPLLEdBAW3d.resize(size);
    m_SPLLEdBBW3d.resize(size);
    m_SPLLEdBCW3d.resize(size);
    m_SPLlogLE3d.resize(size);
    m_OASPL3d.resize(size);
    m_OASPLA3d.resize(size);
    m_OASPLB3d.resize(size);
    m_OASPLC3d.resize(size);
    m_SPLALOG3d.resize(size);
    m_SPLSLOG3d.resize(size);
    m_SPLPLOG3d.resize(size);

    m_SPLadB3d_rotor.resize(size);
    m_SPLsdB3d_rotor.resize(size);
    m_SPLpdB3d_rotor.resize(size);
    m_SPLdB3d_rotor.resize(size);
    m_SPLdBAW3d_rotor.resize(size);
    m_SPLdBBW3d_rotor.resize(size);
    m_SPLdBCW3d_rotor.resize(size);
    m_SPL_LEdB3d_rotor.resize(size);
    m_SPL_LEdBAW3d_rotor.resize(size);
    m_SPL_LEdBBW3d_rotor.resize(size);
    m_SPL_LEdBCW3d_rotor.resize(size);
    m_SPLdB3d_final_rotor.resize(size);
    m_SPLdBAW3d_final_rotor.resize(size);
    m_SPLdBBW3d_final_rotor.resize(size);
    m_SPLdBCW3d_final_rotor.resize(size);
    m_SPLpdB3d_final_rotor.resize(size);
    m_SPLsdB3d_final_rotor.resize(size);
    m_SPLadB3d_final_rotor.resize(size);
    m_SPL_LEdB3d_final_rotor.resize(size);
    m_SPL_LEdBAW3d_final_rotor.resize(size);
    m_SPL_LEdBBW3d_final_rotor.resize(size);
    m_SPL_LEdBCW3d_final_rotor.resize(size);
    m_SPLLEdB3d_rotor.resize(size);
    m_SPLLEdBAW3d_rotor.resize(size);
    m_SPLLEdBBW3d_rotor.resize(size);
    m_SPLLEdBCW3d_rotor.resize(size);
    m_SPLlogLE3d_rotor.resize(size);
    m_OASPL3d_rotor.resize(size);
    m_OASPLA3d_rotor.resize(size);
    m_OASPLB3d_rotor.resize(size);
    m_OASPLC3d_rotor.resize(size);
    m_SPLALOG3d_rotor.resize(size);
    m_SPLSLOG3d_rotor.resize(size);
    m_SPLPLOG3d_rotor.resize(size);

    m_SPLadB3d_rotor_loops.resize(size);
    m_SPLsdB3d_rotor_loops.resize(size);
    m_SPLpdB3d_rotor_loops.resize(size);
    m_SPLdB3d_rotor_loops.resize(size);
    m_SPLdBAW3d_rotor_loops.resize(size);
    m_SPLdBBW3d_rotor_loops.resize(size);
    m_SPLdBCW3d_rotor_loops.resize(size);
    m_SPL_LEdB3d_rotor_loops.resize(size);
    m_SPL_LEdBAW3d_rotor_loops.resize(size);
    m_SPL_LEdBBW3d_rotor_loops.resize(size);
    m_SPL_LEdBCW3d_rotor_loops.resize(size);
    m_SPLdB3d_final_rotor_loops.resize(size);
    m_SPLdBAW3d_final_rotor_loops.resize(size);
    m_SPLdBBW3d_final_rotor_loops.resize(size);
    m_SPLdBCW3d_final_rotor_loops.resize(size);
    m_SPLpdB3d_final_rotor_loops.resize(size);
    m_SPLsdB3d_final_rotor_loops.resize(size);
    m_SPLadB3d_final_rotor_loops.resize(size);
    m_SPL_LEdB3d_final_rotor_loops.resize(size);
    m_SPL_LEdBAW3d_final_rotor_loops.resize(size);
    m_SPL_LEdBBW3d_final_rotor_loops.resize(size);
    m_SPL_LEdBCW3d_final_rotor_loops.resize(size);
    m_SPLLEdB3d_rotor_loops.resize(size);
    m_SPLLEdBAW3d_rotor_loops.resize(size);
    m_SPLLEdBBW3d_rotor_loops.resize(size);
    m_SPLLEdBCW3d_rotor_loops.resize(size);
    m_SPLlogLE3d_rotor_loops.resize(size);
    m_OASPL3d_rotor_loops.resize(size);
    m_OASPLA3d_rotor_loops.resize(size);
    m_OASPLB3d_rotor_loops.resize(size);
    m_OASPLC3d_rotor_loops.resize(size);
    m_SPLALOG3d_rotor_loops.resize(size);
    m_SPLSLOG3d_rotor_loops.resize(size);
    m_SPLPLOG3d_rotor_loops.resize(size);

    m_SPLadB3d_4d.resize(size);
    m_SPLsdB3d_4d.resize(size);
    m_SPLpdB3d_4d.resize(size);
    m_SPLdB3d_4d.resize(size);
    m_SPLdBAW3d_4d.resize(size);
    m_SPLdBBW3d_4d.resize(size);
    m_SPLdBCW3d_4d.resize(size);
    m_SPL_LEdB3d_4d.resize(size);
    m_SPL_LEdBAW3d_4d.resize(size);
    m_SPL_LEdBBW3d_4d.resize(size);
    m_SPL_LEdBCW3d_4d.resize(size);

    m_SPLadB3d_4d_blade.resize(size);
    m_SPLsdB3d_4d_blade.resize(size);
    m_SPLpdB3d_4d_blade.resize(size);
    m_SPLdB3d_4d_blade.resize(size);
    m_SPLdBAW3d_4d_blade.resize(size);
    m_SPLdBBW3d_4d_blade.resize(size);
    m_SPLdBCW3d_4d_blade.resize(size);
    m_SPL_LEdB3d_4d_blade.resize(size);
    m_SPL_LEdBAW3d_4d_blade.resize(size);
    m_SPL_LEdBBW3d_4d_blade.resize(size);
    m_SPL_LEdBCW3d_4d_blade.resize(size);

    for (unsigned int w = 0; w < size; ++w){
        m_SPLadB3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdB3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdB3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBAW3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBBW3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBCW3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdB3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW3d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdB3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBAW3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBBW3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBCW3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdB3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLadB3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdB3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW3d_final[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW3d_final[w].resize(FREQUENCY_TABLE_SIZE);

        m_SPLadB3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdB3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdB3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBAW3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBBW3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBCW3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdB3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW3d_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdB3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBAW3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBBW3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBCW3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdB3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLadB3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdB3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW3d_final_rotor[w].resize(FREQUENCY_TABLE_SIZE);

        m_SPLadB3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdB3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdB3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBAW3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBBW3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBCW3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdB3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW3d_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdB3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBAW3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBBW3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBCW3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdB3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLadB3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdB3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW3d_final_rotor_loops[w].resize(FREQUENCY_TABLE_SIZE);

        m_SPLadB3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdB3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdB3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBAW3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBBW3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBCW3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdB3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW3d_4d[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW3d_4d[w].resize(FREQUENCY_TABLE_SIZE);

        m_SPLadB3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLsdB3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLpdB3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdB3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBAW3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBBW3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPLdBCW3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdB3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBAW3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBBW3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
        m_SPL_LEdBCW3d_4d_blade[w].resize(FREQUENCY_TABLE_SIZE);
    }

QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
unsigned int blades_num = pbem->m_pBData->blades; //number of blades

for (unsigned int w = 0; w < size; ++w){
for (unsigned int s = 0; s < FREQUENCY_TABLE_SIZE; ++s){
        m_SPLadB3d_4d[w][s].resize(blades_num);
        m_SPLsdB3d_4d[w][s].resize(blades_num);
        m_SPLpdB3d_4d[w][s].resize(blades_num);
        m_SPLdB3d_4d[w][s].resize(blades_num);
        m_SPLdBAW3d_4d[w][s].resize(blades_num);
        m_SPLdBBW3d_4d[w][s].resize(blades_num);
        m_SPLdBCW3d_4d[w][s].resize(blades_num);
        m_SPL_LEdB3d_4d[w][s].resize(blades_num);
        m_SPL_LEdBAW3d_4d[w][s].resize(blades_num);
        m_SPL_LEdBBW3d_4d[w][s].resize(blades_num);
        m_SPL_LEdBCW3d_4d[w][s].resize(blades_num);

        m_SPLadB3d_4d_blade[w][s].resize(blades_num);
        m_SPLsdB3d_4d_blade[w][s].resize(blades_num);
        m_SPLpdB3d_4d_blade[w][s].resize(blades_num);
        m_SPLdB3d_4d_blade[w][s].resize(blades_num);
        m_SPLdBAW3d_4d_blade[w][s].resize(blades_num);
        m_SPLdBBW3d_4d_blade[w][s].resize(blades_num);
        m_SPLdBCW3d_4d_blade[w][s].resize(blades_num);
        m_SPL_LEdB3d_4d_blade[w][s].resize(blades_num);
        m_SPL_LEdBAW3d_4d_blade[w][s].resize(blades_num);
        m_SPL_LEdBBW3d_4d_blade[w][s].resize(blades_num);
        m_SPL_LEdBCW3d_4d_blade[w][s].resize(blades_num);
}}

int angles_num;
int anglesteps;
int number_of_rotations;

if (m_parameter->rotation_type==0){
//    angle based
number_of_rotations = m_parameter->number_loops;
anglesteps=m_parameter->anglesteps;
}else{
//    time based
number_of_rotations = m_parameter->time/(60./m_parameter->rot_speed);
anglesteps=m_parameter->timesteps*60.*360./(m_parameter->rot_speed*1000.);
}

angles_num = 360./anglesteps*number_of_rotations;

for (unsigned int w = 0; w < size; ++w){
for (unsigned int s = 0; s < FREQUENCY_TABLE_SIZE; ++s){
for (unsigned int x = 0; x < blades_num; ++x){
    m_SPLadB3d_4d[w][s][x].resize(angles_num);
    m_SPLsdB3d_4d[w][s][x].resize(angles_num);
    m_SPLpdB3d_4d[w][s][x].resize(angles_num);
    m_SPLdB3d_4d[w][s][x].resize(angles_num);
    m_SPLdBAW3d_4d[w][s][x].resize(angles_num);
    m_SPLdBBW3d_4d[w][s][x].resize(angles_num);
    m_SPLdBCW3d_4d[w][s][x].resize(angles_num);
    m_SPL_LEdB3d_4d[w][s][x].resize(angles_num);
    m_SPL_LEdBAW3d_4d[w][s][x].resize(angles_num);
    m_SPL_LEdBBW3d_4d[w][s][x].resize(angles_num);
    m_SPL_LEdBCW3d_4d[w][s][x].resize(angles_num);

    m_SPLadB3d_4d_blade[w][s][x].resize(angles_num);
    m_SPLsdB3d_4d_blade[w][s][x].resize(angles_num);
    m_SPLpdB3d_4d_blade[w][s][x].resize(angles_num);
    m_SPLdB3d_4d_blade[w][s][x].resize(angles_num);
    m_SPLdBAW3d_4d_blade[w][s][x].resize(angles_num);
    m_SPLdBBW3d_4d_blade[w][s][x].resize(angles_num);
    m_SPLdBCW3d_4d_blade[w][s][x].resize(angles_num);
    m_SPL_LEdB3d_4d_blade[w][s][x].resize(angles_num);
    m_SPL_LEdBAW3d_4d_blade[w][s][x].resize(angles_num);
    m_SPL_LEdBBW3d_4d_blade[w][s][x].resize(angles_num);
    m_SPL_LEdBCW3d_4d_blade[w][s][x].resize(angles_num);
}}}

if (sizea<size){
    for (unsigned int i=sizea;i<size;++i){
        for (unsigned int w = 0; w < FREQUENCY_TABLE_SIZE; ++w){
            m_SPLadB3d[i][w]=0;
            m_SPLsdB3d[i][w]=0;
            m_SPLpdB3d[i][w]=0;
            m_SPLdB3d[i][w]=0;
            m_SPLdBAW3d[i][w]=0;
            m_SPLdBBW3d[i][w]=0;
            m_SPLdBCW3d[i][w]=0;
            m_SPL_LEdB3d[i][w]=0;
            m_SPL_LEdBAW3d[i][w]=0;
            m_SPL_LEdBBW3d[i][w]=0;
            m_SPL_LEdBCW3d[i][w]=0;
            m_SPLdB3d_final[i][w]=0;
            m_SPLdBAW3d_final[i][w]=0;
            m_SPLdBBW3d_final[i][w]=0;
            m_SPLdBCW3d_final[i][w]=0;
            m_SPLpdB3d_final[i][w]=0;
            m_SPLsdB3d_final[i][w]=0;
            m_SPLadB3d_final[i][w]=0;
            m_SPL_LEdB3d_final[i][w]=0;
            m_SPL_LEdBAW3d_final[i][w]=0;
            m_SPL_LEdBBW3d_final[i][w]=0;
            m_SPL_LEdBCW3d_final[i][w]=0;

            m_SPLadB3d_rotor[i][w]=0;
            m_SPLsdB3d_rotor[i][w]=0;
            m_SPLpdB3d_rotor[i][w]=0;
            m_SPLdB3d_rotor[i][w]=0;
            m_SPLdBAW3d_rotor[i][w]=0;
            m_SPLdBBW3d_rotor[i][w]=0;
            m_SPLdBCW3d_rotor[i][w]=0;
            m_SPL_LEdB3d_rotor[i][w]=0;
            m_SPL_LEdBAW3d_rotor[i][w]=0;
            m_SPL_LEdBBW3d_rotor[i][w]=0;
            m_SPL_LEdBCW3d_rotor[i][w]=0;
            m_SPLdB3d_final_rotor[i][w]=0;
            m_SPLdBAW3d_final_rotor[i][w]=0;
            m_SPLdBBW3d_final_rotor[i][w]=0;
            m_SPLdBCW3d_final_rotor[i][w]=0;
            m_SPLpdB3d_final_rotor[i][w]=0;
            m_SPLsdB3d_final_rotor[i][w]=0;
            m_SPLadB3d_final_rotor[i][w]=0;
            m_SPL_LEdB3d_final_rotor[i][w]=0;
            m_SPL_LEdBAW3d_final_rotor[i][w]=0;
            m_SPL_LEdBBW3d_final_rotor[i][w]=0;
            m_SPL_LEdBCW3d_final_rotor[i][w]=0;

            m_SPLadB3d_rotor_loops[i][w]=0;
            m_SPLsdB3d_rotor_loops[i][w]=0;
            m_SPLpdB3d_rotor_loops[i][w]=0;
            m_SPLdB3d_rotor_loops[i][w]=0;
            m_SPLdBAW3d_rotor_loops[i][w]=0;
            m_SPLdBBW3d_rotor_loops[i][w]=0;
            m_SPLdBCW3d_rotor_loops[i][w]=0;
            m_SPL_LEdB3d_rotor_loops[i][w]=0;
            m_SPL_LEdBAW3d_rotor_loops[i][w]=0;
            m_SPL_LEdBBW3d_rotor_loops[i][w]=0;
            m_SPL_LEdBCW3d_rotor_loops[i][w]=0;
            m_SPLdB3d_final_rotor_loops[i][w]=0;
            m_SPLdBAW3d_final_rotor_loops[i][w]=0;
            m_SPLdBBW3d_final_rotor_loops[i][w]=0;
            m_SPLdBCW3d_final_rotor_loops[i][w]=0;
            m_SPLpdB3d_final_rotor_loops[i][w]=0;
            m_SPLsdB3d_final_rotor_loops[i][w]=0;
            m_SPLadB3d_final_rotor_loops[i][w]=0;
            m_SPL_LEdB3d_final_rotor_loops[i][w]=0;
            m_SPL_LEdBAW3d_final_rotor_loops[i][w]=0;
            m_SPL_LEdBBW3d_final_rotor_loops[i][w]=0;
            m_SPL_LEdBCW3d_final_rotor_loops[i][w]=0;



            for (unsigned int j = 0; j < blades_num; ++j){
                for (int k = 0; k < angles_num; ++k){
                    m_SPLadB3d_4d[i][w][j][k]=0;
                    m_SPLsdB3d_4d[i][w][j][k]=0;
                    m_SPLpdB3d_4d[i][w][j][k]=0;
                    m_SPLdB3d_4d[i][w][j][k]=0;
                    m_SPLdBAW3d_4d[i][w][j][k]=0;
                    m_SPLdBBW3d_4d[i][w][j][k]=0;
                    m_SPLdBCW3d_4d[i][w][j][k]=0;
                    m_SPL_LEdB3d_4d[i][w][j][k]=0;
                    m_SPL_LEdBAW3d_4d[i][w][j][k]=0;
                    m_SPL_LEdBBW3d_4d[i][w][j][k]=0;
                    m_SPL_LEdBCW3d_4d[i][w][j][k]=0;

                    m_SPLadB3d_4d_blade[i][w][j][k]=0;
                    m_SPLsdB3d_4d_blade[i][w][j][k]=0;
                    m_SPLpdB3d_4d_blade[i][w][j][k]=0;
                    m_SPLdB3d_4d_blade[i][w][j][k]=0;
                    m_SPLdBAW3d_4d_blade[i][w][j][k]=0;
                    m_SPLdBBW3d_4d_blade[i][w][j][k]=0;
                    m_SPLdBCW3d_4d_blade[i][w][j][k]=0;
                    m_SPL_LEdB3d_4d_blade[i][w][j][k]=0;
                    m_SPL_LEdBAW3d_4d_blade[i][w][j][k]=0;
                    m_SPL_LEdBBW3d_4d_blade[i][w][j][k]=0;
                    m_SPL_LEdBCW3d_4d_blade[i][w][j][k]=0;
            }}
        }
    }
}}

//Sara
double NoiseCalculation::calcXRS(double a, double XB, double YB){
double XRS = XB*cos(qDegreesToRadians(a))+YB*sin(qDegreesToRadians(a));
return XRS;
}

double NoiseCalculation::calcYRS(double a, double XB, double YB){
double YRS = -XB*sin(qDegreesToRadians(a))+YB*cos(qDegreesToRadians(a));
return YRS;
}

double NoiseCalculation::calcZRS(double ZB, double r_0, double r_1){
double ZRS = ZB-(r_1-r_0)/2.;
return ZRS;
}

double NoiseCalculation::calcInt_a(double YRS, double c_0, double c_1){
double calc_int_a = YRS-0.75*(c_1-c_0)/2.;
return calc_int_a;
}

double NoiseCalculation::calcXRT(double XRS){
double XRT = XRS;
return XRT;
}

double NoiseCalculation::calcYRT(double b, double calc_int_a, double ZRS){
double YRT = cos(qDegreesToRadians(b))*calc_int_a+sin(qDegreesToRadians(b))*ZRS;
return YRT;
}

double NoiseCalculation::calcZRT(double b, double calc_int_a, double ZRS){
double ZRT = -sin(qDegreesToRadians(b))*calc_int_a+cos(qDegreesToRadians(b))*ZRS;
return ZRT;
}

double NoiseCalculation::calcR_e(double XRT, double YRT, double ZRT){
double r_e=sqrt(pow(XRT,2)+pow(YRT,2)+pow(ZRT,2));
return r_e;
}

double NoiseCalculation::calcTheta_e(double YRT, double ZRT){
double theta_e=qRadiansToDegrees(qAtan(ZRT/YRT));
return theta_e;
}

double NoiseCalculation::calcPhi_e(double XRT, double ZRT){
double phi_e=qRadiansToDegrees(qAtan(XRT/ZRT));
return phi_e;
}

double NoiseCalculation::calcFirstTerm(double Mach, double L, double D, double D_starred, double dist_obs){
double first_term=10.*log10(pow(Mach,5)*L*D*D_starred/pow(dist_obs,2));
return first_term;
}

double NoiseCalculation::calcDh(double Mach, double theta_e, double phi_e,double EddyMach){
double Dh=(2.*pow(sin(qDegreesToRadians(theta_e/2.)),2)*pow(sin(qDegreesToRadians(phi_e)),2))/pow(1+Mach*cos(qDegreesToRadians(theta_e))*(1.+(Mach-Mach*EddyMach)*cos(qDegreesToRadians(phi_e))),2);
return Dh;
}

double NoiseCalculation::calcDl(double Mach, double theta_e, double phi_e){
double Dl=(2.*pow(sin(qDegreesToRadians(theta_e)),2.)*pow(sin(qDegreesToRadians(phi_e)),2.))/pow((1.+(Mach*cos(qDegreesToRadians(theta_e)))),4.);
return Dl;
}

//Sara
double NoiseCalculation::getInputWindSpeed(int blade, int E, int section, double TSR){
    double windspeed=0;
    SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
        double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
        double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
        double z=lstart;

    QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
foreach(BData * bdata, pbem->m_pBEMData->GetBData()){
        if (z==TSR){
    if(m_parameter->state_ss_us==0){
//steady
    if(m_parameter->shear_check){
//wind shear effect
double hub_height = m_parameter->tower_height+m_parameter->tower_to_hub_distance;
int section_radius = bdata->m_pos.value(section);

double anglesteps;

if (m_parameter->rotation_type==0){
//    angle based
anglesteps=m_parameter->anglesteps;
}else{
//    time based
anglesteps=(m_parameter->timesteps-1)*60.*360./(m_parameter->rot_speed*1000.);
}

int blades_num = bdata->blades;
double angle_between_blades=360./blades_num;
double initial_azimuth=m_parameter->initial_azimuth;
double E_o=initial_azimuth;
E_o=E_o+anglesteps;
double azimuthal=(E_o+angle_between_blades*blade)+E*anglesteps;
double section_height = hub_height+(section_radius*sin(qDegreesToRadians(azimuthal)))*cos(qDegreesToRadians(m_parameter->yaw_angle));
double m_meanWindSpeed = pbem->dlg_windspeed;

if(section_height>100){
    // calculated with power law wind profile. https://en.wikipedia.org/wiki/Wind_profile_power_law
windspeed = m_meanWindSpeed*pow((section_height/m_parameter->shear_height),(1./7.));
}else{
    // calculated with log wind profile. Should not be used with heigth above 100m (see wikipedia)
windspeed = m_meanWindSpeed*log((section_height)/m_parameter->shear_roughness)/log(m_parameter->shear_height/m_parameter->shear_roughness);
}}
else {
//no wind shear effect
        windspeed = bdata->m_Windspeed.value(section);
}}
 else {
//unsteady
//Sara
        double hub_height = m_parameter->tower_height+m_parameter->tower_to_hub_distance;
        int section_radius = bdata->m_pos.value(section);
        CVector windspeed_windfield;

        double anglesteps=g_windFieldModule->getShownWindField()->getNumberOfTimesteps()*60.*360./(m_parameter->rot_speed*1000.);

        int blades_num = bdata->blades;
        double angle_between_blades=360./blades_num;
        double initial_azimuth=m_parameter->initial_azimuth;
        double E_o=initial_azimuth;
        E_o=E_o+anglesteps;
        double azimuthal=(E_o+angle_between_blades*blade)+E*anglesteps;
        double time = azimuthal/anglesteps*g_windFieldModule->getShownWindField()->getLengthOfTimestep();

const double X = section_radius*sin(qDegreesToRadians(m_parameter->yaw_angle));
const double Y = section_radius*cos(qDegreesToRadians(azimuthal));
const double Z = hub_height+section_radius*sin(qDegreesToRadians(azimuthal))*cos(qDegreesToRadians(m_parameter->yaw_angle));

        CVector vec (X,Y,Z);

        windspeed = g_windFieldModule->getShownWindField()->getWindspeed(vec,time,0).VAbs();
//qDebug() << "X: " << X;
//qDebug() << "Y: " << Y;
//qDebug() << "Z: " << Z;
//qDebug() << "time: " << time;
//qDebug() << "vel: " << windspeed;
    }
}
z=z+ldelta;
}
return windspeed;
}

//Sara
double NoiseCalculation::getInputMach(double windspeed, int section, double TSR){
double Mach=0;
SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
    double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
    double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
    double z=lstart;

QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
foreach(BData * bdata, pbem->m_pBEMData->GetBData()){
    if (z==TSR){
    if(m_parameter->state_ss_us==0){
//steady

    if(m_parameter->shear_check){
//wind shear effect
double Vrel2 = pow(windspeed,2);
Mach = pow(Vrel2,0.5)/sqrt(bdata->k_air*bdata->r_air*bdata->temp);
}
else {
//no wind shear effect
Mach = bdata->m_Mach.value(section);
}}
 else {
//unsteady
double Vrel2 = pow(windspeed,2);
Mach = pow(Vrel2,0.5)/sqrt(bdata->k_air*bdata->r_air*bdata->temp);
    }
}
    z=z+ldelta;
}
return Mach;
}

double NoiseCalculation::getInputReynolds(double windspeed, int section, double TSR){
    double Reynolds=0;
    SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
        double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
        double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
        double z=lstart;

    QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
foreach(BData * bdata, pbem->m_pBEMData->GetBData()){
        if (z==TSR){
    if(m_parameter->state_ss_us==0){
//steady

    if(m_parameter->shear_check){
//wind shear effect
Reynolds = pow((pow(windspeed*(1-bdata->m_a_axial.value(section)),2)+pow(windspeed*bdata->m_lambda_local.value(section)*(1+bdata->m_a_tangential.value(section)),2)),0.5)*bdata->m_c_local.value(section)/bdata->visc;
}
else {
//no wind shear effect
Reynolds = bdata->m_Reynolds.value(section);
}}
 else {
//unsteady
Reynolds = pow((pow(windspeed*(1-bdata->m_a_axial.value(section)),2)+pow(windspeed*bdata->m_lambda_local.value(section)*(1+bdata->m_a_tangential.value(section)),2)),0.5)*bdata->m_c_local.value(section)/bdata->visc;
    }
}
z=z+ldelta;
}
return Reynolds;
}

//vector for blade
void NoiseCalculation::calculateqs3d_graphics(int blade, int E, double TSR) {
  QList<NoiseOpPoint*> noiseOpPoints = m_parameter->prepareNoiseOpPointList();

SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
    double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
    double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
    double z=lstart;

QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
foreach(BData * bdata, pbem->m_pBEMData->GetBData()){

    if (z==TSR){

int number_of_segments = pbem->dlg_elements;
    double approaxing_wind_speed = m_parameter->u_wind_speed;

        double blade_pitch=pbem->m_pctrlFixedPitch->getValue();
        double rho = pbem->dlg_rho;
        double lambda = pbem->dlg_lambda;
        int mpos_size = pbem->dlg_elements;
        double finalradius = bdata->m_pos.value(mpos_size-1);
        double nom_tg_speed = approaxing_wind_speed*lambda;
        double omega = nom_tg_speed/finalradius;
        double c_0_le = 34000;
        double c_const_rd_le = 19./6.;
        double d_const_rd_le = 65.95;
        double c_const_vk_le = 7./3.;
        double d_const_vk_le = 58.4;
        double c_const_le=0;
        double d_const_le=0;
        double u_le;
        double c_le;
        double I_le;
        double lambda_le;
        double beta_le;
        double beta_le_rotor;
        double L_le;
        double K_le;
        double S_le;
        double LFC_le;
        double LFC_le_rotor;

        //quasi 3d rotor calculation
        double rot_speed = m_parameter->rot_speed;
        double omega_rotor;
        double XLT = m_parameter->obs_x_pos_rotor;
        double YLT = m_parameter->obs_y_pos_rotor;
        double ZLT = m_parameter->obs_z_pos_rotor;
        double HR=pbem->m_pBlade->m_HubRadius; //hub radius
        double H=m_parameter->tower_height;//tower height
        double Y=m_parameter->yaw_angle;//yaw angle

        //definitions
        double axial_ind_fact[number_of_segments];
        double axial_ind_fact_n[number_of_segments];
        double axial_velocity[number_of_segments];
        double tangential_speed[number_of_segments];
        double resultant_local_speed[number_of_segments];
        double chord[number_of_segments];
        double Reynolds[number_of_segments];
        double Reynolds_rotor[number_of_segments];
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
        double theta_BEM[number_of_segments];
        double r_R[number_of_segments];
        double c_Rx[number_of_segments];
        double D_starred_C_HT[number_of_segments];
        double D_starred_HT[number_of_segments];
        double D_starred_C_N[number_of_segments];
        double D_starred_C_HT_rotor[number_of_segments];
        double D_starred_HT_rotor[number_of_segments];
        double D_starred_C_N_rotor[number_of_segments];
        double D_starred_N[number_of_segments];
        double D_starred_N_rotor[number_of_segments];
        double Dh[number_of_segments];
        double Dl[number_of_segments];
        double Dl_le[number_of_segments];
        double Dh_rotor[number_of_segments];
        double Dl_rotor[number_of_segments];
        double Dl_le_rotor[number_of_segments];
        double aux0_le[number_of_segments];
        double aux0_le_rotor[number_of_segments];
        double Mach[number_of_segments];
        double Mach_rotor[number_of_segments];
        double corr_fact[number_of_segments];
        double D_starred_HT_S[number_of_segments];
        double D_starred_HT_P[number_of_segments];
        double D_starred_N_S[number_of_segments];
        double D_starred_N_P[number_of_segments];
        double D_starred_HT_S_rotor[number_of_segments];
        double D_starred_HT_P_rotor[number_of_segments];
        double D_starred_N_S_rotor[number_of_segments];
        double D_starred_N_P_rotor[number_of_segments];
        double L[number_of_segments];
        double SwAlpha[number_of_segments];
        double SwAlpha_1[number_of_segments];
        double SwAlpha_2[number_of_segments];
        double gamma[number_of_segments];
        double gamma0[number_of_segments];
        double beta[number_of_segments];
        double beta0[number_of_segments];
        double gamma0_gamma_min[number_of_segments];
        double gamma0_gamma_plus[number_of_segments];
        double gamma_rotor[number_of_segments];
        double beta_rotor[number_of_segments];
        double beta0_rotor[number_of_segments];
        double gamma0_gamma_min_rotor[number_of_segments];
        double gamma0_gamma_plus_rotor[number_of_segments];
        double K1[number_of_segments];
        double K2[number_of_segments];
        double K1_rotor[number_of_segments];
        double K2_rotor[number_of_segments];
        double EddyMach_calc[number_of_segments];
        double EddyMach_calc_rotor[number_of_segments];
        double dist_obs[number_of_segments];
        double dist_obs_rotor[number_of_segments];
        double D_starred_S[number_of_segments];
        double D_starred_P[number_of_segments];
        double D_starred_S_rotor[number_of_segments];
        double D_starred_P_rotor[number_of_segments];
        double first_term_Dh_S[number_of_segments];
        double first_term_Dl_S[number_of_segments];
        double first_term_Dh_P[number_of_segments];
        double first_term_Dh_S_rotor[number_of_segments];
        double first_term_Dl_S_rotor[number_of_segments];
        double first_term_Dh_P_rotor[number_of_segments];
        double St1[number_of_segments];
        double St1_rotor[number_of_segments];
        double St2[number_of_segments];
        double St2_rotor[number_of_segments];
        double b0[number_of_segments];
        double b0_rotor[number_of_segments];
        double B_min_b0[number_of_segments];
        double B_min_b0_rotor[number_of_segments];
        double B_max_b0[number_of_segments];
        double B_max_b0_rotor[number_of_segments];
        double BR_b0[number_of_segments];
        double BR_b0_rotor[number_of_segments];
        double RCmod[number_of_segments];
        double RCmod_rotor[number_of_segments];
        double ao_Rc[number_of_segments];
        double ao_Rc_rotor[number_of_segments];
        double A_min_ao[number_of_segments];
        double A_max_ao[number_of_segments];
        double A_min_ao_rotor[number_of_segments];
        double A_max_ao_rotor[number_of_segments];
        double K1_3[number_of_segments];
        double K1_3_rotor[number_of_segments];
        double AR_ao[number_of_segments];
        double AR_ao_rotor[number_of_segments];
        double St1_bar_rotor[number_of_segments];
        double Re_disp_thick[number_of_segments];
        double Re_disp_thick_rotor[number_of_segments];
        double delta_K1[number_of_segments];
        double delta_K1_rotor[number_of_segments];
        double b[number_of_segments];
        double a[number_of_segments];
        double XRS_rotor[number_of_segments];
        double YRS_rotor[number_of_segments];
        double ZRS_rotor[number_of_segments];
        double XRS[number_of_segments];
        double YRS[number_of_segments];
        double ZRS[number_of_segments];
        double XRT[number_of_segments];
        double YRT[number_of_segments];
        double ZRT[number_of_segments];
        double XRT_rotor[number_of_segments];
        double YRT_rotor[number_of_segments];
        double ZRT_rotor[number_of_segments];
        double theta_e[number_of_segments];
        double theta_e_rotor[number_of_segments];
        double phi_e[number_of_segments];
        double phi_e_rotor[number_of_segments];
        double calc_int_a[number_of_segments];
        double calc_int_a_rotor[number_of_segments];
        double r_e[number_of_segments];
        double r_e_rotor[number_of_segments];
        double r_1[number_of_segments];
        double c_1[number_of_segments];
        double r_0[number_of_segments];
        double c_0[number_of_segments];
        double local_twist[number_of_segments];
        double r_e_le[number_of_segments];
        double aux1_le[number_of_segments];
        double aux4_le[number_of_segments];
        double aux5_le[number_of_segments];
        double r_e_le_rotor[number_of_segments];
        double aux1_le_rotor[number_of_segments];
        double aux4_le_rotor[number_of_segments];
        double aux5_le_rotor[number_of_segments];

        double DStarXFoilS[number_of_segments];
        double DStarXFoilP[number_of_segments];

        double vel[number_of_segments];
        double vel_rotor[number_of_segments];

        double aux_m_SPLadB3d=0;
        double aux_m_SPLsdB3d=0;
        double aux_m_SPLpdB3d=0;
        double aux_m_SPLdB3d=0;
        double aux_m_SPLdBAW3d=0;
        double aux_m_SPLdBBW3d=0;
        double aux_m_SPLdBCW3d=0;
        double aux_m_SPL_LEdB3d=0;
        double aux_m_SPL_LEdBAW3d=0;
        double aux_m_SPL_LEdBBW3d=0;
        double aux_m_SPL_LEdBCW3d=0;

        double aux_m_SPLadB3d_rotor=0;
        double aux_m_SPLsdB3d_rotor=0;
        double aux_m_SPLpdB3d_rotor=0;
        double aux_m_SPLdB3d_rotor=0;
        double aux_m_SPLdBAW3d_rotor=0;
        double aux_m_SPLdBBW3d_rotor=0;
        double aux_m_SPLdBCW3d_rotor=0;
        double aux_m_SPL_LEdB3d_rotor=0;
        double aux_m_SPL_LEdBAW3d_rotor=0;
        double aux_m_SPL_LEdBBW3d_rotor=0;
        double aux_m_SPL_LEdBCW3d_rotor=0;

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
        double splog_LedB=0;
        double splog_LedBAW=0;
        double splog_LedBBW=0;
        double splog_LedBCW=0;

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

        for (int i = 0; i < number_of_segments; ++i) {

            int w=FREQUENCY_TABLE_SIZE;

                double SPL_LedBCW[w];
                double SPL_LedBBW[w];
                double SPL_LedBAW[w];
                double SPL_LedB[w];
                double SPL_C[w];
                double SPL_B[w];
                double SPL_A[w];
                double SPL_dB[w];
                double SPL_P[w];
                double SPL_S[w];
                double SPL_alpha[w];
                double dBC_P[w];
                double dBB_P[w];
                double dBA_P[w];
                double SPL_dB_P[w];
                double SPL_LedBCW_rotor[w];
                double SPL_LedBBW_rotor[w];
                double SPL_LedBAW_rotor[w];
                double SPL_LedB_rotor[w];
                double SPL_C_rotor[w];
                double SPL_B_rotor[w];
                double SPL_A_rotor[w];
                double SPL_dB_rotor[w];
                double SPL_P_rotor[w];
                double SPL_S_rotor[w];
                double SPL_alpha_rotor[w];
                double dBC_P_rotor[w];
                double dBB_P_rotor[w];
                double dBA_P_rotor[w];
                double SPL_dB_P_rotor[w];
                double A_a_P[w];
                double A_a_P_rotor[w];
                double A_max_P[w];
                double A_min_P[w];
                double A_max_P_rotor[w];
                double A_min_P_rotor[w];
                double a_P[w];
                double a_P_rotor[w];
                double Sts[w];
                double Sts_rotor[w];
                double b_alpha[w];
                double b_alpha_rotor[w];
                double B_min[w];
                double B_max[w];
                double B_min_rotor[w];
                double B_max_rotor[w];
                double B_b[w];
                double B_b_rotor[w];
                double a_alpha[w];
                double a_alpha_rotor[w];
                double A_min_alpha[w];
                double A_max_alpha[w];
                double A_min_alpha_rotor[w];
                double A_max_alpha_rotor[w];
                double Alin_a[w];
                double Alin_a_rotor[w];
                double SPL_alpha_min0[w];
                double SPL_alpha_big0[w];
                double dBA_alpha_min0[w];
                double dBA_alpha_big0[w];
                double dBB_alpha_min0[w];
                double dBC_alpha_big0[w];
                double dBC_alpha_min0[w];
                double dBB_alpha_big0[w];
                double SPL_alpha_min0_rotor[w];
                double SPL_alpha_big0_rotor[w];
                double dBA_alpha_min0_rotor[w];
                double dBA_alpha_big0_rotor[w];
                double dBB_alpha_min0_rotor[w];
                double dBC_alpha_big0_rotor[w];
                double dBC_alpha_min0_rotor[w];
                double dBB_alpha_big0_rotor[w];
                double Stp_P[w];
                double Stp_P_rotor[w];
                double Sts_St1_bar[w];
                double Sts_St1_bar_rotor[w];
                double dBC_S[w];
                double dBB_S[w];
                double dBA_S[w];
                double dBA_S_rotor[w];
                double dBC_S_rotor[w];
                double dBB_S_rotor[w];
                double SPL_dB_S[w];
                double SPL_dB_S_rotor[w];
                double A_a_S[w];
                double A_a_S_rotor[w];
                double A_max_S[w];
                double A_min_S[w];
                double A_max_S_rotor[w];
                double A_min_S_rotor[w];
                double a_S[w];
                double a_S_rotor[w];
                double St1_bar[w];

                vel_rotor[i]=getInputWindSpeed(blade, E, i, TSR);
                Reynolds_rotor[i]=getInputReynolds(vel_rotor[i], i, TSR);
                Mach_rotor[i]=getInputMach(vel_rotor[i], i, TSR);

                vel[i]=bdata->m_Windspeed.value(i);
                Reynolds[i]=bdata->m_Reynolds.value(i);
                Mach[i]=bdata->m_Mach.value(i);

     for (int j = 0; j < w; ++j) {
            // definitions
            axial_ind_fact[i] = bdata->m_a_axial.value(i);
           if (i==(number_of_segments-1)){axial_ind_fact_n[i] = bdata->m_a_axial.value(i);}else {axial_ind_fact_n[i] = bdata->m_a_axial.value(i+1);}

            if (i<number_of_segments/2) {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact[i]);}
            else {axial_velocity[i] = approaxing_wind_speed*(1.f-axial_ind_fact_n[i]);}

            tangential_speed[i] = omega*bdata->m_pos.value(i)*(1.+bdata->m_a_tangential.value(i));
            resultant_local_speed[i] = qSqrt(pow(axial_velocity[i],2)+pow(tangential_speed[i],2));
            chord[i] = bdata->m_c_local.value(i);

            Reynolds_BEM[i]=Reynolds[i];
            Reynolds_polar[i]=noiseOpPoints[i]->getReynolds();
            Reynolds_error[i]=qFabs(Reynolds_polar[i]-Reynolds_BEM[i])/Reynolds_BEM[i]*100.;
            Mach_polar[i]=noiseOpPoints[i]->getMach();

            Mach_BEM[i] = Mach[i];
            Mach_error[i]=qFabs(Mach_polar[i]-Mach_BEM[i])/Mach_BEM[i]*100.;
            alpha_BEM[i] = bdata->m_alpha.value(i);
            alpha[i]=alpha_BEM[i];
            theta_BEM[i] = bdata->m_theta.value(i);
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
            if (r_R[i] == r_R1) {c_Rx[i] = c_R1;}
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

if (alpha[i]==0.){
D_starred_HT_S[i]=D_starred_HT[i];
D_starred_HT_P[i]=D_starred_HT[i];
D_starred_N_S[i]=D_starred_N[i];
D_starred_N_P[i]=D_starred_N[i];
}
else{
//alpha !=0 pressure side
if (alpha[i]!=0.){
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
DStarXFoilS[i]=m_DStarInterpolatedS3d[i];
DStarXFoilP[i]=m_DStarInterpolatedP3d[i];

//heavy tripping

if (Reynolds_rotor[i]>300000){
    D_starred_C_HT_rotor[i]=pow(10,(3.411-1.5397*log10(Reynolds_rotor[i])+0.1059*pow(log10(Reynolds_rotor[i]),2)));
}
else {D_starred_C_HT_rotor[i]=0.0601*(pow(Reynolds_rotor[i],(-0.114)));}

D_starred_HT_rotor[i]=chord[i]*D_starred_C_HT_rotor[i];

//natural transition
D_starred_C_N_rotor[i]=pow(10,(3.0187-1.5397*log10(Reynolds_rotor[i])+0.1059*pow(log10(Reynolds_rotor[i]),2)));
D_starred_N_rotor[i]=D_starred_C_N_rotor[i]*chord[i];

if (alpha[i]==0.){
D_starred_HT_S_rotor[i]=D_starred_HT_rotor[i];
D_starred_HT_P_rotor[i]=D_starred_HT_rotor[i];
D_starred_N_S_rotor[i]=D_starred_N_rotor[i];
D_starred_N_P_rotor[i]=D_starred_N_rotor[i];
}
else{
//alpha !=0 pressure side
if (alpha[i]!=0){
D_starred_HT_P_rotor[i]=D_starred_HT_rotor[i]*corr_fact[i];
D_starred_N_P_rotor[i]=D_starred_N_rotor[i]*corr_fact[i];
}

//alpha !=0 suction side heavy tripping
if ((alpha[i]>0) & (alpha[i]<=5)){
corr_fact[i]=pow(10,(0.0679*alpha[i]));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
D_starred_HT_S_rotor[i]=D_starred_HT_rotor[i]*corr_fact[i];
}

if ((alpha[i]>5) & (alpha[i]<=12.5)){
corr_fact[i]=0.381*(pow(10,(0.1516*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
D_starred_HT_S_rotor[i]=D_starred_HT_rotor[i]*corr_fact[i];
}

if ((alpha[i]>12.5) & (alpha[i]<=25)){
corr_fact[i]=14.296*(pow(10,(0.0258*alpha[i])));
D_starred_HT_S[i]=D_starred_HT[i]*corr_fact[i];
D_starred_HT_S_rotor[i]=D_starred_HT_rotor[i]*corr_fact[i];
}

//alpha !=0 suction side natural transition
if ((alpha[i]>0) & (alpha[i]<=7.5)){
corr_fact[i]=pow(10.,(0.0679*alpha[i]));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
D_starred_N_S_rotor[i]=D_starred_N_rotor[i]*corr_fact[i];
}

if ((alpha[i]>7.5) & (alpha[i]<=12.5)){
corr_fact[i]=0.0162*(pow(10.,(0.3066*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
D_starred_N_S_rotor[i]=D_starred_N_rotor[i]*corr_fact[i];
}

if ((alpha[i]>12.5) & (alpha[i]<=25)){
corr_fact[i]=54.42*(pow(10.,(0.0258*alpha[i])));
D_starred_N_S[i]=D_starred_N[i]*corr_fact[i];
D_starred_N_S_rotor[i]=D_starred_N_rotor[i]*corr_fact[i];
}}

//For D* Xfoil
DStarXFoilS[i]=m_DStarInterpolatedS3d[i];
DStarXFoilP[i]=m_DStarInterpolatedP3d[i];

//Length of Wetted Trailing Edge
if (i==(number_of_segments-1)){L[i]=0;}else{
L[i]=bdata->m_pos.value(i+1)-bdata->m_pos.value(i);}

//Calculate the Switching Angle
SwAlpha_1[i]=23.43*Mach[i]+4.651;
SwAlpha_2[i]=12.5;

if (SwAlpha_1[i]<SwAlpha_2[i]){SwAlpha[i]=SwAlpha_1[i];}
else {SwAlpha[i]=SwAlpha_2[i];}

double EddyMach = m_parameter->eddyConvectionMach;

gamma[i]=27.094*Mach[i]+3.32;
gamma_rotor[i]=27.094*Mach_rotor[i]+3.32;

gamma0[i]=SwAlpha_1[i];

beta[i]=72.65*Mach[i]+10.74;
beta_rotor[i]=72.65*Mach_rotor[i]+10.74;

beta0[i]=-34.19*Mach[i]-13.82;
beta0_rotor[i]=-34.19*Mach_rotor[i]-13.82;

gamma0_gamma_min[i]=gamma0[i]-gamma[i];
gamma0_gamma_plus[i]=gamma0[i]+gamma[i];

gamma0_gamma_min_rotor[i]=gamma0[i]-gamma_rotor[i];
gamma0_gamma_plus_rotor[i]=gamma0[i]+gamma_rotor[i];

if (Reynolds[i]<247000)
{K1[i]=-4.31*log10(Reynolds[i])+156.3;}
else if (Reynolds[i]>800000)
{K1[i]=128.5;}
else {K1[i]=-9.*log10(Reynolds[i])+181.6;}

if (alpha[i]<gamma0_gamma_min[i])
{K2[i]=K1[i]-1000.;}
else if (alpha[i]>gamma0_gamma_plus[i])
{K2[i]=K1[i]-12.;}
else
{K2[i]=K1[i]+(sqrt(pow(beta[i],2)-pow((beta[i]/gamma[i]),2)*pow((alpha[i]-gamma0[i]),2)))+beta0[i];}

if (Reynolds_rotor[i]<247000)
{K1_rotor[i]=-4.31*log10(Reynolds_rotor[i])+156.3;}
else if (Reynolds_rotor[i]>800000)
{K1_rotor[i]=128.5;}
else {K1_rotor[i]=-9.*log10(Reynolds_rotor[i])+181.6;}

if (alpha[i]<gamma0_gamma_min[i])
{K2_rotor[i]=K1_rotor[i]-1000.;}
else if (alpha[i]>gamma0_gamma_plus[i])
{K2_rotor[i]=K1_rotor[i]-12.;}
else
{K2_rotor[i]=K1_rotor[i]+(sqrt(pow(beta[i],2)-pow((beta[i]/gamma[i]),2)*pow((alpha[i]-gamma0[i]),2)))+beta0[i];}


double EddyMach_perc=EddyMach;

EddyMach_calc[i]=Mach[i]*EddyMach_perc;
EddyMach_calc_rotor[i]=Mach_rotor[i]*EddyMach_perc;

//delta starred type, if natural transition or heavy-tripping
if (m_parameter->dstar_type==0){
//XFoil calculation
    D_starred_S[i]=DStarXFoilS[i];
    D_starred_P[i]=DStarXFoilP[i];
    D_starred_S_rotor[i]=DStarXFoilS[i];
    D_starred_P_rotor[i]=DStarXFoilP[i];
}
else if (m_parameter->dstar_type==1){
//    BPM calculation
    FoilPolarDlg *pFoilPolarDlg = (FoilPolarDlg *) g_mainFrame->m_pctrlXDirectWidget;
        double TopTrip=pFoilPolarDlg->m_XTopTr;
        double BotTrip=pFoilPolarDlg->m_XBotTr;
    if((TopTrip==0.) & (BotTrip==0.)) {
    //        natural transition
        D_starred_S[i]=D_starred_N_S[i];
        D_starred_P[i]=D_starred_N_P[i];
        D_starred_S_rotor[i]=D_starred_N_S_rotor[i];
        D_starred_P_rotor[i]=D_starred_N_P_rotor[i];
    }
    else {
    //heavy tripping
        D_starred_S[i]=D_starred_HT_S[i];
        D_starred_P[i]=D_starred_HT_P[i];
        D_starred_S_rotor[i]=D_starred_HT_S_rotor[i];
        D_starred_P_rotor[i]=D_starred_HT_P_rotor[i];
}}
else if (m_parameter->dstar_type==2){
//user
    NoiseParameter *pNoiseParameter = (NoiseParameter *) g_mainFrame->m_pSimuWidget;
    D_starred_S[i]=pNoiseParameter->D_starred_S_user[i];
    D_starred_P[i]=pNoiseParameter->D_starred_P_user[i];
    D_starred_S_rotor[i]=pNoiseParameter->D_starred_S_user[i];
    D_starred_P_rotor[i]=pNoiseParameter->D_starred_P_user[i];
}
//Sara

//double B=0;
double XB=0;
double YB=0;
double ZB=0;

//    re phi_e and theta_e calculation p 77 C_Project_Log_Text_Jan_16.pdf
//    Input X e , Y e , Z e
//    Attribute their respective values to X B , Y B , Z B
    XB=m_parameter->obs_x_pos;
    YB=m_parameter->obs_y_pos;
    ZB=m_parameter->obs_z_pos;

//rotor
    omega_rotor = 2.*M_PI*rot_speed/60.; //rotor
    double XUT=XLT;
    double YUT=YLT;
    double ZUT=ZLT-H; //Sara it was +H

    double XYB=XUT*cos(qDegreesToRadians(Y))+YUT*sin(qDegreesToRadians(Y));
    double YYB=XUT*-sin(qDegreesToRadians(Y))+YUT*cos(qDegreesToRadians(Y));
    double ZYB=ZUT;

    double TTH=m_parameter->tower_to_hub_distance;//tower to hub distance
    int blades_num = bdata->blades;
    double anglesteps;

    if (m_parameter->rotation_type==0){
    //    angle based
    anglesteps=m_parameter->anglesteps;
    }else{
    //    time based
    anglesteps=m_parameter->timesteps*60.*360./(m_parameter->rot_speed*1000.);
    }

    double angle_between_blades=360./blades_num;
    double initial_azimuth=m_parameter->initial_azimuth; //initial azimuth;
    double E_o=initial_azimuth;
//    double E_f;
    int number_of_rotations;
    E_o=E_o+anglesteps;
    double azimuthal=(E_o+angle_between_blades*blade)+E*anglesteps;

    double XHF=XYB+TTH;//Sara it was -TTH
    double YHF=YYB-HR*sin(qDegreesToRadians(azimuthal)); //Sara it was YYB
    double ZHF=ZYB-HR*(1+cos(qDegreesToRadians(azimuthal))); //Sara it was ZYB

    if (m_parameter->rotation_type==0){
    //    angle based
    number_of_rotations = m_parameter->number_loops;
    anglesteps=m_parameter->anglesteps;
    }else{
    //    time based
    number_of_rotations = m_parameter->time/(60./m_parameter->rot_speed);
    anglesteps=m_parameter->timesteps*60.*360./(m_parameter->rot_speed*1000.);
    }

    double XRR=XHF;
    double YRR=YHF*cos(qDegreesToRadians(azimuthal))+XHF*sin(qDegreesToRadians(azimuthal));
    double ZRR=YHF*(-sin(qDegreesToRadians(azimuthal)))+ZHF*cos(qDegreesToRadians(azimuthal));

    double XB_rotor=XRR;
    double YB_rotor=YRR;
    double ZB_rotor=ZRR; //Sara it was ZRR+HR
//rotor

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

//    the angle a is the total angle between the YB ZB blade reference system plane and the local midsection chord line p 75 handout
    a[i]=local_twist[i]+blade_pitch;

    XRS[i]=calcXRS(a[i],XB,YB);
    YRS[i]=calcYRS(a[i],XB,YB);
    ZRS[i]=calcZRS(ZB,r_0[i],r_1[i]);

//rotor
    XRS_rotor[i]=calcXRS(a[i],XB_rotor,YB_rotor);
    YRS_rotor[i]=calcYRS(a[i],XB_rotor,YB_rotor);
    ZRS_rotor[i]=calcZRS(ZB_rotor,r_0[i],r_1[i]);
//rotor

    calc_int_a[i]=calcInt_a(YRS[i],c_0[i],c_1[i]);

    calc_int_a_rotor[i]=calcInt_a(YRS_rotor[i],c_0[i],c_1[i]);//rotor

    XRT[i]=calcXRT(XRS[i]);
    YRT[i]=calcYRT(b[i],calc_int_a[i],ZRS[i]);
    ZRT[i]=calcZRT(b[i],calc_int_a[i],ZRS[i]);

//rotor
    XRT_rotor[i]=calcXRT(XRS_rotor[i]);
    YRT_rotor[i]=calcYRT(b[i],calc_int_a_rotor[i],ZRS_rotor[i]);
    ZRT_rotor[i]=calcZRT(b[i],calc_int_a_rotor[i],ZRS_rotor[i]);
//rotor

    r_e[i]=calcR_e(XRT[i],YRT[i],ZRT[i]);
    theta_e[i]=calcTheta_e(YRT[i],ZRT[i]);
    phi_e[i]=calcPhi_e(XRT[i], ZRT[i]);
    dist_obs[i]=r_e[i];

//rotor
    r_e_rotor[i]=calcR_e(XRT_rotor[i],YRT_rotor[i],ZRT_rotor[i]);
    theta_e_rotor[i]=calcTheta_e(YRT_rotor[i],ZRT_rotor[i]);
    phi_e_rotor[i]=calcPhi_e(XRT_rotor[i], ZRT_rotor[i]);
    dist_obs_rotor[i]=r_e_rotor[i];
//rotor

alpha_error[i]=qFabs(alpha_polar[i]-alpha_BEM[i])/alpha_BEM[i]*100.;

St1[i]=0.02*(pow(Mach[i],-0.6));
St1_rotor[i]=0.02*(pow(Mach_rotor[i],-0.6));

if(alpha[i]<1.33){
    St2[i]=St1[i];
    St2_rotor[i]=St1_rotor[i];}
else if(alpha[i]>12.5){
    St2[i]=St1[i]*4.72;
    St2_rotor[i]=St1_rotor[i]*4.72;}
else {St2[i]=St1[i]*pow(10.,(0.0054*pow((alpha[i]-1.33),2)));
St2_rotor[i]=St1_rotor[i]*pow(10.,(0.0054*pow((alpha[i]-1.33),2)));}

if (Reynolds[i]<95200)
{b0[i]= 0.3;}
else if (Reynolds[i]>857000)
{b0[i]= 0.56;}
else {b0[i]=-4.48*pow(10,-13)*(pow((Reynolds[i]-857000.),2)+0.56);}

if (Reynolds_rotor[i]<95200)
{b0_rotor[i]= 0.3;}
else if (Reynolds_rotor[i]>857000)
{b0_rotor[i]= 0.56;}
else {b0_rotor[i]=-4.48*pow(10,-13)*(pow((Reynolds_rotor[i]-857000.),2)+0.56);}

if (b0[i]<0.13)
{B_min_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
else if (b0[i]>0.145)
{B_min_b0[i]=-817.81*pow(b0[i],3)+335.21*pow(b0[i],2)-135.024*b0[i]+10.619;}
else{B_min_b0[i]=-83.607*b0[i]+8.138;}

if (b0_rotor[i]<0.13)
{B_min_b0_rotor[i]=sqrt(16.888-886.788*pow(b0_rotor[i],2))-4.109;}
else if (b0_rotor[i]>0.145)
{B_min_b0_rotor[i]=-817.81*pow(b0_rotor[i],3)+335.21*pow(b0[i],2)-135.024*b0_rotor[i]+10.619;}
else{B_min_b0_rotor[i]=-83.607*b0_rotor[i]+8.138;}

if (b0[i]<0.1)
{B_max_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
else if (b0[i]>0.187)
{B_max_b0[i]=-80.541*pow(b0[i],3)+44.174*pow(b0[i],2)-39.381*b0[i]+2.344;}
else {B_max_b0[i]=-31.33*b0[i]+1.854;}

if (b0_rotor[i]<0.1)
{B_max_b0_rotor[i]=sqrt(16.888-886.788*pow(b0_rotor[i],2))-4.109;}
else if (b0_rotor[i]>0.187)
{B_max_b0_rotor[i]=-80.541*pow(b0_rotor[i],3)+44.174*pow(b0_rotor[i],2)-39.381*b0_rotor[i]+2.344;}
else {B_max_b0_rotor[i]=-31.33*b0_rotor[i]+1.854;}

BR_b0[i]=(-20-B_min_b0[i])/(B_max_b0[i]-B_min_b0[i]);
BR_b0_rotor[i]=(-20-B_min_b0_rotor[i])/(B_max_b0_rotor[i]-B_min_b0_rotor[i]);

RCmod[i]=3*Reynolds[i];
RCmod_rotor[i]=3*Reynolds_rotor[i];

if (RCmod[i]<95200){ao_Rc[i]=0.57;}
else if (RCmod[i]>857000){ao_Rc[i]=1.13;}
else {ao_Rc[i]=-9.57*pow(10,-13)*(pow((RCmod[i]-857000),2)+1.13);}

if (RCmod_rotor[i]<95200){ao_Rc_rotor[i]=0.57;}
else if (RCmod_rotor[i]>857000){ao_Rc_rotor[i]=1.13;}
else {ao_Rc_rotor[i]=-9.57*pow(10,-13)*(pow((RCmod_rotor[i]-857000),2)+1.13);}

if(ao_Rc[i]<0.204){A_min_ao[i]=sqrt(67.552-886.788*pow(ao_Rc[i],2))-8.219;}
else if (ao_Rc[i]>0.244){A_min_ao[i]=-142.795*pow(ao_Rc[i],3)+103.656*pow(ao_Rc[i],2)-57.757*ao_Rc[i]+6.006;}
else {A_min_ao[i]=-32.665*ao_Rc[i]+3.981;}

if(ao_Rc_rotor[i]<0.204){A_min_ao_rotor[i]=sqrt(67.552-886.788*pow(ao_Rc_rotor[i],2))-8.219;}
else if (ao_Rc_rotor[i]>0.244){A_min_ao_rotor[i]=-142.795*pow(ao_Rc_rotor[i],3)+103.656*pow(ao_Rc_rotor[i],2)-57.757*ao_Rc_rotor[i]+6.006;}
else {A_min_ao_rotor[i]=-32.665*ao_Rc_rotor[i]+3.981;}

if (ao_Rc[i]<0.13){A_max_ao[i]=sqrt(67.552-886.788*pow(ao_Rc[i],2))-8.219;}
else if(ao_Rc[i]>0.321){A_max_ao[i]=-4.669*pow(ao_Rc[i],3)+3.491*pow(ao_Rc[i],2)-16.699*ao_Rc[i]+1.149;}
else {A_max_ao[i]=-15.901*ao_Rc[i]+1.098;}

if (ao_Rc_rotor[i]<0.13){A_max_ao_rotor[i]=sqrt(67.552-886.788*pow(ao_Rc_rotor[i],2))-8.219;}
else if(ao_Rc_rotor[i]>0.321){A_max_ao_rotor[i]=-4.669*pow(ao_Rc_rotor[i],3)+3.491*pow(ao_Rc_rotor[i],2)-16.699*ao_Rc_rotor[i]+1.149;}
else {A_max_ao_rotor[i]=-15.901*ao_Rc_rotor[i]+1.098;}

K1_3[i]=K1[i]-3.;
K1_3_rotor[i]=K1_rotor[i]-3.;

AR_ao[i]=(-20-A_min_ao[i])/(A_max_ao[i]-A_min_ao[i]);
AR_ao_rotor[i]=(-20-A_min_ao_rotor[i])/(A_max_ao_rotor[i]-A_min_ao_rotor[i]);

St1_bar[i]=(St1[i]+St2[i])/2.;
St1_bar_rotor[i]=(St1_rotor[i]+St2_rotor[i])/2.;

Re_disp_thick[i]=rho*vel[i]*D_starred_P[i]/(0.0000178);
Re_disp_thick_rotor[i]=rho*vel_rotor[i]*D_starred_P_rotor[i]/(0.0000178);

if (Re_disp_thick[i]>5000){delta_K1[i]=0;}
else {delta_K1[i]=alpha[i]*(1.43*log10(Re_disp_thick[i])-5.29);}

if (Re_disp_thick_rotor[i]>5000){delta_K1_rotor[i]=0;}
else {delta_K1_rotor[i]=alpha[i]*(1.43*log10(Re_disp_thick_rotor[i])-5.29);}

    if (Reynolds[i]<95200)
    {b0[i]= 0.3;}
    else if (Reynolds[i]>857000)
    {b0[i]= 0.56;}
    else {b0[i]=-4.48*pow(10.,-13)*(pow((Reynolds[i]-857000.),2)+0.56);}

    if (Reynolds_rotor[i]<95200)
    {b0_rotor[i]= 0.3;}
    else if (Reynolds_rotor[i]>857000)
    {b0_rotor[i]= 0.56;}
    else {b0_rotor[i]=-4.48*pow(10.,-13)*(pow((Reynolds_rotor[i]-857000.),2)+0.56);}

    if (b0[i]<0.13)
    {B_min_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
    else if (b0[i]>0.145)
    {B_min_b0[i]=-817.81*pow(b0[i],3)+335.21*pow(b0[i],2)-135.024*b0[i]+10.619;}
    else{B_min_b0[i]=-83.607*b0[i]+8.138;}

    if (b0_rotor[i]<0.13)
    {B_min_b0_rotor[i]=sqrt(16.888-886.788*pow(b0_rotor[i],2))-4.109;}
    else if (b0_rotor[i]>0.145)
    {B_min_b0_rotor[i]=-817.81*pow(b0_rotor[i],3)+335.21*pow(b0_rotor[i],2)-135.024*b0_rotor[i]+10.619;}
    else{B_min_b0_rotor[i]=-83.607*b0_rotor[i]+8.138;}

    if (b0[i]<0.1)
    {B_max_b0[i]=sqrt(16.888-886.788*pow(b0[i],2))-4.109;}
    else if (b0[i]>0.187)
    {B_max_b0[i]=-80.541*pow(b0[i],3)+44.174*pow(b0[i],2)-39.381*b0[i]+2.344;}
    else {B_max_b0[i]=-31.33*b0[i]+1.854;}

    if (b0_rotor[i]<0.1)
    {B_max_b0_rotor[i]=sqrt(16.888-886.788*pow(b0_rotor[i],2))-4.109;}
    else if (b0_rotor[i]>0.187)
    {B_max_b0_rotor[i]=-80.541*pow(b0_rotor[i],3)+44.174*pow(b0_rotor[i],2)-39.381*b0_rotor[i]+2.344;}
    else {B_max_b0_rotor[i]=-31.33*b0_rotor[i]+1.854;}

    BR_b0[i]=(-20-B_min_b0[i])/(B_max_b0[i]-B_min_b0[i]);
    BR_b0_rotor[i]=(-20-B_min_b0_rotor[i])/(B_max_b0_rotor[i]-B_min_b0_rotor[i]);

Sts[j]=CENTRAL_BAND_FREQUENCY[j]*D_starred_S[i]/vel[i];
Sts_rotor[j]=CENTRAL_BAND_FREQUENCY[j]*D_starred_S_rotor[i]/vel_rotor[i];

b_alpha[j]=qFabs(log10(Sts[j]/St2[i]));
b_alpha_rotor[j]=qFabs(log10(Sts_rotor[j]/St2_rotor[i]));

    if (b_alpha[j]<0.13)
    {B_min[j]=sqrt(16.888-886.788*pow(b_alpha[j],2));}
    else if(b_alpha[j]>0.145)
    {B_min[j]=-817.81*pow(b_alpha[j],3)+335.21*pow(b_alpha[j],2)-135.024*b_alpha[j]+10.619;}
    else {B_min[j]=-83.607*b_alpha[j]+8.138;}

    if (b_alpha_rotor[j]<0.13)
    {B_min_rotor[j]=sqrt(16.888-886.788*pow(b_alpha_rotor[j],2));}
    else if(b_alpha_rotor[j]>0.145)
    {B_min_rotor[j]=-817.81*pow(b_alpha_rotor[j],3)+335.21*pow(b_alpha_rotor[j],2)-135.024*b_alpha_rotor[j]+10.619;}
    else {B_min_rotor[j]=-83.607*b_alpha_rotor[j]+8.138;}

    if (b_alpha[j]<0.1)
    {B_max[j]=sqrt(16.888-886.788*pow(b_alpha[j],2))-4.109;}
    else if(b_alpha[j]>0.187)
    {B_max[j]=-80.541*pow(b_alpha[j],3)+44.174*pow(b_alpha[j],2)-39.381*b_alpha[j]+2.344;}
    else {B_max[j]=-31.33*b_alpha[j]+1.854;}

    if (b_alpha_rotor[j]<0.1)
    {B_max_rotor[j]=sqrt(16.888-886.788*pow(b_alpha_rotor[j],2))-4.109;}
    else if(b_alpha_rotor[j]>0.187)
    {B_max_rotor[j]=-80.541*pow(b_alpha_rotor[j],3)+44.174*pow(b_alpha_rotor[j],2)-39.381*b_alpha_rotor[j]+2.344;}
    else {B_max_rotor[j]=-31.33*b_alpha_rotor[j]+1.854;}

    B_b[j]=B_min[j]+BR_b0[i]*(B_max[j]-B_min[j]);
    if(qIsInf(B_b[j])){B_b[j]=0;}

    B_b_rotor[j]=B_min_rotor[j]+BR_b0_rotor[i]*(B_max_rotor[j]-B_min_rotor[j]);
    if(qIsInf(B_b_rotor[j])){B_b_rotor[j]=0;}
    a_alpha_rotor[j]=qFabs(log10(Sts_rotor[j]/St2_rotor[i]));

    if (a_alpha[j]<0.204)
    {A_min_alpha[j]=sqrt(67.552-886.788*pow(a_alpha[j],2))-8.219;}
    else if(a_alpha[j]>0.244)
    {A_min_alpha[j]=-142.795*pow(a_alpha[j],3)+103.656*pow(a_alpha[j],2)-57.757*a_alpha[j]+6.006;}
    else {A_min_alpha[j]=-32.665*a_alpha[j]+3.981;}

    if (a_alpha_rotor[j]<0.204)
    {A_min_alpha_rotor[j]=sqrt(67.552-886.788*pow(a_alpha_rotor[j],2))-8.219;}
    else if(a_alpha_rotor[j]>0.244)
    {A_min_alpha_rotor[j]=-142.795*pow(a_alpha_rotor[j],3)+103.656*pow(a_alpha_rotor[j],2)-57.757*a_alpha_rotor[j]+6.006;}
    else {A_min_alpha_rotor[j]=-32.665*a_alpha_rotor[j]+3.981;}

    if (a_alpha[j]<0.13)
    {A_max_alpha[j]=sqrt(67.552-886.788*pow(a_alpha[j],2))-8.219;}
    else if(a_alpha[j]>0.321)
    {A_max_alpha[j]=-4.669*pow(a_alpha[j],3)+3.491*pow(a_alpha[j],2)-16.699*a_alpha[j]+1.149;}
    else {A_max_alpha[j]=-15.901*a_alpha[j]+1.098;}

    if (a_alpha_rotor[j]<0.13)
    {A_max_alpha_rotor[j]=sqrt(67.552-886.788*pow(a_alpha_rotor[j],2))-8.219;}
    else if(a_alpha_rotor[j]>0.321)
    {A_max_alpha_rotor[j]=-4.669*pow(a_alpha_rotor[j],3)+3.491*pow(a_alpha_rotor[j],2)-16.699*a_alpha_rotor[j]+1.149;}
    else {A_max_alpha_rotor[j]=-15.901*a_alpha_rotor[j]+1.098;}
    Alin_a[j]=A_min_alpha[j]+AR_ao[i]*(A_max_alpha[j]-A_min_alpha[j]);
    Alin_a_rotor[j]=A_min_alpha_rotor[j]+AR_ao_rotor[i]*(A_max_alpha_rotor[j]-A_min_alpha_rotor[j]);

    dBA_alpha_min0[j]=SPL_alpha_min0[j]+AWeighting[j];
    dBA_alpha_big0[j]=SPL_alpha_big0[j]+AWeighting[j];
    dBB_alpha_min0[j]=SPL_alpha_min0[j]+BWeighting[j];
    dBC_alpha_big0[j]=SPL_alpha_big0[j]+CWeighting[j];
    dBC_alpha_min0[j]=SPL_alpha_min0[j]+CWeighting[j];
    dBB_alpha_big0[j]=SPL_alpha_big0[j]+BWeighting[j];

    dBA_alpha_min0_rotor[j]=SPL_alpha_min0_rotor[j]+AWeighting[j];
    dBA_alpha_big0_rotor[j]=SPL_alpha_big0_rotor[j]+AWeighting[j];
    dBB_alpha_min0_rotor[j]=SPL_alpha_min0_rotor[j]+BWeighting[j];
    dBC_alpha_big0_rotor[j]=SPL_alpha_big0_rotor[j]+CWeighting[j];
    dBC_alpha_min0_rotor[j]=SPL_alpha_min0_rotor[j]+CWeighting[j];
    dBB_alpha_big0_rotor[j]=SPL_alpha_big0_rotor[j]+BWeighting[j];

    St1_bar[j]=(St1[i]+St2[i])/2.;
    St1_bar_rotor[j]=(St1_rotor[i]+St2_rotor[i])/2.;

    a_S[j]=qFabs(log10(Sts[j]/St1_bar[j]));
    a_S_rotor[j]=qFabs(log10(Sts_rotor[j]/St1_bar_rotor[j]));

    if (a_S[j]<0.204){A_min_S[j]=sqrt(67.552-886.788*pow(a_S[j],2))-8.219;}
    else if (a_S[j]>0.244){A_min_S[j]=-142.795*pow(a_S[j],3)+103.656*pow(a_S[j],2)-57.757*a_S[j]+6.006;}
    else {A_min_S[j]=-32.665*a_S[j]+3.981;}

    if (a_S_rotor[j]<0.204){A_min_S_rotor[j]=sqrt(67.552-886.788*pow(a_S_rotor[j],2))-8.219;}
    else if (a_S_rotor[j]>0.244){A_min_S_rotor[j]=-142.795*pow(a_S_rotor[j],3)+103.656*pow(a_S_rotor[j],2)-57.757*a_S_rotor[j]+6.006;}
    else {A_min_S_rotor[j]=-32.665*a_S_rotor[j]+3.981;}

    if (a_S[j]<0.13){A_max_S[j]=sqrt(67.552-886.788*pow(a_S[j],2))-8.219;}
    else if (a_S[j]>0.321){A_max_S[j]=-4.669*pow(a_S[j],3)+3.491*pow(a_S[j],2)-16.699*a_S[j]+1.149;}
    else {A_max_S[j]=-15.901*a_S[j]+1.098;}

    if (a_S_rotor[j]<0.13){A_max_S_rotor[j]=sqrt(67.552-886.788*pow(a_S_rotor[j],2))-8.219;}
    else if (a_S_rotor[j]>0.321){A_max_S_rotor[j]=-4.669*pow(a_S_rotor[j],3)+3.491*pow(a_S_rotor[j],2)-16.699*a_S_rotor[j]+1.149;}
    else {A_max_S_rotor[j]=-15.901*a_S_rotor[j]+1.098;}

    A_a_S[j]=A_min_S[j]+AR_ao[i]*(A_max_S[j]-A_min_S[j]);
    A_a_S_rotor[j]=A_min_S_rotor[j]+AR_ao_rotor[i]*(A_max_S_rotor[j]-A_min_S_rotor[j]);

    dBA_S[j]=SPL_dB_S[j]+AWeighting[j];
    dBB_S[j]=SPL_dB_S[j]+BWeighting[j];
    dBC_S[j]=SPL_dB_S[j]+CWeighting[j];

    dBA_S_rotor[j]=SPL_dB_S_rotor[j]+AWeighting[j];
    dBB_S_rotor[j]=SPL_dB_S_rotor[j]+BWeighting[j];
    dBC_S_rotor[j]=SPL_dB_S_rotor[j]+CWeighting[j];

    Sts_St1_bar[j]=Sts[j]/St1_bar[i];
    Sts_St1_bar_rotor[j]=Sts_rotor[j]/St1_bar_rotor[i];

Stp_P[j]=CENTRAL_BAND_FREQUENCY[j]*D_starred_P[i]/vel[i];
Stp_P_rotor[j]=CENTRAL_BAND_FREQUENCY[j]*D_starred_P_rotor[i]/vel_rotor[i];

    a_P[j]=qFabs(log10(Stp_P[j]/St1[i]));
    a_P_rotor[j]=qFabs(log10(Stp_P_rotor[j]/St1_rotor[i]));

    if(a_P[j]<0.204){A_min_P[j]=sqrt(67.552-886.788*pow(a_P[j],2))-8.219;}
    else if (a_P[j]>0.244){A_min_P[j]=-142.795*pow(a_P[j],3)+103.656*pow(a_P[j],2)-57.757*a_P[j]+6.006;}
    else {A_min_P[j]=-32.665*a_P[j]+3.981;}

    if(a_P_rotor[j]<0.204){A_min_P_rotor[j]=sqrt(67.552-886.788*pow(a_P_rotor[j],2))-8.219;}
    else if (a_P_rotor[j]>0.244){A_min_P_rotor[j]=-142.795*pow(a_P_rotor[j],3)+103.656*pow(a_P_rotor[j],2)-57.757*a_P_rotor[j]+6.006;}
    else {A_min_P_rotor[j]=-32.665*a_P_rotor[j]+3.981;}

    if(a_P[j]<0.13){A_max_P[j]=sqrt(67.552-886.788*pow(a_P[j],2))-8.219;}
    else if (a_P[j]>0.321){A_max_P[j]=-4.669*pow(a_P[j],3)+3.491*pow(a_P[j],2)-16.699*a_P[j]+1.149;}
    else {A_max_P[j]=-15.901*a_P[j]+1.098;}

    if(a_P_rotor[j]<0.13){A_max_P_rotor[j]=sqrt(67.552-886.788*pow(a_P_rotor[j],2))-8.219;}
    else if (a_P_rotor[j]>0.321){A_max_P_rotor[j]=-4.669*pow(a_P_rotor[j],3)+3.491*pow(a_P_rotor[j],2)-16.699*a_P_rotor[j]+1.149;}
    else {A_max_P_rotor[j]=-15.901*a_P_rotor[j]+1.098;}

    A_a_P[j]=A_min_P[j]+AR_ao[i]*(A_max_P[j]-A_min_P[j]);
    A_a_P_rotor[j]=A_min_P_rotor[j]+AR_ao_rotor[i]*(A_max_P_rotor[j]-A_min_P_rotor[j]);

    SPL_A[j]=0;
    SPL_B[j]=0;
    SPL_C[j]=0;
    dBA_P[j]=SPL_dB_P[j]+AWeighting[j];
    dBB_P[j]=SPL_dB_P[j]+BWeighting[j];
    dBC_P[j]=SPL_dB_P[j]+CWeighting[j];
    SPL_LedB[w]=0;
    SPL_LedBAW[w]=0;
    SPL_LedBBW[w]=0;
    SPL_LedBCW[w]=0;

    //rotor
    SPL_A_rotor[j]=0;
    SPL_B_rotor[j]=0;
    SPL_C_rotor[j]=0;
    dBA_P_rotor[j]=SPL_dB_P_rotor[j]+AWeighting[j];
    dBB_P_rotor[j]=SPL_dB_P_rotor[j]+BWeighting[j];
    dBC_P_rotor[j]=SPL_dB_P_rotor[j]+CWeighting[j];
    SPL_LedB_rotor[w]=0;
    SPL_LedBAW_rotor[w]=0;
    SPL_LedBBW_rotor[w]=0;
    SPL_LedBCW_rotor[w]=0;
    //rotor

    u_le=approaxing_wind_speed*100.;
    c_le=100.*chord[i];
    I_le=m_parameter->TurbulenceIntensity;
    lambda_le=100.*m_parameter->IntegralLengthScale;
    r_e_le[i]=100.*dist_obs[i];
    r_e_le_rotor[i]=100.*dist_obs_rotor[i];//rotor

    beta_le=sqrt(1-pow(Mach[i],2));
    beta_le_rotor=sqrt(1-pow(Mach_rotor[i],2));
    L_le=100.*L[i];
    K_le=M_PI*CENTRAL_BAND_FREQUENCY[j]*c_le/u_le;
    S_le=sqrt(pow((2.*M_PI*K_le/(pow(beta_le, 2)))+(pow((1+(2.4*K_le/pow(beta_le,2))),-1)),-1));
    LFC_le = 10.*Mach[i]*pow(S_le*K_le/beta_le,2);
    LFC_le_rotor = 10.*Mach_rotor[i]*pow(S_le*K_le/beta_le,2);

    if(m_parameter->Lowson_type==2){
        c_const_le = c_const_rd_le;
        d_const_le = d_const_rd_le;
     }
    else if(m_parameter->Lowson_type==1){
            c_const_le = c_const_vk_le;
            d_const_le = d_const_vk_le;
         }
        else {
            c_const_le=0;
            d_const_le = 0;
        }

if (m_parameter->lowFreq) {
Dh[i]=calcDh(Mach[i],theta_e[i],phi_e[i],EddyMach);
Dh_rotor[i]=calcDh(Mach_rotor[i],theta_e_rotor[i],phi_e_rotor[i],EddyMach);
}
else{Dh[i]=1; Dh_rotor[i]=1;}
if (m_parameter->lowFreq) {
Dl[i]=calcDl(Mach[i],theta_e[i],phi_e[i]);
Dl_rotor[i]=calcDl(Mach_rotor[i],theta_e_rotor[i],phi_e_rotor[i]);
}
else{Dl[i]=1; Dl_rotor[i]=1;}

Dl_le[i]=0.5*Dl[i];

Dl_le_rotor[i]=0.5*Dl_rotor[i];

aux0_le[i]=0.5*(lambda_le*L_le*pow(rho/1000., 2)*pow(c_0_le, 2)*pow(u_le, 2)*pow(Mach[i], 3)*pow(I_le, 2)*Dl_le[i])/(pow(r_e_le[i], 2));
aux1_le[i]=10.*log10(pow(LFC_le/(1+LFC_le), 2))+d_const_le;
aux4_le[i]=pow(K_le,3)/pow(1+(pow(K_le,2)),c_const_le);
aux5_le[i]=10.*log10(aux0_le[i]*aux4_le[i]);

//rotor
aux0_le_rotor[i]=0.5*(lambda_le*L_le*pow(rho/1000., 2)*pow(c_0_le, 2)*pow(u_le, 2)*pow(Mach_rotor[i], 3)*pow(I_le, 2)*Dl_le_rotor[i])/(pow(r_e_le_rotor[i], 2));
aux1_le_rotor[i]=10.*log10(pow(LFC_le_rotor/(1+LFC_le_rotor), 2))+d_const_le;
aux4_le_rotor[i]=pow(K_le,3)/pow(1+(pow(K_le,2)),c_const_le);
aux5_le_rotor[i]=10.*log10(aux0_le_rotor[i]*aux4_le_rotor[i]);
//rotor

//Validation:

//BPM validation:
//p 17 C_Project_Log_Text_15_jan_16
//if(((((alpha[i]<=19.8 & Reynolds[i]<3*pow(10,6)) & (Mach[i]<0.21)) & (Reynolds[i]>0)) & (Mach[i]>0))){
//if((alpha[i]<=19.8) & (Mach[i]<0.21) & (Reynolds[i]>0) & (Mach[i]>0)){
BPM_validation=true;
//}
//else{
//SPL_alpha[j]=-999999999999.;
//SPL_S[j]=-999999999999.;
//SPL_P[j]=-999999999999.;
//SPL_dB[j]=-999999999999.;
//SPL_A[j]=-999999999999.;
//SPL_B[j]=-999999999999.;
//SPL_C[j]=-999999999999.;
//BPM_validation=false;
//}

//Lowson validation:
//if (((((m_parameter.Lowson_type!=0) & (Mach[i]<=0.18)) & (Mach[i]>0) & (Reynolds[i]<=(6.*pow(10,5)))) & (Reynolds[i]>0))){
//if (((((m_parameter.Lowson_type!=0) & (Mach[i]<=0.18)) & (Mach[i]>0)))){
SPL_LedB[j]=10.*log10(pow(10,(aux1_le[i]+aux5_le[i])/10.));
SPL_LedBAW[j]=SPL_LedB[j]+AWeighting[j];
SPL_LedBBW[j]=SPL_LedB[j]+BWeighting[j];
SPL_LedBCW[j]=SPL_LedB[j]+CWeighting[j];
//rotor
SPL_LedB_rotor[j]=10.*log10(pow(10,(aux1_le_rotor[i]+aux5_le_rotor[i])/10.));
SPL_LedBAW_rotor[j]=SPL_LedB_rotor[j]+AWeighting[j];
SPL_LedBBW_rotor[j]=SPL_LedB_rotor[j]+BWeighting[j];
SPL_LedBCW_rotor[j]=SPL_LedB_rotor[j]+CWeighting[j];
//rotor
LE_validation=true;

//}
//else{
//SPL_LedB[j]=-999999999999.;
//SPL_LedBAW[j]=-999999999999.;
//SPL_LedBBW[j]=-999999999999.;
//SPL_LedBCW[j]=-999999999999.;
//LE_validation=false;
//}

first_term_Dl_S[i]=calcFirstTerm(Mach[i],L[i],Dl[i],D_starred_S[i],dist_obs[i]);
first_term_Dh_P[i]=calcFirstTerm(Mach[i],L[i],Dh[i],D_starred_P[i],dist_obs[i]);
first_term_Dh_S[i]=calcFirstTerm(Mach[i],L[i],Dh[i],D_starred_S[i],dist_obs[i]);

//rotor
first_term_Dl_S_rotor[i]=calcFirstTerm(Mach_rotor[i],L[i],Dl_rotor[i],D_starred_S_rotor[i],dist_obs_rotor[i]);
first_term_Dh_P_rotor[i]=calcFirstTerm(Mach_rotor[i],L[i],Dh_rotor[i],D_starred_P_rotor[i],dist_obs_rotor[i]);
first_term_Dh_S_rotor[i]=calcFirstTerm(Mach_rotor[i],L[i],Dh_rotor[i],D_starred_S_rotor[i],dist_obs_rotor[i]);
//rotor

SPL_dB_P[j]=delta_K1[i]+A_a_P[j]+K1_3[i]+first_term_Dh_P[i];
SPL_dB_P_rotor[j]=delta_K1_rotor[i]+A_a_P_rotor[j]+K1_3_rotor[i]+first_term_Dh_P_rotor[i];

if (alpha[i]<SwAlpha[i]){SPL_P[j]=SPL_dB_P[j]; SPL_P_rotor[j]=SPL_dB_P_rotor[j];}
else {SPL_P[j]=-999999999999.; SPL_P_rotor[j]=-999999999999.;}

SPL_dB_S[j]=first_term_Dh_S[i]+A_a_S[j]+K1_3[i];
SPL_dB_S_rotor[j]=first_term_Dh_S_rotor[i]+A_a_S_rotor[j]+K1_3_rotor[i];

if (alpha[i]<SwAlpha[i]){SPL_S[j]=SPL_dB_S[j]; SPL_S_rotor[j]=SPL_dB_S_rotor[j];}
else {SPL_S[j]=-999999999999.; SPL_S_rotor[j]=-999999999999.;}

SPL_alpha_min0[j]=first_term_Dh_S[i]+K2[i]+B_b[j];
SPL_alpha_big0[j]=first_term_Dl_S[i]+K2[i]+Alin_a[j];

SPL_alpha_min0_rotor[j]=first_term_Dh_S_rotor[i]+K2_rotor[i]+B_b_rotor[j];
SPL_alpha_big0_rotor[j]=first_term_Dl_S_rotor[i]+K2_rotor[i]+Alin_a_rotor[j];

if (alpha[i]<SwAlpha[i]){SPL_alpha[j]=SPL_alpha_min0[j]; SPL_alpha_rotor[j]=SPL_alpha_min0_rotor[j];}
else {SPL_alpha[j]=SPL_alpha_big0[j]; SPL_alpha_rotor[j]=SPL_alpha_big0_rotor[j];}

if(m_parameter->Lowson_type!=0){SPL_dB[j]=10.*log10(pow(10.,(SPL_alpha[j]/10.))+pow(10.,(SPL_S[j]/10.))+pow(10.,(SPL_P[j]/10.))+pow(10.,(SPL_LedB[j]/10.)));
SPL_dB_rotor[j]=10.*log10(pow(10.,(SPL_alpha_rotor[j]/10.))+pow(10.,(SPL_S_rotor[j]/10.))+pow(10.,(SPL_P_rotor[j]/10.))+pow(10.,(SPL_LedB_rotor[j]/10.)));
}
else{SPL_dB[j]=10.*log10(pow(10.,(SPL_alpha[j]/10.))+pow(10.,(SPL_S[j]/10.))+pow(10.,(SPL_P[j]/10.)));
SPL_dB_rotor[j]=10.*log10(pow(10.,(SPL_alpha_rotor[j]/10.))+pow(10.,(SPL_S_rotor[j]/10.))+pow(10.,(SPL_P_rotor[j]/10.)));
}

SPL_A[j]=SPL_dB[j]+AWeighting[j];
SPL_B[j]=SPL_dB[j]+BWeighting[j];
SPL_C[j]=SPL_dB[j]+CWeighting[j];

SPL_A_rotor[j]=SPL_dB_rotor[j]+AWeighting[j];
SPL_B_rotor[j]=SPL_dB_rotor[j]+BWeighting[j];
SPL_C_rotor[j]=SPL_dB_rotor[j]+CWeighting[j];

if(BPM_validation){
aux_m_SPLadB3d=SPL_alpha[j];
aux_m_SPLsdB3d=SPL_S[j];
aux_m_SPLpdB3d=SPL_P[j];
aux_m_SPLdB3d=SPL_dB[j];
aux_m_SPLdBAW3d=SPL_A[j];
aux_m_SPLdBBW3d=SPL_B[j];
aux_m_SPLdBCW3d=SPL_C[j];

aux_m_SPLadB3d_rotor=SPL_alpha_rotor[j];
aux_m_SPLsdB3d_rotor=SPL_S_rotor[j];
aux_m_SPLpdB3d_rotor=SPL_P_rotor[j];
aux_m_SPLdB3d_rotor=SPL_dB_rotor[j];
aux_m_SPLdBAW3d_rotor=SPL_A_rotor[j];
aux_m_SPLdBBW3d_rotor=SPL_B_rotor[j];
aux_m_SPLdBCW3d_rotor=SPL_C_rotor[j];
}else{
aux_m_SPLadB3d=0;
aux_m_SPLsdB3d=0;
aux_m_SPLpdB3d=0;
aux_m_SPLdB3d=0;
aux_m_SPLdBAW3d=0;
aux_m_SPLdBBW3d=0;
aux_m_SPLdBCW3d=0;

aux_m_SPLadB3d_rotor=0;
aux_m_SPLsdB3d_rotor=0;
aux_m_SPLpdB3d_rotor=0;
aux_m_SPLdB3d_rotor=0;
aux_m_SPLdBAW3d_rotor=0;
aux_m_SPLdBBW3d_rotor=0;
aux_m_SPLdBCW3d_rotor=0;
}
if (LE_validation!=0){
aux_m_SPL_LEdB3d=SPL_LedB[j];
aux_m_SPL_LEdBAW3d=SPL_LedBAW[j];
aux_m_SPL_LEdBBW3d=SPL_LedBBW[j];
aux_m_SPL_LEdBCW3d=SPL_LedBCW[j];

aux_m_SPL_LEdB3d_rotor=SPL_LedB_rotor[j];
aux_m_SPL_LEdBAW3d_rotor=SPL_LedBAW_rotor[j];
aux_m_SPL_LEdBBW3d_rotor=SPL_LedBBW_rotor[j];
aux_m_SPL_LEdBCW3d_rotor=SPL_LedBCW_rotor[j];
    }
    else{
aux_m_SPL_LEdB3d=0;
aux_m_SPL_LEdBAW3d=0;
aux_m_SPL_LEdBBW3d=0;
aux_m_SPL_LEdBCW3d=0;

aux_m_SPL_LEdB3d_rotor=0;
aux_m_SPL_LEdBAW3d_rotor=0;
aux_m_SPL_LEdBBW3d_rotor=0;
aux_m_SPL_LEdBCW3d_rotor=0;
    }

if(qIsInf(aux_m_SPLadB3d) || qIsNaN(aux_m_SPLadB3d)){aux_m_SPLadB3d=0;}
if(qIsInf(aux_m_SPLsdB3d) || qIsNaN(aux_m_SPLsdB3d)){aux_m_SPLsdB3d=0;}
if(qIsInf(aux_m_SPLpdB3d) || qIsNaN(aux_m_SPLpdB3d)){aux_m_SPLpdB3d=0;}
if(qIsInf(aux_m_SPLdB3d) || qIsNaN(aux_m_SPLdB3d)){aux_m_SPLdB3d=0;}
if(qIsInf(aux_m_SPLdBAW3d) || qIsNaN(aux_m_SPLdBAW3d)){aux_m_SPLdBAW3d=0;}
if(qIsInf(aux_m_SPLdBBW3d) || qIsNaN(aux_m_SPLdBBW3d)){aux_m_SPLdBBW3d=0;}
if(qIsInf(aux_m_SPLdBCW3d) || qIsNaN(aux_m_SPLdBCW3d)){aux_m_SPLdBCW3d=0;}
if(qIsInf(aux_m_SPL_LEdB3d) || qIsNaN(aux_m_SPL_LEdB3d)){aux_m_SPL_LEdB3d=0;}
if(qIsInf(aux_m_SPL_LEdBAW3d) || qIsNaN(aux_m_SPL_LEdBAW3d)){aux_m_SPL_LEdBAW3d=0;}
if(qIsInf(aux_m_SPL_LEdBBW3d) || qIsNaN(aux_m_SPL_LEdBBW3d)){aux_m_SPL_LEdBBW3d=0;}
if(qIsInf(aux_m_SPL_LEdBCW3d) || qIsNaN(aux_m_SPL_LEdBCW3d)){aux_m_SPL_LEdBCW3d=0;}

if(qIsInf(aux_m_SPLadB3d_rotor) || qIsNaN(aux_m_SPLadB3d_rotor)){aux_m_SPLadB3d_rotor=0;}
if(qIsInf(aux_m_SPLsdB3d_rotor) || qIsNaN(aux_m_SPLsdB3d_rotor)){aux_m_SPLsdB3d_rotor=0;}
if(qIsInf(aux_m_SPLpdB3d_rotor) || qIsNaN(aux_m_SPLpdB3d_rotor)){aux_m_SPLpdB3d_rotor=0;}
if(qIsInf(aux_m_SPLdB3d_rotor) || qIsNaN(aux_m_SPLdB3d_rotor)){aux_m_SPLdB3d_rotor=0;}
if(qIsInf(aux_m_SPLdBAW3d_rotor) || qIsNaN(aux_m_SPLdBAW3d_rotor)){aux_m_SPLdBAW3d_rotor=0;}
if(qIsInf(aux_m_SPLdBBW3d_rotor) || qIsNaN(aux_m_SPLdBBW3d_rotor)){aux_m_SPLdBBW3d_rotor=0;}
if(qIsInf(aux_m_SPLdBCW3d_rotor) || qIsNaN(aux_m_SPLdBCW3d_rotor)){aux_m_SPLdBCW3d_rotor=0;}
if(qIsInf(aux_m_SPL_LEdB3d_rotor) || qIsNaN(aux_m_SPL_LEdB3d_rotor)){aux_m_SPL_LEdB3d_rotor=0;}
if(qIsInf(aux_m_SPL_LEdBAW3d_rotor) || qIsNaN(aux_m_SPL_LEdBAW3d_rotor)){aux_m_SPL_LEdBAW3d_rotor=0;}
if(qIsInf(aux_m_SPL_LEdBBW3d_rotor) || qIsNaN(aux_m_SPL_LEdBBW3d_rotor)){aux_m_SPL_LEdBBW3d_rotor=0;}
if(qIsInf(aux_m_SPL_LEdBCW3d_rotor) || qIsNaN(aux_m_SPL_LEdBCW3d_rotor)){aux_m_SPL_LEdBCW3d_rotor=0;}

//multi 3D curves
if((blade==0) & (E==0)){
m_SPLadB3d[i][j]=aux_m_SPLadB3d;
m_SPLsdB3d[i][j]=aux_m_SPLsdB3d;
m_SPLpdB3d[i][j]=aux_m_SPLpdB3d;
m_SPLdB3d[i][j]=aux_m_SPLdB3d;
m_SPLdBAW3d[i][j]=aux_m_SPLdBAW3d;
m_SPLdBBW3d[i][j]=aux_m_SPLdBBW3d;
m_SPLdBCW3d[i][j]=aux_m_SPLdBCW3d;

if(m_parameter->Lowson_type!=0){
m_SPL_LEdB3d[i][j]=aux_m_SPL_LEdB3d;
m_SPL_LEdBAW3d[i][j]=aux_m_SPL_LEdBAW3d;
m_SPL_LEdBBW3d[i][j]=aux_m_SPL_LEdBBW3d;
m_SPL_LEdBCW3d[i][j]=aux_m_SPL_LEdBCW3d;
}
else{
m_SPL_LEdB3d[i][j]=0;
m_SPL_LEdBAW3d[i][j]=0;
m_SPL_LEdBBW3d[i][j]=0;
m_SPL_LEdBCW3d[i][j]=0;
}}

m_SPLadB3d_4d_blade[i][j][blade][E]=aux_m_SPLadB3d;
m_SPLsdB3d_4d_blade[i][j][blade][E]=aux_m_SPLsdB3d;
m_SPLpdB3d_4d_blade[i][j][blade][E]=aux_m_SPLpdB3d;
m_SPLdB3d_4d_blade[i][j][blade][E]=aux_m_SPLdB3d;
m_SPLdBAW3d_4d_blade[i][j][blade][E]=aux_m_SPLdBAW3d;
m_SPLdBBW3d_4d_blade[i][j][blade][E]=aux_m_SPLdBBW3d;
m_SPLdBCW3d_4d_blade[i][j][blade][E]=aux_m_SPLdBCW3d;

m_SPLadB3d_4d[i][j][blade][E]=aux_m_SPLadB3d_rotor;
m_SPLsdB3d_4d[i][j][blade][E]=aux_m_SPLsdB3d_rotor;
m_SPLpdB3d_4d[i][j][blade][E]=aux_m_SPLpdB3d_rotor;
m_SPLdB3d_4d[i][j][blade][E]=aux_m_SPLdB3d_rotor;
m_SPLdBAW3d_4d[i][j][blade][E]=aux_m_SPLdBAW3d_rotor;
m_SPLdBBW3d_4d[i][j][blade][E]=aux_m_SPLdBBW3d_rotor;
m_SPLdBCW3d_4d[i][j][blade][E]=aux_m_SPLdBCW3d_rotor;

if(m_parameter->Lowson_type!=0){
m_SPL_LEdB3d_4d_blade[i][j][blade][E]=aux_m_SPL_LEdB3d;
m_SPL_LEdBAW3d_4d_blade[i][j][blade][E]=aux_m_SPL_LEdBAW3d;
m_SPL_LEdBBW3d_4d_blade[i][j][blade][E]=aux_m_SPL_LEdBBW3d;
m_SPL_LEdBCW3d_4d_blade[i][j][blade][E]=aux_m_SPL_LEdBCW3d;

m_SPL_LEdB3d_4d[i][j][blade][E]=aux_m_SPL_LEdB3d_rotor;
m_SPL_LEdBAW3d_4d[i][j][blade][E]=aux_m_SPL_LEdBAW3d_rotor;
m_SPL_LEdBBW3d_4d[i][j][blade][E]=aux_m_SPL_LEdBBW3d_rotor;
m_SPL_LEdBCW3d_4d[i][j][blade][E]=aux_m_SPL_LEdBCW3d_rotor;
}
else{
m_SPL_LEdB3d_4d[i][j][blade][E]=0;
m_SPL_LEdBAW3d_4d[i][j][blade][E]=0;
m_SPL_LEdBBW3d_4d[i][j][blade][E]=0;
m_SPL_LEdBCW3d_4d[i][j][blade][E]=0;

m_SPL_LEdB3d_4d_blade[i][j][blade][E]=0;
m_SPL_LEdBAW3d_4d_blade[i][j][blade][E]=0;
m_SPL_LEdBBW3d_4d_blade[i][j][blade][E]=0;
m_SPL_LEdBCW3d_4d_blade[i][j][blade][E]=0;
}
}}
    }
z=z+ldelta;
}}

//calculation for rotor in loop
void NoiseCalculation::calculateqs3d_graphics_loops(){
QList<NoiseOpPoint*> noiseOpPoints = m_parameter->prepareNoiseOpPointList();
setupVectorsqs3d();

//D star interpolated
SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
    double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
    double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
    double z=lstart;

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    int number_of_segments = pBEM->m_pBData->m_pos.size();
    double chord[number_of_segments];
    m_DStarInterpolatedS3d.resize(number_of_segments+1);
    m_DStarInterpolatedP3d.resize(number_of_segments+1);
    double chordmax=0;
    double TSR = m_parameter->TSRtd;

for (int posOpPoint = 0; posOpPoint < noiseOpPoints.size(); ++posOpPoint) {
    NoiseOpPoint *nop = noiseOpPoints[posOpPoint];

foreach(BData * bdata, pBEM->m_pBEMData->GetBData()){
    if (z==TSR){

for (int j = 0; j < number_of_segments; ++j) {
chord[j] = bdata->m_c_local.value(j);
if(chord[j]>chordmax){chordmax=chord[j];}

bool dStarOrder = false;

//When angle is negative D* search must be inverted
if(nop->getAlphaDegree() < 0){
    dStarOrder = true;
}

if (m_parameter->opPointSource == NoiseParameter::OnePolar ||
    m_parameter->opPointSource == NoiseParameter::MultiplePolars)
{

m_DStarInterpolatedS3d[j] = getDStarInterpolated3d(dStarOrder,(chord[j]/(chordmax)),nop);
//qDebug() << "D* S 3d:" << m_DStarInterpolatedS3d[j];
m_DStarInterpolatedP3d[j] = getDStarInterpolated3d(!dStarOrder,(chord[j]/(chordmax)),nop);
//qDebug() << "D* P 3d:" << m_DStarInterpolatedP3d;
}}}
z=z+ldelta;
}}
//finish D star interpolated

QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
int blades_num = pbem->m_pBData->blades;
double anglesteps;
int number_of_rotations;

if (m_parameter->rotation_type==0){
//    angle based
number_of_rotations = m_parameter->number_loops;
anglesteps=m_parameter->anglesteps;
}else{
//    time based
number_of_rotations = m_parameter->time/(60./m_parameter->rot_speed);
anglesteps=m_parameter->timesteps*60.*360./(m_parameter->rot_speed*1000.);
}

int angles_num=360./anglesteps*number_of_rotations;

for (int blade=0;blade<blades_num;++blade){
for (int E=0;E<angles_num;++E){
        calculateqs3d_graphics(blade,E,TSR);
}}

ProgressBar(2);//Sara
}

//calculation for blade
void NoiseCalculation::calculateqs3d_blade() {
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter->prepareNoiseOpPointList();
    QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
    unsigned int number_of_segments = pbem->dlg_elements;

        double aux_m_SPLadB3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPLsdB3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPLpdB3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPLdB3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPLdBAW3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPLdBBW3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPLdBCW3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPL_LEdB3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPL_LEdBAW3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPL_LEdBBW3d_final[FREQUENCY_TABLE_SIZE];
        double aux_m_SPL_LEdBCW3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPLadB3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPLsdB3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPLpdB3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPLdB3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPLdBAW3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPLdBBW3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPLdBCW3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPL_LEdB3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPL_LEdBAW3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPL_LEdBBW3d_final[FREQUENCY_TABLE_SIZE];
        double auxa_m_SPL_LEdBCW3d_final[FREQUENCY_TABLE_SIZE];

        double Final_qs3d_alpha_aux=0;
        double Final_qs3d_S_aux=0;
        double Final_qs3d_P_aux=0;
        double Final_qs3d_LE_aux=0;
        double Final_qs3d_aux=0;

        // Resize vectors acording to OpPoints total
        unsigned int sizea = 0;
        unsigned int size;
        if (m_parameter->opPointSource == NoiseParameter::OnePolar ||
            m_parameter->opPointSource == NoiseParameter::MultiplePolars)
        {
            sizea = m_parameter->analyzedOpPoints.size();
         } else if (m_parameter->opPointSource == NoiseParameter::OriginalBpm) {
            sizea = 1;
        }

    if (number_of_segments>sizea){size = number_of_segments;} else {size = sizea;}
        double auxa_m_OASPL3d[size];
        double auxa_m_OASPLA3d[size];
        double auxa_m_OASPLB3d[size];
        double auxa_m_OASPLC3d[size];
        double auxa_m_SPLALOG3d[size];
        double auxa_m_SPLSLOG3d[size];
        double auxa_m_SPLPLOG3d[size];
        double auxa_m_SPLLEdBAW3d[size];
        double auxa_m_SPLLEdBBW3d[size];
        double auxa_m_SPLLEdBCW3d[size];
        double auxa_m_SPLlogLE3d[size];

    for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
        aux_m_SPLadB3d_final[j]=0;
        aux_m_SPLsdB3d_final[j]=0;
        aux_m_SPLpdB3d_final[j]=0;
        aux_m_SPLdB3d_final[j]=0;
        aux_m_SPLdBAW3d_final[j]=0;
        aux_m_SPLdBBW3d_final[j]=0;
        aux_m_SPLdBCW3d_final[j]=0;
        aux_m_SPL_LEdB3d_final[j]=0;
        aux_m_SPL_LEdBAW3d_final[j]=0;
        aux_m_SPL_LEdBBW3d_final[j]=0;
        aux_m_SPL_LEdBCW3d_final[j]=0;
        auxa_m_SPLadB3d_final[j]=0;
        auxa_m_SPLsdB3d_final[j]=0;
        auxa_m_SPLpdB3d_final[j]=0;
        auxa_m_SPLdB3d_final[j]=0;
        auxa_m_SPLdBAW3d_final[j]=0;
        auxa_m_SPLdBBW3d_final[j]=0;
        auxa_m_SPLdBCW3d_final[j]=0;
        auxa_m_SPL_LEdB3d_final[j]=0;
        auxa_m_SPL_LEdBAW3d_final[j]=0;
        auxa_m_SPL_LEdBBW3d_final[j]=0;
        auxa_m_SPL_LEdBCW3d_final[j]=0;

    for (unsigned int i = 0; i < number_of_segments; ++i) {
        if (SPLadB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLadB3d_final[j] += 0;} else {auxa_m_SPLadB3d_final[j] += pow(10.,(SPLadB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLsdB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLsdB3d_final[j] += 0;} else {auxa_m_SPLsdB3d_final[j] += pow(10.,(SPLsdB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLpdB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLpdB3d_final[j] += 0;} else {auxa_m_SPLpdB3d_final[j] += pow(10.,(SPLpdB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLdB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLdB3d_final[j] += 0;} else {auxa_m_SPLdB3d_final[j] += pow(10.,(SPLdB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLdBAW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLdBAW3d_final[j] += 0;} else {auxa_m_SPLdBAW3d_final[j] += pow(10.,(SPLdBAW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLdBBW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLdBBW3d_final[j] += 0;} else {auxa_m_SPLdBBW3d_final[j] += pow(10.,(SPLdBBW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLdBCW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLdBCW3d_final[j] += 0;} else {auxa_m_SPLdBCW3d_final[j] += pow(10.,(SPLdBCW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPL_LEdB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPL_LEdB3d_final[j] += 0;} else {auxa_m_SPL_LEdB3d_final[j] += pow(10.,(SPL_LEdB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPL_LEdBAW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPL_LEdBAW3d_final[j] += 0;} else {auxa_m_SPL_LEdBAW3d_final[j] += pow(10.,(SPL_LEdBAW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPL_LEdBBW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPL_LEdBBW3d_final[j] += 0;} else {auxa_m_SPL_LEdBBW3d_final[j] += pow(10.,(SPL_LEdBBW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPL_LEdBCW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPL_LEdBCW3d_final[j] += 0;} else {auxa_m_SPL_LEdBCW3d_final[j] += pow(10.,(SPL_LEdBCW3d_4d_blade()[i][j][0][0]/10.));}
    }

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    unsigned int number_of_segments = pBEM->dlg_elements;

    if (number_of_segments>sizea){size = number_of_segments;} else {size = sizea;}

    for (unsigned int i=0;i<size;++i){
    for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
    //    3d curves
        m_SPLadB3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLadB3d_final[j]);
        m_SPLsdB3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLsdB3d_final[j]);
        m_SPLpdB3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLpdB3d_final[j]);
        m_SPLdB3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdB3d_final[j]);
        m_SPLdBAW3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBAW3d_final[j]);
        m_SPLdBBW3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBBW3d_final[j]);
        m_SPLdBCW3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBCW3d_final[j]);
        if (m_parameter->Lowson_type!=0){
        m_SPL_LEdB3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdB3d_final[j]);
        m_SPL_LEdBAW3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBAW3d_final[j]);
        m_SPL_LEdBBW3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBBW3d_final[j]);
        m_SPL_LEdBCW3d_final[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBCW3d_final[j]);
        }
        else{
        m_SPL_LEdB3d_final[i][j]=0;
        m_SPL_LEdBAW3d_final[i][j]=0;
        m_SPL_LEdBBW3d_final[i][j]=0;
        m_SPL_LEdBCW3d_final[i][j]=0;
        }
    }}}

    //OASPL complete for quasi 3d
    int i=number_of_segments-1;
    for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
    Final_qs3d_alpha_aux += pow(10.,(m_SPLadB3d_final[i][j])/10.);
    Final_qs3d_S_aux += pow(10.,(m_SPLsdB3d_final[i][j])/10.);
    Final_qs3d_P_aux += pow(10.,(m_SPLpdB3d_final[i][j])/10.);
    Final_qs3d_LE_aux += pow(10.,(m_SPL_LEdB3d_final[i][j])/10.);
    Final_qs3d_aux += pow(10.,(m_SPLdB3d_final[i][j])/10.);
        }

    Final_qs3d_alpha = 10*log10(Final_qs3d_alpha_aux);
    Final_qs3d_S = 10*log10(Final_qs3d_S_aux);
    Final_qs3d_P = 10*log10(Final_qs3d_P_aux);
    Final_qs3d_LE =  10*log10(Final_qs3d_LE_aux);
    Final_qs3d = 10*log10(Final_qs3d_aux);

    //calculation for the OASPL for the csv output file
    for (unsigned int i=0;i<size;++i){
        auxa_m_OASPL3d[i]=0;
        auxa_m_OASPLA3d[i]=0;
        auxa_m_OASPLB3d[i]=0;
        auxa_m_OASPLC3d[i]=0;
        auxa_m_SPLALOG3d[i]=0;
        auxa_m_SPLSLOG3d[i]=0;
        auxa_m_SPLPLOG3d[i]=0;
        auxa_m_SPLLEdBAW3d[i]=0;
        auxa_m_SPLLEdBBW3d[i]=0;
        auxa_m_SPLLEdBCW3d[i]=0;
        auxa_m_SPLlogLE3d[i]=0;

    for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
        if (SPLadB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLALOG3d[i] += 0;} else {auxa_m_SPLALOG3d[i] += pow(10.,(SPLadB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLsdB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLSLOG3d[i] += 0;} else {auxa_m_SPLSLOG3d[i] += pow(10.,(SPLsdB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLpdB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLPLOG3d[i] += 0;} else {auxa_m_SPLPLOG3d[i] += pow(10.,(SPLpdB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLdB3d_4d_blade()[i][j][0][0]==0.){auxa_m_OASPL3d[i] += 0;} else {auxa_m_OASPL3d[i]+= pow(10.,(SPLdB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLdBAW3d_4d_blade()[i][j][0][0]==0.){auxa_m_OASPLA3d[i] += 0;} else {auxa_m_OASPLA3d[i] += pow(10.,(SPLdBAW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLdBBW3d_4d_blade()[i][j][0][0]==0.){auxa_m_OASPLB3d[i] += 0;} else {auxa_m_OASPLB3d[i] += pow(10.,(SPLdBBW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPLdBCW3d_4d_blade()[i][j][0][0]==0.){auxa_m_OASPLC3d[i] += 0;} else {auxa_m_OASPLC3d[i] += pow(10.,(SPLdBCW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPL_LEdB3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLlogLE3d[i] += 0;} else {auxa_m_SPLlogLE3d[i] += pow(10.,(SPL_LEdB3d_4d_blade()[i][j][0][0]/10.));}
        if (SPL_LEdBAW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLLEdBAW3d[i] += 0;} else {auxa_m_SPLLEdBAW3d[i] += pow(10.,(SPL_LEdBAW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPL_LEdBBW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLLEdBBW3d[i] += 0;} else {auxa_m_SPLLEdBBW3d[i] += pow(10.,(SPL_LEdBBW3d_4d_blade()[i][j][0][0]/10.));}
        if (SPL_LEdBCW3d_4d_blade()[i][j][0][0]==0.){auxa_m_SPLLEdBCW3d[i] += 0;} else {auxa_m_SPLLEdBCW3d[i] += pow(10.,(SPL_LEdBCW3d_4d_blade()[i][j][0][0]/10.));}
    }
    m_OASPL3d[i]=10*log10(auxa_m_OASPL3d[i]);
    m_OASPLA3d[i]=10*log10(auxa_m_OASPLA3d[i]);
    m_OASPLB3d[i]=10*log10(auxa_m_OASPLB3d[i]);
    m_OASPLC3d[i]=10*log10(auxa_m_OASPLC3d[i]);
    m_SPLALOG3d[i]=10*log10(auxa_m_SPLALOG3d[i]);
    m_SPLSLOG3d[i]=10*log10(auxa_m_SPLSLOG3d[i]);
    m_SPLPLOG3d[i]=10*log10(auxa_m_SPLPLOG3d[i]);
    if (m_parameter->Lowson_type!=0){
    m_SPLLEdBAW3d[i]=10*log10(auxa_m_SPLLEdBAW3d[i]);
    m_SPLLEdBBW3d[i]=10*log10(auxa_m_SPLLEdBBW3d[i]);
    m_SPLLEdBCW3d[i]=10*log10(auxa_m_SPLLEdBCW3d[i]);
    m_SPLlogLE3d[i]=10*log10(auxa_m_SPLlogLE3d[i]);
    }
    else{
    m_SPLLEdBAW3d[i]=0;
    m_SPLLEdBBW3d[i]=0;
    m_SPLLEdBCW3d[i]=0;
    m_SPLlogLE3d[i]=0;
    }}

ProgressBar(3);//Sara
}

//calculation for one blade placed in rotor coordinates
void NoiseCalculation::calculateqs3d_rotor() {

QList<NoiseOpPoint*> noiseOpPoints = m_parameter->prepareNoiseOpPointList();
QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
unsigned int number_of_segments = pbem->dlg_elements;

    double aux_m_SPLadB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLsdB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLpdB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLdB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLdBAW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLdBBW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLdBCW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPL_LEdB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPL_LEdBAW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPL_LEdBBW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double aux_m_SPL_LEdBCW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLadB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLsdB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLpdB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLdB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLdBAW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLdBBW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLdBCW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPL_LEdB3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPL_LEdBAW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPL_LEdBBW3d_final_rotor[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPL_LEdBCW3d_final_rotor[FREQUENCY_TABLE_SIZE];

    double Final_qs3d_alpha_aux_rotor=0;
    double Final_qs3d_S_aux_rotor=0;
    double Final_qs3d_P_aux_rotor=0;
    double Final_qs3d_LE_aux_rotor=0;
    double Final_qs3d_aux_rotor=0;

    // Resize vectors acording to OpPoints total
    unsigned int sizea = 0;
    unsigned int size;
    if (m_parameter->opPointSource == NoiseParameter::OnePolar ||
        m_parameter->opPointSource == NoiseParameter::MultiplePolars)
    {
        sizea = m_parameter->analyzedOpPoints.size();
     } else if (m_parameter->opPointSource == NoiseParameter::OriginalBpm) {
        sizea = 1;
    }

if (number_of_segments>sizea){size = number_of_segments;} else {size = sizea;}
    double auxa_m_OASPL3d_rotor[size];
    double auxa_m_OASPLA3d_rotor[size];
    double auxa_m_OASPLB3d_rotor[size];
    double auxa_m_OASPLC3d_rotor[size];
    double auxa_m_SPLALOG3d_rotor[size];
    double auxa_m_SPLSLOG3d_rotor[size];
    double auxa_m_SPLPLOG3d_rotor[size];
    double auxa_m_SPLLEdBAW3d_rotor[size];
    double auxa_m_SPLLEdBBW3d_rotor[size];
    double auxa_m_SPLLEdBCW3d_rotor[size];
    double auxa_m_SPLlogLE3d_rotor[size];

for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
    aux_m_SPLadB3d_final_rotor[j]=0;
    aux_m_SPLsdB3d_final_rotor[j]=0;
    aux_m_SPLpdB3d_final_rotor[j]=0;
    aux_m_SPLdB3d_final_rotor[j]=0;
    aux_m_SPLdBAW3d_final_rotor[j]=0;
    aux_m_SPLdBBW3d_final_rotor[j]=0;
    aux_m_SPLdBCW3d_final_rotor[j]=0;
    aux_m_SPL_LEdB3d_final_rotor[j]=0;
    aux_m_SPL_LEdBAW3d_final_rotor[j]=0;
    aux_m_SPL_LEdBBW3d_final_rotor[j]=0;
    aux_m_SPL_LEdBCW3d_final_rotor[j]=0;
    auxa_m_SPLadB3d_final_rotor[j]=0;
    auxa_m_SPLsdB3d_final_rotor[j]=0;
    auxa_m_SPLpdB3d_final_rotor[j]=0;
    auxa_m_SPLdB3d_final_rotor[j]=0;
    auxa_m_SPLdBAW3d_final_rotor[j]=0;
    auxa_m_SPLdBBW3d_final_rotor[j]=0;
    auxa_m_SPLdBCW3d_final_rotor[j]=0;
    auxa_m_SPL_LEdB3d_final_rotor[j]=0;
    auxa_m_SPL_LEdBAW3d_final_rotor[j]=0;
    auxa_m_SPL_LEdBBW3d_final_rotor[j]=0;
    auxa_m_SPL_LEdBCW3d_final_rotor[j]=0;

    for (unsigned int i = 0; i < number_of_segments; ++i) {
        if (SPLadB3d_4d()[i][j][0][0]==0.){auxa_m_SPLadB3d_final_rotor[j] += 0;} else {auxa_m_SPLadB3d_final_rotor[j] += pow(10.,(SPLadB3d_4d()[i][j][0][0]/10.));}
        if (SPLsdB3d_4d()[i][j][0][0]==0.){auxa_m_SPLsdB3d_final_rotor[j] += 0;} else {auxa_m_SPLsdB3d_final_rotor[j] += pow(10.,(SPLsdB3d_4d()[i][j][0][0]/10.));}
        if (SPLpdB3d_4d()[i][j][0][0]==0.){auxa_m_SPLpdB3d_final_rotor[j] += 0;} else {auxa_m_SPLpdB3d_final_rotor[j] += pow(10.,(SPLpdB3d_4d()[i][j][0][0]/10.));}
        if (SPLdB3d_4d()[i][j][0][0]==0.){auxa_m_SPLdB3d_final_rotor[j] += 0;} else {auxa_m_SPLdB3d_final_rotor[j] += pow(10.,(SPLdB3d_4d()[i][j][0][0]/10.));}
        if (SPLdBAW3d_4d()[i][j][0][0]==0.){auxa_m_SPLdBAW3d_final_rotor[j] += 0;} else {auxa_m_SPLdBAW3d_final_rotor[j] += pow(10.,(SPLdBAW3d_4d()[i][j][0][0]/10.));}
        if (SPLdBBW3d_4d()[i][j][0][0]==0.){auxa_m_SPLdBBW3d_final_rotor[j] += 0;} else {auxa_m_SPLdBBW3d_final_rotor[j] += pow(10.,(SPLdBBW3d_4d()[i][j][0][0]/10.));}
        if (SPLdBCW3d_4d()[i][j][0][0]==0.){auxa_m_SPLdBCW3d_final_rotor[j] += 0;} else {auxa_m_SPLdBCW3d_final_rotor[j] += pow(10.,(SPLdBCW3d_4d()[i][j][0][0]/10.));}
        if (SPL_LEdB3d_4d()[i][j][0][0]==0.){auxa_m_SPL_LEdB3d_final_rotor[j] += 0;} else {auxa_m_SPL_LEdB3d_final_rotor[j] += pow(10.,(SPL_LEdB3d_4d()[i][j][0][0]/10.));}
        if (SPL_LEdBAW3d_4d()[i][j][0][0]==0.){auxa_m_SPL_LEdBAW3d_final_rotor[j] += 0;} else {auxa_m_SPL_LEdBAW3d_final_rotor[j] += pow(10.,(SPL_LEdBAW3d_4d()[i][j][0][0]/10.));}
        if (SPL_LEdBBW3d_4d()[i][j][0][0]==0.){auxa_m_SPL_LEdBBW3d_final_rotor[j] += 0;} else {auxa_m_SPL_LEdBBW3d_final_rotor[j] += pow(10.,(SPL_LEdBBW3d_4d()[i][j][0][0]/10.));}
        if (SPL_LEdBCW3d_4d()[i][j][0][0]==0.){auxa_m_SPL_LEdBCW3d_final_rotor[j] += 0;} else {auxa_m_SPL_LEdBCW3d_final_rotor[j] += pow(10.,(SPL_LEdBCW3d_4d()[i][j][0][0]/10.));}
    }

    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    unsigned int number_of_segments = pBEM->dlg_elements;

    if (number_of_segments>sizea){size = number_of_segments;} else {size = sizea;}

    for (unsigned int i=0;i<size;++i){
    for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
    //    3d curves
        m_SPLadB3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLadB3d_final_rotor[j]);
        m_SPLsdB3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLsdB3d_final_rotor[j]);
        m_SPLpdB3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLpdB3d_final_rotor[j]);
        m_SPLdB3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdB3d_final_rotor[j]);
        m_SPLdBAW3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBAW3d_final_rotor[j]);
        m_SPLdBBW3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBBW3d_final_rotor[j]);
        m_SPLdBCW3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBCW3d_final_rotor[j]);
        if (m_parameter->Lowson_type!=0){
        m_SPL_LEdB3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdB3d_final_rotor[j]);
        m_SPL_LEdBAW3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBAW3d_final_rotor[j]);
        m_SPL_LEdBBW3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBBW3d_final_rotor[j]);
        m_SPL_LEdBCW3d_final_rotor[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBCW3d_final_rotor[j]);
        }
        else{
        m_SPL_LEdB3d_final_rotor[i][j]=0;
        m_SPL_LEdBAW3d_final_rotor[i][j]=0;
        m_SPL_LEdBBW3d_final_rotor[i][j]=0;
        m_SPL_LEdBCW3d_final_rotor[i][j]=0;
        }
    }
    }}

    //OASPL complete for quasi 3d
    int i=number_of_segments-1;
    for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
    Final_qs3d_alpha_aux_rotor += pow(10.,(m_SPLadB3d_final_rotor[i][j])/10.);
    Final_qs3d_S_aux_rotor += pow(10.,(m_SPLsdB3d_final_rotor[i][j])/10.);
    Final_qs3d_P_aux_rotor += pow(10.,(m_SPLpdB3d_final_rotor[i][j])/10.);
    Final_qs3d_LE_aux_rotor += pow(10.,(m_SPL_LEdB3d_final_rotor[i][j])/10.);
    Final_qs3d_aux_rotor += pow(10.,(m_SPLdB3d_final_rotor[i][j])/10.);
        }

    Final_qs3d_alpha_rotor = 10*log10(Final_qs3d_alpha_aux_rotor);
    Final_qs3d_S_rotor = 10*log10(Final_qs3d_S_aux_rotor);
    Final_qs3d_P_rotor = 10*log10(Final_qs3d_P_aux_rotor);
    Final_qs3d_LE_rotor =  10*log10(Final_qs3d_LE_aux_rotor);
    Final_qs3d_rotor = 10*log10(Final_qs3d_aux_rotor);

    //calculation for the OASPL for the csv output file
    for (unsigned int i=0;i<size;++i){
        auxa_m_OASPL3d_rotor[i]=0;
        auxa_m_OASPLA3d_rotor[i]=0;
        auxa_m_OASPLB3d_rotor[i]=0;
        auxa_m_OASPLC3d_rotor[i]=0;
        auxa_m_SPLALOG3d_rotor[i]=0;
        auxa_m_SPLSLOG3d_rotor[i]=0;
        auxa_m_SPLPLOG3d_rotor[i]=0;
        auxa_m_SPLLEdBAW3d_rotor[i]=0;
        auxa_m_SPLLEdBBW3d_rotor[i]=0;
        auxa_m_SPLLEdBCW3d_rotor[i]=0;
        auxa_m_SPLlogLE3d_rotor[i]=0;

    for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
        if (SPLadB3d_4d()[i][j][0][0]==0.){auxa_m_SPLALOG3d_rotor[i] += 0;} else {auxa_m_SPLALOG3d_rotor[i] += pow(10.,(SPLadB3d_4d()[i][j][0][0]/10.));}
        if (SPLsdB3d_4d()[i][j][0][0]==0.){auxa_m_SPLSLOG3d_rotor[i] += 0;} else {auxa_m_SPLSLOG3d_rotor[i] += pow(10.,(SPLsdB3d_4d()[i][j][0][0]/10.));}
        if (SPLpdB3d_4d()[i][j][0][0]==0.){auxa_m_SPLPLOG3d_rotor[i] += 0;} else {auxa_m_SPLPLOG3d_rotor[i] += pow(10.,(SPLpdB3d_4d()[i][j][0][0]/10.));}
        if (SPLdB3d_4d()[i][j][0][0]==0.){auxa_m_OASPL3d_rotor[i] += 0;} else {auxa_m_OASPL3d_rotor[i]+= pow(10.,(SPLdB3d_4d()[i][j][0][0]/10.));}
        if (SPLdBAW3d_4d()[i][j][0][0]==0.){auxa_m_OASPLA3d_rotor[i] += 0;} else {auxa_m_OASPLA3d_rotor[i] += pow(10.,(SPLdBAW3d_4d()[i][j][0][0]/10.));}
        if (SPLdBBW3d_4d()[i][j][0][0]==0.){auxa_m_OASPLB3d_rotor[i] += 0;} else {auxa_m_OASPLB3d_rotor[i] += pow(10.,(SPLdBBW3d_4d()[i][j][0][0]/10.));}
        if (SPLdBCW3d_4d()[i][j][0][0]==0.){auxa_m_OASPLC3d_rotor[i] += 0;} else {auxa_m_OASPLC3d_rotor[i] += pow(10.,(SPLdBCW3d_4d()[i][j][0][0]/10.));}
        if (SPL_LEdB3d_4d()[i][j][0][0]==0.){auxa_m_SPLlogLE3d_rotor[i] += 0;} else {auxa_m_SPLlogLE3d_rotor[i] += pow(10.,(SPL_LEdB3d_4d()[i][j][0][0]/10.));}
        if (SPL_LEdBAW3d_4d()[i][j][0][0]==0.){auxa_m_SPLLEdBAW3d_rotor[i] += 0;} else {auxa_m_SPLLEdBAW3d_rotor[i] += pow(10.,(SPL_LEdBAW3d_4d()[i][j][0][0]/10.));}
        if (SPL_LEdBBW3d_4d()[i][j][0][0]==0.){auxa_m_SPLLEdBBW3d_rotor[i] += 0;} else {auxa_m_SPLLEdBBW3d_rotor[i] += pow(10.,(SPL_LEdBBW3d_4d()[i][j][0][0]/10.));}
        if (SPL_LEdBCW3d_4d()[i][j][0][0]==0.){auxa_m_SPLLEdBCW3d_rotor[i] += 0;} else {auxa_m_SPLLEdBCW3d_rotor[i] += pow(10.,(SPL_LEdBCW3d_4d()[i][j][0][0]/10.));}
    }
    m_OASPL3d_rotor[i]=10*log10(auxa_m_OASPL3d_rotor[i]);
    m_OASPLA3d_rotor[i]=10*log10(auxa_m_OASPLA3d_rotor[i]);
    m_OASPLB3d_rotor[i]=10*log10(auxa_m_OASPLB3d_rotor[i]);
    m_OASPLC3d_rotor[i]=10*log10(auxa_m_OASPLC3d_rotor[i]);
    m_SPLALOG3d_rotor[i]=10*log10(auxa_m_SPLALOG3d_rotor[i]);
    m_SPLSLOG3d_rotor[i]=10*log10(auxa_m_SPLSLOG3d_rotor[i]);
    m_SPLPLOG3d_rotor[i]=10*log10(auxa_m_SPLPLOG3d_rotor[i]);
    if (m_parameter->Lowson_type!=0){
    m_SPLLEdBAW3d_rotor[i]=10*log10(auxa_m_SPLLEdBAW3d_rotor[i]);
    m_SPLLEdBBW3d_rotor[i]=10*log10(auxa_m_SPLLEdBBW3d_rotor[i]);
    m_SPLLEdBCW3d_rotor[i]=10*log10(auxa_m_SPLLEdBCW3d_rotor[i]);
    m_SPLlogLE3d_rotor[i]=10*log10(auxa_m_SPLlogLE3d_rotor[i]);
    }
    else{
    m_SPLLEdBAW3d_rotor[i]=0;
    m_SPLLEdBBW3d_rotor[i]=0;
    m_SPLLEdBCW3d_rotor[i]=0;
    m_SPLlogLE3d_rotor[i]=0;
    }}

ProgressBar(4);//Sara
}

//calculation for all blades in rotation movement
void NoiseCalculation::calculateqs3d_rotor_loops() {
    QList<NoiseOpPoint*> noiseOpPoints = m_parameter->prepareNoiseOpPointList();
    QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
    int blades_num = pbem->m_pBData->blades;
    double anglesteps;
    int number_of_rotations;

    if (m_parameter->rotation_type==0){
    //    angle based
    number_of_rotations = m_parameter->number_loops;
    anglesteps=m_parameter->anglesteps;
    }else{
    //    time based
    number_of_rotations = m_parameter->time/(60./m_parameter->rot_speed);
    anglesteps=m_parameter->timesteps*60.*360./(m_parameter->rot_speed*1000.);
    }

    int angles_num=360./anglesteps*number_of_rotations;

    unsigned int number_of_segments = pbem->dlg_elements;

    double aux_m_SPLadB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLsdB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLpdB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLdB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLdBAW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLdBBW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPLdBCW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPL_LEdB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPL_LEdBAW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPL_LEdBBW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double aux_m_SPL_LEdBCW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLadB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLsdB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLpdB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLdB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLdBAW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLdBBW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPLdBCW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPL_LEdB3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPL_LEdBAW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPL_LEdBBW3d_final_4d[FREQUENCY_TABLE_SIZE];
    double auxa_m_SPL_LEdBCW3d_final_4d[FREQUENCY_TABLE_SIZE];

    double Final_qs3d_alpha_aux_4d=0;
    double Final_qs3d_S_aux_4d=0;
    double Final_qs3d_P_aux_4d=0;
    double Final_qs3d_LE_aux_4d=0;
    double Final_qs3d_aux_4d=0;

        // Resize vectors acording to OpPoints total
        unsigned int sizea = 0;
        unsigned int size;
        if (m_parameter->opPointSource == NoiseParameter::OnePolar ||
            m_parameter->opPointSource == NoiseParameter::MultiplePolars)
        {
            sizea = m_parameter->analyzedOpPoints.size();
         } else if (m_parameter->opPointSource == NoiseParameter::OriginalBpm) {
            sizea = 1;
        }

    if (number_of_segments>sizea){size = number_of_segments;} else {size = sizea;}

    double auxa_m_OASPL3d_rotor_loops[size];
    double auxa_m_OASPLA3d_rotor_loops[size];
    double auxa_m_OASPLB3d_rotor_loops[size];
    double auxa_m_OASPLC3d_rotor_loops[size];
    double auxa_m_SPLALOG3d_rotor_loops[size];
    double auxa_m_SPLSLOG3d_rotor_loops[size];
    double auxa_m_SPLPLOG3d_rotor_loops[size];
    double auxa_m_SPLLEdBAW3d_rotor_loops[size];
    double auxa_m_SPLLEdBBW3d_rotor_loops[size];
    double auxa_m_SPLLEdBCW3d_rotor_loops[size];
    double auxa_m_SPLlogLE3d_rotor_loops[size];

    for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
        aux_m_SPLadB3d_final_4d[j]=0;
        aux_m_SPLsdB3d_final_4d[j]=0;
        aux_m_SPLpdB3d_final_4d[j]=0;
        aux_m_SPLdB3d_final_4d[j]=0;
        aux_m_SPLdBAW3d_final_4d[j]=0;
        aux_m_SPLdBBW3d_final_4d[j]=0;
        aux_m_SPLdBCW3d_final_4d[j]=0;
        aux_m_SPL_LEdB3d_final_4d[j]=0;
        aux_m_SPL_LEdBAW3d_final_4d[j]=0;
        aux_m_SPL_LEdBBW3d_final_4d[j]=0;
        aux_m_SPL_LEdBCW3d_final_4d[j]=0;
        auxa_m_SPLadB3d_final_4d[j]=0;
        auxa_m_SPLsdB3d_final_4d[j]=0;
        auxa_m_SPLpdB3d_final_4d[j]=0;
        auxa_m_SPLdB3d_final_4d[j]=0;
        auxa_m_SPLdBAW3d_final_4d[j]=0;
        auxa_m_SPLdBBW3d_final_4d[j]=0;
        auxa_m_SPLdBCW3d_final_4d[j]=0;
        auxa_m_SPL_LEdB3d_final_4d[j]=0;
        auxa_m_SPL_LEdBAW3d_final_4d[j]=0;
        auxa_m_SPL_LEdBBW3d_final_4d[j]=0;
        auxa_m_SPL_LEdBCW3d_final_4d[j]=0;

        for (int blade=0;blade<blades_num;++blade){
        for (int E=0;E<angles_num;++E){
       for (unsigned int i = 0; i < number_of_segments; ++i) {
           if (SPLadB3d_4d()[i][j][blade][E]==0.){auxa_m_SPLadB3d_final_4d[j] += 0;} else {auxa_m_SPLadB3d_final_4d[j] += pow(10.,(SPLadB3d_4d()[i][j][blade][E]/10.));}
           if (SPLsdB3d_4d()[i][j][blade][E]==0.){auxa_m_SPLsdB3d_final_4d[j] += 0;} else {auxa_m_SPLsdB3d_final_4d[j] += pow(10.,(SPLsdB3d_4d()[i][j][blade][E]/10.));}
           if (SPLpdB3d_4d()[i][j][blade][E]==0.){auxa_m_SPLpdB3d_final_4d[j] += 0;} else {auxa_m_SPLpdB3d_final_4d[j] += pow(10.,(SPLpdB3d_4d()[i][j][blade][E]/10.));}
           if (SPLdB3d_4d()[i][j][blade][E]==0.){auxa_m_SPLdB3d_final_4d[j] += 0;} else {auxa_m_SPLdB3d_final_4d[j] += pow(10.,(SPLdB3d_4d()[i][j][blade][E]/10.));}
           if (SPLdBAW3d_4d()[i][j][blade][E]==0.){auxa_m_SPLdBAW3d_final_4d[j] += 0;} else {auxa_m_SPLdBAW3d_final_4d[j] += pow(10.,(SPLdBAW3d_4d()[i][j][blade][E]/10.));}
           if (SPLdBBW3d_4d()[i][j][blade][E]==0.){auxa_m_SPLdBBW3d_final_4d[j] += 0;} else {auxa_m_SPLdBBW3d_final_4d[j] += pow(10.,(SPLdBBW3d_4d()[i][j][blade][E]/10.));}
           if (SPLdBCW3d_4d()[i][j][blade][E]==0.){auxa_m_SPLdBCW3d_final_4d[j] += 0;} else {auxa_m_SPLdBCW3d_final_4d[j] += pow(10.,(SPLdBCW3d_4d()[i][j][blade][E]/10.));}
           if (SPL_LEdB3d_4d()[i][j][blade][E]==0.){auxa_m_SPL_LEdB3d_final_4d[j] += 0;} else {auxa_m_SPL_LEdB3d_final_4d[j] += pow(10.,(SPL_LEdB3d_4d()[i][j][blade][E]/10.));}
           if (SPL_LEdBAW3d_4d()[i][j][blade][E]==0.){auxa_m_SPL_LEdBAW3d_final_4d[j] += 0;} else {auxa_m_SPL_LEdBAW3d_final_4d[j] += pow(10.,(SPL_LEdBAW3d_4d()[i][j][blade][E]/10.));}
           if (SPL_LEdBBW3d_4d()[i][j][blade][E]==0.){auxa_m_SPL_LEdBBW3d_final_4d[j] += 0;} else {auxa_m_SPL_LEdBBW3d_final_4d[j] += pow(10.,(SPL_LEdBBW3d_4d()[i][j][blade][E]/10.));}
           if (SPL_LEdBCW3d_4d()[i][j][blade][E]==0.){auxa_m_SPL_LEdBCW3d_final_4d[j] += 0;} else {auxa_m_SPL_LEdBCW3d_final_4d[j] += pow(10.,(SPL_LEdBCW3d_4d()[i][j][blade][E]/10.));}
        }}}}

        if (number_of_segments>sizea){size = number_of_segments;} else {size = sizea;}

        for (unsigned int i=0;i<size;++i){
        for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
        //    3d curves
            m_SPLadB3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLadB3d_final_4d[j]);
            m_SPLsdB3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLsdB3d_final_4d[j]);
            m_SPLpdB3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLpdB3d_final_4d[j]);
            m_SPLdB3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdB3d_final_4d[j]);
            m_SPLdBAW3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBAW3d_final_4d[j]);
            m_SPLdBBW3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBBW3d_final_4d[j]);
            m_SPLdBCW3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPLdBCW3d_final_4d[j]);
            if (m_parameter->Lowson_type!=0){
            m_SPL_LEdB3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdB3d_final_4d[j]);
            m_SPL_LEdBAW3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBAW3d_final_4d[j]);
            m_SPL_LEdBBW3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBBW3d_final_4d[j]);
            m_SPL_LEdBCW3d_final_rotor_loops[i][j]=10*log10((1./number_of_segments)*auxa_m_SPL_LEdBCW3d_final_4d[j]);
            }
            else{
            m_SPL_LEdB3d_final_rotor_loops[i][j]=0;
            m_SPL_LEdBAW3d_final_rotor_loops[i][j]=0;
            m_SPL_LEdBBW3d_final_rotor_loops[i][j]=0;
            m_SPL_LEdBCW3d_final_rotor_loops[i][j]=0;
            }
        }}

        //OASPL complete for quasi 3d
        int i=number_of_segments-1;
        for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
            Final_qs3d_alpha_aux_4d += pow(10.,(m_SPLadB3d_final_rotor_loops[i][j])/10.);
            Final_qs3d_S_aux_4d += pow(10.,(m_SPLsdB3d_final_rotor_loops[i][j])/10.);
            Final_qs3d_P_aux_4d += pow(10.,(m_SPLpdB3d_final_rotor_loops[i][j])/10.);
            Final_qs3d_LE_aux_4d += pow(10.,(m_SPL_LEdB3d_final_rotor_loops[i][j])/10.);
            Final_qs3d_aux_4d += pow(10.,(m_SPLdB3d_final_rotor_loops[i][j])/10.);
                }

            Final_qs3d_alpha_rotor_loops = 10*log10(Final_qs3d_alpha_aux_4d);
            Final_qs3d_S_rotor_loops = 10*log10(Final_qs3d_S_aux_4d);
            Final_qs3d_P_rotor_loops = 10*log10(Final_qs3d_P_aux_4d);
            Final_qs3d_LE_rotor_loops =  10*log10(Final_qs3d_LE_aux_4d);
            Final_qs3d_rotor_loops = 10*log10(Final_qs3d_aux_4d);

        //calculation for the OASPL for the csv output file
        for (unsigned int i=0;i<size;++i){
            auxa_m_OASPL3d_rotor_loops[i]=0;
            auxa_m_OASPLA3d_rotor_loops[i]=0;
            auxa_m_OASPLB3d_rotor_loops[i]=0;
            auxa_m_OASPLC3d_rotor_loops[i]=0;
            auxa_m_SPLALOG3d_rotor_loops[i]=0;
            auxa_m_SPLSLOG3d_rotor_loops[i]=0;
            auxa_m_SPLPLOG3d_rotor_loops[i]=0;
            auxa_m_SPLLEdBAW3d_rotor_loops[i]=0;
            auxa_m_SPLLEdBBW3d_rotor_loops[i]=0;
            auxa_m_SPLLEdBCW3d_rotor_loops[i]=0;
            auxa_m_SPLlogLE3d_rotor_loops[i]=0;

            for (int blade=0;blade<blades_num;++blade){
            for (int E=0;E<angles_num;++E){
        for (int j= 0; j< FREQUENCY_TABLE_SIZE;++j){
            if (SPLadB3d_4d()[i][j][blade][E]==0.){auxa_m_SPLALOG3d_rotor_loops[i] += 0;} else {auxa_m_SPLALOG3d_rotor_loops[i] += pow(10.,(SPLadB3d_4d()[i][j][blade][E]/10.));}
            if (SPLsdB3d_4d()[i][j][blade][E]==0.){auxa_m_SPLSLOG3d_rotor_loops[i] += 0;} else {auxa_m_SPLSLOG3d_rotor_loops[i] += pow(10.,(SPLsdB3d_4d()[i][j][blade][E]/10.));}
            if (SPLpdB3d_4d()[i][j][blade][E]==0.){auxa_m_SPLPLOG3d_rotor_loops[i] += 0;} else {auxa_m_SPLPLOG3d_rotor_loops[i] += pow(10.,(SPLpdB3d_4d()[i][j][blade][E]/10.));}
            if (SPLdB3d_4d()[i][j][blade][E]==0.){auxa_m_OASPL3d_rotor_loops[i] += 0;} else {auxa_m_OASPL3d_rotor_loops[i]+= pow(10.,(SPLdB3d_4d()[i][j][blade][E]/10.));}
            if (SPLdBAW3d_4d()[i][j][blade][E]==0.){auxa_m_OASPLA3d_rotor_loops[i] += 0;} else {auxa_m_OASPLA3d_rotor_loops[i] += pow(10.,(SPLdBAW3d_4d()[i][j][blade][E]/10.));}
            if (SPLdBBW3d_4d()[i][j][blade][E]==0.){auxa_m_OASPLB3d_rotor_loops[i] += 0;} else {auxa_m_OASPLB3d_rotor_loops[i] += pow(10.,(SPLdBBW3d_4d()[i][j][blade][E]/10.));}
            if (SPLdBCW3d_4d()[i][j][blade][E]==0.){auxa_m_OASPLC3d_rotor_loops[i] += 0;} else {auxa_m_OASPLC3d_rotor_loops[i] += pow(10.,(SPLdBCW3d_4d()[i][j][blade][E]/10.));}
            if (SPL_LEdB3d_4d()[i][j][blade][E]==0.){auxa_m_SPLlogLE3d_rotor_loops[i] += 0;} else {auxa_m_SPLlogLE3d_rotor_loops[i] += pow(10.,(SPL_LEdB3d_4d()[i][j][blade][E]/10.));}
            if (SPL_LEdBAW3d_4d()[i][j][blade][E]==0.){auxa_m_SPLLEdBAW3d_rotor_loops[i] += 0;} else {auxa_m_SPLLEdBAW3d_rotor_loops[i] += pow(10.,(SPL_LEdBAW3d_4d()[i][j][blade][E]/10.));}
            if (SPL_LEdBBW3d_4d()[i][j][blade][E]==0.){auxa_m_SPLLEdBBW3d_rotor_loops[i] += 0;} else {auxa_m_SPLLEdBBW3d_rotor_loops[i] += pow(10.,(SPL_LEdBBW3d_4d()[i][j][blade][E]/10.));}
            if (SPL_LEdBCW3d_4d()[i][j][blade][E]==0.){auxa_m_SPLLEdBCW3d_rotor_loops[i] += 0;} else {auxa_m_SPLLEdBCW3d_rotor_loops[i] += pow(10.,(SPL_LEdBCW3d_4d()[i][j][blade][E]/10.));}
        }}}
            m_OASPL3d_rotor_loops[i]=10*log10(auxa_m_OASPL3d_rotor_loops[i]);
            m_OASPLA3d_rotor_loops[i]=10*log10(auxa_m_OASPLA3d_rotor_loops[i]);
            m_OASPLB3d_rotor_loops[i]=10*log10(auxa_m_OASPLB3d_rotor_loops[i]);
            m_OASPLC3d_rotor_loops[i]=10*log10(auxa_m_OASPLC3d_rotor_loops[i]);
            m_SPLALOG3d_rotor_loops[i]=10*log10(auxa_m_SPLALOG3d_rotor_loops[i]);
            m_SPLSLOG3d_rotor_loops[i]=10*log10(auxa_m_SPLSLOG3d_rotor_loops[i]);
            m_SPLPLOG3d_rotor_loops[i]=10*log10(auxa_m_SPLPLOG3d_rotor_loops[i]);
            if (m_parameter->Lowson_type!=0){
            m_SPLLEdBAW3d_rotor_loops[i]=10*log10(auxa_m_SPLLEdBAW3d_rotor_loops[i]);
            m_SPLLEdBBW3d_rotor_loops[i]=10*log10(auxa_m_SPLLEdBBW3d_rotor_loops[i]);
            m_SPLLEdBCW3d_rotor_loops[i]=10*log10(auxa_m_SPLLEdBCW3d_rotor_loops[i]);
            m_SPLlogLE3d_rotor_loops[i]=10*log10(auxa_m_SPLlogLE3d_rotor_loops[i]);
            }
            else{
            m_SPLLEdBAW3d_rotor_loops[i]=0;
            m_SPLLEdBBW3d_rotor_loops[i]=0;
            m_SPLLEdBCW3d_rotor_loops[i]=0;
            m_SPLlogLE3d_rotor_loops[i]=0;
        }}

ProgressBar(5);//Sara
}

//Sara
void NoiseCalculation::setInitialValues(){
QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
double outer_radius=pbem->m_pTData->OuterRadius;

 //obs x pos
         if ((m_parameter->obs_x_pos==0.) & (m_parameter->obs_y_pos==0.) & (m_parameter->obs_z_pos==0.)){
         m_parameter->obs_x_pos=10;
         m_parameter->obs_y_pos=10;
         double hub_radius=pbem->m_pBlade->m_HubRadius;
         double blade_radius=(outer_radius-hub_radius);
         m_parameter->obs_z_pos=blade_radius/2.;
         }

 //x pos rotor
         if((m_parameter->obs_x_pos_rotor==0.) & (m_parameter->obs_y_pos_rotor==0.) & (m_parameter->obs_z_pos_rotor==0.)){
           m_parameter->obs_x_pos_rotor=0;
           m_parameter->obs_y_pos_rotor=0;
           m_parameter->obs_z_pos_rotor=1.5*outer_radius;
         }

 //time
         if(m_parameter->rotation_type==0){m_parameter->time=60./m_parameter->rot_speed;}
         if (m_parameter->time==0){m_parameter->time=60./m_parameter->rot_speed;}

 //timesteps
         if(m_parameter->rotation_type==1){m_parameter->timesteps=(m_parameter->time/(60./m_parameter->rot_speed))/360.*m_parameter->anglesteps*1000.;}
         if (m_parameter->timesteps==0){m_parameter->timesteps=(m_parameter->time/(60./m_parameter->rot_speed))/360.*45.*1000.;}

 //anglesteps
         if(m_parameter->rotation_type==1){m_parameter->anglesteps=m_parameter->timesteps*60.*360./(m_parameter->rot_speed*1000.);}

 //number of loops
 if(m_parameter->rotation_type==1){m_parameter->number_loops=m_parameter->time/(60./m_parameter->rot_speed);}

     //TSR w and u calculation
     double m_TSR_calc=2.*PI*m_parameter->rot_speed/60.*outer_radius/m_parameter->u_wind_speed;

     double m_rot_speed_calc=m_parameter->TSRtd*m_parameter->u_wind_speed*60./(2.*PI*outer_radius);

     double m_u_wind_speed_calc=2.*PI*m_parameter->rot_speed/60.*outer_radius/m_parameter->TSRtd;

 //    calculation for non sets
     if(!m_parameter->u_wind_speed_check){m_parameter->u_wind_speed=m_u_wind_speed_calc;}

     if(!m_parameter->rot_speed_check){m_parameter->rot_speed=m_rot_speed_calc;}

     if(!m_parameter->TSR_check){
     SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
     double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
     double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
     double lend  =   pSimuWidget->m_pctrlLELineEdit->getValue();

 if(m_TSR_calc>lend){
     m_TSR_calc=lend;
     m_parameter->TSRtd=m_TSR_calc;
     m_parameter->rot_speed=m_rot_speed_calc;}

 else if(m_TSR_calc<lstart){
     m_TSR_calc=lstart;
     m_parameter->TSRtd=m_TSR_calc;
     m_parameter->rot_speed=m_rot_speed_calc;}

     else{
 int f=(lend-lstart)/ldelta;
 double m=0;
 double x=0;

     for (double i=0;i<=f;++i){
         x=lstart+i*ldelta;
         m=m_TSR_calc-x;
         if(m<ldelta/2.){break;}
     }

 if(m_TSR_calc==x){
     m_parameter->TSRtd=m_TSR_calc;}
 else{
     double m_TSR_calc_aux=x;

     double m_rot_speed_calc_aux=m_TSR_calc_aux*m_parameter->u_wind_speed*60./(2.*PI*outer_radius);

     double m_u_wind_speed_calc_aux=2.*PI*m_parameter->rot_speed/60.*outer_radius/m_TSR_calc_aux;

 m_parameter->TSRtd=m_TSR_calc_aux;
 double delta_u=qAbs(m_parameter->u_wind_speed-m_u_wind_speed_calc_aux);
 double delta_w=qAbs(m_parameter->rot_speed-m_rot_speed_calc_aux);

 if(delta_w<=delta_u){m_parameter->rot_speed=m_rot_speed_calc_aux;}else{m_parameter->u_wind_speed=m_u_wind_speed_calc_aux;}
 }}}
}

void NoiseCalculation::ProgressBar(int index){
    NoiseSimulation *pNoiseSimulation = (NoiseSimulation *) g_mainFrame->m_pBEM;
    pNoiseSimulation->progress_dlg_canceled=false;
        int w=0;
        if (m_parameter->qs3DSim==0){w=1000000/5;}
        else if (m_parameter->qs3DSim==1){w=1000000/3;}
        else {w=1000000;}
        NoiseCreatorDialog *pNoiseCreatorDialog = (NoiseCreatorDialog *) g_mainFrame->m_pBEM;
if(pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){
    pNoiseCreatorDialog->m_progress_dlg->cancel();
}else{
if (index == 1) {
    for(int j = 0; j<=w; ++j)
        {
            pNoiseCreatorDialog->m_progress_dlg->setValue(j);
            if(pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){
                pNoiseSimulation->progress_dlg_canceled=true;
                pNoiseCreatorDialog->m_progress_dlg->cancel();
                break;
            }
        }}
else {
        for(int j = (index-1)*w+1; j<=index*w; ++j)
            {
                pNoiseCreatorDialog->m_progress_dlg->setValue(j);
                if(pNoiseCreatorDialog->m_progress_dlg->wasCanceled()){
                    pNoiseSimulation->progress_dlg_canceled=true;
                    pNoiseCreatorDialog->m_progress_dlg->cancel();
                    break;
                }
            }
}}}
    //Sara
