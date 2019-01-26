#include "NoiseSimulation.h"

#include "../ParameterViewer.h"
#include "../Store.h"
#include "../Objects/Polar.h"
#include "../Objects/OpPoint.h"
#include "../Graph/NewCurve.h"
#include "NoiseModule.h"
#include "../ColorManager.h"
#include "NoiseOpPoint.h"
#include "../XBEM/BEM.h"


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
//	curve->setAllPoints(xVector.data(), yVector.data(), xVector.size());
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

    stream <<   qSetFieldWidth(14)  <<
              "Sect"  << ";" <<
              "Radius [m]"  << ";" <<
              "r/R"  << ";" <<
//              "Blade Section" << ";" <<
              "Chord [m]" << ";" <<
              "Theta [deg]" << ";" <<
              "Axial Ind. Fact. (a)"  << ";" <<
              "Axial Velocity [m/s]" << ";" <<
              "Tg. Ind. Fact. [a']"  << ";" <<
              "Tg. Speed [m/s]" << ";" <<
              "Res. Local Speed Calc [m/s]" << ";" <<
              "Res. Local Speed BEM [m/s]" << ";" <<
              "Re BEM"  << ";" <<
              "Re calc"  << ";" <<
              "Mach BEM"  << ";" <<
              "Mach calc"  << ";" <<
              "(cl/cd) max"   << ";" <<
              "cl"   << ";" <<
              "cd"   << ";" <<
              "(cl/cd) max angle"    <<";"   <<
              "c/R  "<<";"   <<
              endl;

        for (int i = 0; i < number_of_segments; ++i) {

            // definitions
            double axial_ind_fact = bdata->m_a_axial.value(i);
            double axial_velocity = approaxing_wind_speed*(1-axial_ind_fact);

            double tangential_speed = omega*bdata->m_pos.value(i)*(1+bdata->m_a_tangential.value(i));
            double resultant_local_speed = qSqrt(pow(axial_velocity,2)+pow(tangential_speed,2));
            double chord = bdata->m_c_local.value(i);
            double Reynolds_calc = rho*resultant_local_speed*chord/dynamic_visc;
            double Mach_calc = resultant_local_speed/sqrt(R_air*K_air*T_std_cond);
            double alpha = bdata->m_alpha.value(i);
            double cl_cd =  bdata->m_LD.value(i);
            double r_R = bdata->m_pos.value(i)/finalradius;

            double c_Rx = 0;
            double r_R0  =  0.05; double c_R0 = 0.05500;
            double r_R1  =  0.25; double c_R1 = 0.07500;
            double r_R2  =  1.00; double c_R2 = 0.02000;

            if (r_R <= r_R0) {c_Rx = c_R0;}
            if (r_R > r_R0 && r_R < r_R1) {c_Rx = (r_R-r_R0)*(c_R1-c_R0)/(0.25-r_R0)+c_R0;}
            if (r_R <= r_R1 && r_R >= r_R1) {c_Rx = c_R1;}
            if (r_R > r_R1 && r_R < r_R2) {c_Rx = (r_R-r_R1)*(c_R2-c_R1)/(r_R2-r_R1)+c_R1;}
            if (r_R >= r_R2) {c_Rx = c_R2;}

            QString c_R= QString::number(c_Rx, 'f', 5);

        stream << qSetFieldWidth(14)  <<
                      (i+1) << ";" <<
                      bdata->m_pos.value(i) << ";" <<
                      r_R << ";" <<
//                      "HOLD" << ";" << //defaultName
                      chord << ";" <<
                      bdata->m_theta.value(i) << ";" <<
                      axial_ind_fact << ";" <<
                      axial_velocity << ";" <<
                      bdata->m_a_tangential.value(i) << ";" <<
                      tangential_speed << ";" <<
                      resultant_local_speed << ";" <<
                      bdata->m_Windspeed.value(i) << ";" <<
                      bdata->m_Reynolds.value(i) << ";" <<
                      Reynolds_calc << ";" <<
                      bdata->m_Mach.value(i) << ";" <<
                      Mach_calc << ";" <<
                      cl_cd << ";" <<
                      bdata->m_CL.value(i) << ";" <<
                      bdata->m_CD.value(i) << ";" <<
                      alpha <<  ";" <<
                      c_R  << ";" << endl;
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

        stream << endl;
        stream << endl;
        stream << "With Prandtl's Tip Loss Correction and Glauert Correction" << endl;
        stream << "BEM Method. Hansen 2008 p. 50" << endl;
        stream << endl;
        stream << "Sect Defined:"   <<  (sectionx+1);
        stream << endl;
        stream << endl;

        stream <<   qSetFieldWidth(14)  <<
        "Sect"  << ";" <<
        "a"  << ";" <<
        "a (a>0.2)"  << ";" <<
        "a'"  << ";" <<
        "a var [%]"  << ";" <<
        "a' var [%]"  << ";" <<
        "(1-a)"  << ";" <<
        "(1-a')"  << ";" <<
        "V0"  << ";" <<
        "Wr"  << ";" <<
        "phi [deg]"  << ";" <<
        "["  << ";" <<
        "theta [deg]"  << ";" <<
        "cl"  << ";" <<
        "cd"  << ";" <<
        "cn"  << ";" <<
        "ct"  << ";" <<
        "?"  << ";" <<
        "f"  << ";" <<
        "F"  << ";" <<
        "K"  << ";" <<
        endl;

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
if (i==-1){
    a_calc[-1]=0;
    a_corrected_calc[-1]=0;
    a_lin_calc[-1]=0;
    a_variation_calc[-1]=0;
    a_lin_variation_calc[-1]=0;
    phi_calc[-1]=(qRadiansToDegrees(qAtan((one_a_calc[-1]*(double)approaxing_wind_speed)/(one_a_lin_calc[-1]*(double)wr))));
    F_calc[-1]=0;
    K_calc[-1]=0;
    one_a_calc[-1]=1;
    one_a_lin_calc[-1]=1;
    theta_calc[-1]=(phi_calc[-1]-lbrac);
    cn_calc[-1]=(cl*qCos(qDegreesToRadians(phi_calc[-1]))+cd*(double)qSin(qDegreesToRadians(phi_calc[-1])));
    ct_calc[-1]=(cl*qSin(qDegreesToRadians(phi_calc[-1]))-cd*(double)qCos(qDegreesToRadians(phi_calc[-1])));
    f_calc[-1]=(3*(double)(50-radiusx)/(double)(2*radiusx*(double)qSin(qDegreesToRadians(phi_calc[-1]))));
    Fma_calc[-1]=(2/M_PI*(1/(double)qCos(qExp(-f_calc[-1]))));
    K_calc[-1]=(4*(double)Fma_calc[-1]*pow(qSin(qDegreesToRadians(phi_calc[-1])),2)/doubt_calc*cn_calc[-1]);
}
else{
        a_lin_calc[i] = (1.f/((4.f*Fma_calc[i-1]*qSin(qDegreesToRadians(phi_calc[i-1]))*qCos(qDegreesToRadians(phi_calc[i-1]))/doubt_calc*ct_calc[i-1])-1.f));

        a_calc[i] = (1/(double)(K_calc[i-1]+1));
        a_corrected_calc[i] = (0.5*(2+K_calc[i-1]*(double)(1-2*ac)-qSqrt(pow((K_calc[i-1]*(1-2*ac)+2),2)+4*(double)(K_calc[i-1]*pow(ac,2)-1))));

        if (i==0){a_variation_calc[i] = ((a_calc[i]-a_calc[i-1])*(double)100);} else {
        a_variation_calc[i] = (((a_calc[i]-a_calc[i-1])/a_calc[i-1])*(double)100);}

        if (i==0){a_lin_variation_calc[i] = ((a_lin_calc[i]-a_lin_calc[i-1])*100);} else {
        a_lin_variation_calc[i] = (((a_lin_calc[i]-a_lin_calc[i-1])/(double)a_lin_calc[i-1])*100);}

        if (a_calc[i]<ac){one_a_calc[i]=(1-a_calc[i]);} else {one_a_calc[i]=(1-a_corrected_calc[i]);}

        one_a_lin_calc[i]=(1+a_lin_calc[i]);

        phi_calc[i]=(qRadiansToDegrees(qAtan((one_a_calc[i]*approaxing_wind_speed)/(double)(one_a_lin_calc[i]*wr))));

        theta_calc[i]=(phi_calc[i]-lbrac);

        cn_calc[i] = (cl*qCos(qDegreesToRadians(phi_calc[i]))+cd*qSin(qDegreesToRadians(phi_calc[i])));

        ct_calc[i] = (cl*qSin(qDegreesToRadians(phi_calc[i]))-cd*qCos(qDegreesToRadians(phi_calc[i])));

        f_calc[i]=(3*(50-radiusx)/(double)(2*radiusx*qSin(qDegreesToRadians(phi_calc[i]))));

        Fma_calc[i]=(2/M_PI*(1/(double)qCos(qExp(-f_calc[i]))));

        K_calc[i]=(4*Fma_calc[i]*pow(qSin(qDegreesToRadians(phi_calc[i])),2)/(double)doubt_calc*cn_calc[i]);
}

        stream << qSetFieldWidth(14)  <<
                  (i+1) << ";" <<
                  a_calc[i] << ";" <<
                  a_corrected_calc[i] << ";" <<
                  a_lin_calc[i] << ";" <<
                  a_variation_calc[i] << ";" <<
                  a_lin_variation_calc[i] << ";" <<
                  one_a_calc[i] << ";" <<
                  one_a_lin_calc[i] << ";" <<
                  approaxing_wind_speed << ";" <<
                  wr << ";" <<
                  phi_calc[i] << ";" <<
                  lbrac << ";" <<
                  theta_calc[i] << ";" <<
                  cl << ";" <<
                  cd << ";" <<
                  cn_calc[i] << ";" <<
                  ct_calc[i] << ";" <<
                  doubt_calc << ";" <<
                  f_calc[i] << ";" <<
                  Fma_calc[i] << ";" <<
                  K_calc[i] << ";" <<
                endl;
}
                stream << endl;
                stream << endl;
                stream << endl;
}
//       delete[] variables;
       qDeleteAll(noiseOpPoints);
}


//teste

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
                if(set) m_parameter.sects = value.toDouble();
                if(m_parameter.sects<13) m_parameter.sects=13;
                else value = m_parameter.sects; break;
        //Sara
	}

	return (set ? QVariant() : value);
}
