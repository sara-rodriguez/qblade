#ifndef NOISESIMULATION_H
#define NOISESIMULATION_H

#include "../StorableObject.h"
#include "../Graph/ShowAsGraphInterface.h"
#include "../ParameterObject.h"
#include "NoiseCalculation.h"
#include "NoiseParameter.h"

#include "../XBEM/BData.h" //Sara

template <class KeyType>
class ParameterViewer;

class NoiseSimulation : public StorableObject, public ShowAsGraphInterface, public ParameterObject<Parameter::NoiseSimulation>
{
public:
	static NoiseSimulation* newBySerialize ();
	NoiseSimulation(ParameterViewer<Parameter::NoiseSimulation> *viewer);
	
	void serialize();
	void restorePointers();
	NewCurve* newCurve (QString /*xAxis*/, QString /*yAxis*/, NewGraph::GraphType /*graphType*/) { return NULL; }
	NewCurve* newCurve (QString xAxis, QString yAxis, NewGraph::GraphType graphType, int opPointIndex);
	QString getObjectName () { return m_objectName; }
	static QStringList getAvailableVariables (NewGraph::GraphType graphType = NewGraph::None);
	static QStringList prepareMissingObjectMessage();
	
	void simulate();  // can throw NoiseException
	void exportCalculation (QTextStream &stream);
    void export3DCalculation (QTextStream &stream);//Sara
	
	void setAnalyzedOpPoints (QVector<OpPoint*> newList);
	void setSelectFrom (NoiseParameter::OpPointSource select) { m_parameter.opPointSource = select; }
	QVector<OpPoint*> getAnalyzedOpPoints () { return m_parameter.analyzedOpPoints; }
	NoiseParameter::OpPointSource getSelectFrom () { return m_parameter.opPointSource; }
	
private:
	NoiseSimulation () { }
	QPen doGetPen (int curveIndex, int highlightedIndex, bool forTheDot);
	QVariant accessParameter(Parameter::NoiseSimulation::Key key, QVariant value = QVariant());
	
	NoiseParameter m_parameter;
	NoiseCalculation m_calculation;

    //Sara
    double m_DStarFinalS;
    double m_DStarFinalP;
    double originalMach;
    double m_slices;
    double x;

    QList <double> m_pos;
    QList <double> m_c_local;       //local chord
    QList <double> m_lambda_local;  //local lambda
    QList <double> m_p_tangential;  //tangential thrust
    QList <double> m_p_normal;      //radial thrust
    QList <double> m_a_axial;       //axial induction factor
    QList <double> m_a_tangential;      //radial induction factor
    QList <double> m_Fa_axial;      //averaged axial induction factor
    QList <double> m_Fa_radial;     //averaged radial induction factor
    QList <double> m_circ;          //circulation
    QList <double> m_theta;         //angles in the wind triangle
    QList <double> m_alpha;         //angles in the wind triangle
    QList <double> m_phi;           //angles in the wind triangle
    QList <double> m_CL;            //lift coeff
    QList <double> m_CD;            //drag coeff
    QList <double> m_LD;            //lift to drag coeff
    QList <double> m_Cn;            //normal coefficient
    QList <double> m_Ct;            //thrust coefficient
    QList <double> m_F;             //Tip Loss Coefficient
    QList <double> m_Reynolds;      //Reynolds Number
    QList <double> m_DeltaReynolds; //Delta between local Re and Re of Polar
    QList <double>  m_Roughness;    //the critical roughness for the blade surface
    QList <double>  m_Windspeed;    //windspeed at section (turbine calc)
    QList <double>  m_Iterations;   //total number of iterations required to converge
    QList <double>  m_Mach;         //Mach Number
    //Sara
};

#endif // NOISESIMULATION_H
