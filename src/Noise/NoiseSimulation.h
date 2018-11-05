#ifndef NOISESIMULATION_H
#define NOISESIMULATION_H

#include "../StorableObject.h"
#include "../Graph/ShowAsGraphInterface.h"
#include "../ParameterObject.h"
#include "NoiseCalculation.h"
#include "NoiseParameter.h"

#include "../XBEM/BEM.h" //Sara

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
    //Sara
};

#endif // NOISESIMULATION_H
