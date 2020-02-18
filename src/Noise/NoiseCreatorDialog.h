#ifndef NOISECREATORDIALOG_H
#define NOISECREATORDIALOG_H

class QGridLayout;
class QScrollArea;

#include "../CreatorDialog.h"
#include "../ParameterViewer.h"
#include "../StoreAssociatedComboBox_include.h"
#include "../XBEM/TBEMData.h" //Sara
#include <QCheckBox>//Sara
#include <QComboBox>//Sara
class NoiseModule;

class NoiseCreatorDialog : public CreatorDialog, public ParameterViewer<Parameter::NoiseSimulation>
{
	Q_OBJECT	
	
public:

    NoiseCreatorDialog(NoiseSimulation *presetSimulation, NoiseModule *module);
//Sara
    double change_TSR;
    int phi_selection;
    double AlphaDelta;

    QVector<int> a_D_starred_index_user;
    QVector<double> a_D_starred_S_user;
    QVector<double> a_D_starred_P_user;
    //Sara

private:
	class OpPointRecord {  // a class representing one row in the OpPoint selection view
	public:
		OpPointRecord (QString airfoil, QString polar, QString name, OpPoint *opPoint, bool readyForNoise)
			: airfoil(airfoil), polar(polar), name(name), opPoint(opPoint), checkBox(NULL), readyForNoise(readyForNoise) { }
		static bool sort(OpPointRecord &first, OpPointRecord &second) {
			if (first.airfoil == second.airfoil) {
				if (first.polar == second.polar) {
					return first.name.toDouble() < second.name.toDouble();
				} else {
					return first.polar < second.polar;
				}
			} else {
				return first.airfoil < second.airfoil;
			}
		}
		QString airfoil, polar, name;
		OpPoint *opPoint;
		QCheckBox *checkBox;
		bool readyForNoise;
	};
	
	void initView ();
	void prepareOpPointRecords(bool allPolars);
	void fillOpPointView();

	NoiseModule *m_module;
	NoiseSimulation *m_editedSimulation;
	QList<OpPointRecord> m_opPointRecords;  // store a sorted list of all OpPoints

    //Sara
    QCheckBox *m_rot_speed_check;
    NumberEdit *m_rot_speed_numberedit;
    QCheckBox *m_TSR_check;
    QDoubleSpinBox *m_TSR_spinbox;
    QCheckBox *m_u_wind_speed_check;
    NumberEdit *m_u_wind_speed_numberedit;
    QPushButton *buttonle;
    QComboBox *dstar_combobox;
    QComboBox *mode_combobox;
//Sara

	QButtonGroup *m_selectFromButtons;
    QButtonGroup *m_selectFromButtonsLE;//Sara
	FoilComboBox *m_airfoilComboBox;
    FoilComboBox *m_airfoilComboBoxtd;//Sara
    FoilComboBox *m_airfoilComboBoxte;//Sara
	PolarComboBox *m_polarComboBox;
	QScrollArea *m_opPointScrollArea;
	QWidget *m_opPointViewWidget;
	QWidget *m_originalBpmWidget;
	
    //Sara
    QList< QList<double> > m_CsvParameters;
    QList<QString> m_CsvFileHeader;
    QLabel *m_CsvFileLabel, *m_StarredDLabel;
    int m_StarredDType;
//Sara

private slots:
	void onSelectButtonsClicked (int id);
	void onPolarBoxChange ();
	void onAllButtonToggled (bool pressed);
	
	void onCreateButtonClicked ();
    void onVerifyDeltaFor3D();
	void onUnitsChanged () { }  // no need for this

    //Sara
        void OnImportStarredD();
        void OnSetDstarButton(int index);
        void OnModeDefine(int index);
        void OnRotSpeedCheck();
        void OnWindSpeedCheck();
        void OnTSRCheck();
        void OnWarningSet3();
        //Sara
};

#endif // NOISECREATORDIALOG_H
