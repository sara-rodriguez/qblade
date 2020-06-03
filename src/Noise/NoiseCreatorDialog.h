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

    //Sara begin
    QCheckBox *m_rot_speed_check;
    QCheckBox *m_shear_check;
    QCheckBox *m_TSR_check;
    QCheckBox *m_u_wind_speed_check;

    QComboBox *dstar_combobox;
    QComboBox *mode_combobox;
    QComboBox *rotation_combobox;
    QComboBox *qs3DSim_combobox;

    QPushButton *buttonle;

    QDoubleSpinBox *m_TSR_spinbox;
    QDoubleSpinBox *m_initial_azimuth_spinbox;
    QDoubleSpinBox *m_yaw_angle_spinbox;

    NumberEdit *m_rot_speed_numberedit;
    NumberEdit *m_anglesteps_numberedit;
    NumberEdit *m_tower_height_numberedit;
    NumberEdit *m_u_wind_speed_numberedit;
    NumberEdit *m_tower_to_hub_distance_numberedit;
    NumberEdit *m_number_loops_numberedit;
    NumberEdit *m_time_numberedit;
    NumberEdit *m_timesteps_numberedit;
    NumberEdit *m_shear_roughness_numberedit;
    NumberEdit *m_shear_height_numberedit;
    NumberEdit *m_shear_speed_numberedit;
//Sara end

    QButtonGroup *m_selectFromButtons;
    QButtonGroup *m_selectFromButtonsLE;
    FoilComboBox *m_airfoilComboBox;
    FoilComboBox *m_airfoilComboBoxtd;
    FoilComboBox *m_airfoilComboBoxte;
    PolarComboBox *m_polarComboBox;
    QScrollArea *m_opPointScrollArea;
    QWidget *m_opPointViewWidget;
    QWidget *m_originalBpmWidget;

    //Sara begin
    QList< QList<double> > m_CsvParameters;
    QList<QString> m_CsvFileHeader;
    QLabel *m_CsvFileLabel, *m_StarredDLabel;
    int m_StarredDType;
//Sara end

private slots:
    void onSelectButtonsClicked (int id);
    void onPolarBoxChange ();
    void onAllButtonToggled (bool pressed);

    void onCreateButtonClicked ();
    void onVerifyDeltaFor3D(); //Sara
    void onVerifyWindfield(); //Sara
    void onUnitsChanged () { }  // no need for this

    //Sara begin
        void OnImportStarredD();
        void OnSetDstarButton(int index);
        void OnModeDefine(int index);
//        void OnStepTypeDefine(int index);
        void OnRotationDefine(int index);
        void OnRotSpeedCheck();
        void OnWindSpeedCheck();
        void OnTSRCheck();
        void OnWarningSet3();
        void OnShearLayerCheck();
        //Sara end
};

#endif // NOISECREATORDIALOG_H
