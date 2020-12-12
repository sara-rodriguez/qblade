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
#include <QProgressDialog>//Sara
#include "../MainFrame.h"//Sara
class NoiseModule;

class NoiseCreatorDialog : public CreatorDialog, public ParameterViewer<Parameter::NoiseSimulation>
{
    Q_OBJECT

public:

    NoiseCreatorDialog(NoiseSimulation *presetSimulation, NoiseModule *module);

//Sara
    void OnProgressDlg();
    QProgressDialog *m_progress_dlg;

    double change_TSR;
    double AlphaDelta;
    double Mach_initial=0.18;
    double u_wind_speed=12;
    double rot_speed_value = 16;
    double TSR_val_in=7;
    double outer_radius=50;
    double blade_radius=46;

    int phi_selection;
    int sum=2;
    int user_sel;
    int foilindex=g_mainFrame->m_pctrlFoil->currentIndex();

    bool rot_speed_in=false;
    bool u_wind_speed_in=true;
    bool TSR_in=true;
    bool shear_in=false;
    bool shear_roughness_in=false;
    bool shear_height_in=false;
    bool shear_speed_in=false;
    bool rot_in=true;
    bool loops_in = true;
    bool time_in = false;
    bool timesteps_in = false;
    bool delta_user_in = false;
    bool anglesteps_in=true;
    bool check_qs3D=false;
    bool blunt_in=false;
    bool tipvortex_in=false;
    bool check_TE;
    bool check_LE;
    bool op_points_qs3d;

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
    QCheckBox *m_rot_speed_check, *m_TSR_check, *m_u_wind_speed_check, *m_valRel_TE_check, *m_valReu_TE_check, *m_valRel_LE_check, *m_valReu_LE_check, *m_TE_a_check, *m_TE_b_check, *m_TE_c_check, *m_valMal_TE_check, *m_valMau_TE_check, *m_valAOAl_TE_check, *m_valAOAu_TE_check, *m_valMal_LE_check, *m_valMau_LE_check, *m_LBLVS_check, *m_autopolars_check, *m_hblunt_check;

    QGroupBox *m_blunt_check, *m_shear_check, *m_tipvortex_check;

    QComboBox *dstar_combobox, *mode_combobox, *rotation_combobox, *qs3DSim_combobox, *w_TSR_combobox;

    QPushButton *buttonle, *button_cancel, *all_op_points;

    QDoubleSpinBox *m_TSR_spinbox, *m_initial_azimuth_spinbox, *m_yaw_angle_spinbox;

    QRadioButton *one_polar_radiobutton, *multi_polars_radiobutton, *BPM_radiobutton;

    NumberEdit *m_rot_speed_numberedit, *m_anglesteps_numberedit, *m_tower_height_numberedit, *m_u_wind_speed_numberedit, *m_tower_to_hub_distance_numberedit, *m_number_loops_numberedit, *m_time_numberedit, *m_timesteps_numberedit, *m_shear_roughness_numberedit, *m_shear_height_numberedit, *m_shear_speed_numberedit, *m_valRel_TE_numberedit, *m_valReu_TE_numberedit, *m_valMal_TE_numberedit, *m_valMau_TE_numberedit, *m_valAOAl_TE_numberedit, *m_valAOAu_TE_numberedit, *m_valRel_LE_numberedit, *m_valReu_LE_numberedit, *m_valMal_LE_numberedit, *m_valMau_LE_numberedit, *m_hblunt_numberedit;

    int i=0;
    bool all_oppoints_checked;//Sara
//Sara end

    QButtonGroup *m_selectFromButtons;
    QButtonGroup *m_selectFromButtonsLE;
    FoilComboBox *m_airfoilComboBox;
    PolarComboBox *m_polarComboBox;
    QScrollArea *m_opPointScrollArea;
    QWidget *m_opPointViewWidget;
    QWidget *m_originalBpmWidget;

    //Sara begin
    QList< QList<double> > m_CsvParameters;
    QList<QString> m_CsvFileHeader;
    QLabel *m_CsvFileLabel, *m_StarredDLabel, *m_valReLabel, *m_valMaLabel, *m_valAOALabel;
    int m_StarredDType;

private slots:
    void onSelectButtonsClicked (int id);
    void onPolarBoxChange ();
    void onFoilBoxChange ();//Sara
    void onAllButtonToggled (bool pressed);

    void onCreateButtonClicked ();
    void onVerifyWindfield(); //Sara
    void onUnitsChanged () { }  // no need for this

    //Sara begin
        void OnImportStarredD();
        void OnSetDstarButton(int index);
        void OnModeDefine(int index);
        void OnRotationDefine(int index);
        void OnRotSpeedCheck(bool index);
        void OnWindSpeedCheck(bool index);
        void OnTSRCheck(bool index);
        void OnWarningSet3();
        void OnShearLayerCheck(bool index);
        void OnBluntCheck(bool index);
        void OnhBluntCheck(bool index);
        void OnTipVortexCheck(bool index);
        void OnTECheck();
        //Sara end
};

#endif // NOISECREATORDIALOG_H
