#include "NoiseCreatorDialog.h"

#include <QGroupBox>
#include <QScrollArea>
#include <QRadioButton>
#include <QButtonGroup>
#include <algorithm>
#include <QMessageBox>

#include "../ParameterGrid.h"
#include "NoiseSimulation.h"
#include "../Store.h"
#include "NoiseModule.h"
#include "../Objects/OpPoint.h"
#include "../Objects/Polar.h"
#include "../Objects/Foil.h"
#include "NoiseException.h"
#include "../Noise/NoiseCreatorDialog.h"

//Sara
#include "../XBEM/BEM.h"
#include "../XDirect/XDirect.h"
#include "../XUnsteadyBEM/WindField.h"
#include "../MainFrame.h"
#include "../XUnsteadyBEM/WindFieldModule.h"
#include "NoiseCalculation.h"
//Sara

typedef Parameter::NoiseSimulation P;

NoiseCreatorDialog::NoiseCreatorDialog(NoiseSimulation *presetSimulation, NoiseModule *module)
    : m_module(module),
      m_editedSimulation(presetSimulation),
      m_opPointViewWidget(NULL)
{
    setWindowTitle("Noise Simulation");	//Sara
    //Sara
    FoilPolarDlg *pFoilPolarDlg = (FoilPolarDlg *) g_mainFrame->m_pBEM;
    double Mach_initial=pFoilPolarDlg->m_pctrlMach->getValue();
    if(Mach_initial<0.01){Mach_initial=0.18;}

//Sara

    QTabWidget *tabWidget = new QTabWidget;
    m_contentVBox->insertWidget(0, tabWidget);

    QWidget *widget = new QWidget;
    tabWidget->addTab(widget, "Parameters");
        QHBoxLayout *hBox = new QHBoxLayout;
        setMinimumSize(950,540);//Sara
        widget->setLayout(hBox);
            QGroupBox *groupBox = new QGroupBox ("Simulation Parameters");
            hBox->addWidget(groupBox);
                ParameterGrid<P> *pGrid = new ParameterGrid<P>(this);
                groupBox->setLayout(pGrid);
                    pGrid->addEdit(P::Name, LineEdit, new QLineEdit, "Name of Simulation:", "Noise Simulation");
                    get<QLineEdit>(P::Name)->setMinimumWidth(150);
                    pGrid->addEdit(P::WettedLength, NumberEditType, new NumberEdit(),
                                  "Length of wetted Trailing-Edge (L) []:", 1, LENGTH);
                    pGrid->addEdit(P::DistanceObsever, NumberEditType, new NumberEdit(),
                                  "Distance from observer to TE (re) []:", 1.22, LENGTH);

                    SimuWidget *pSimuWidget = (SimuWidget *) g_mainFrame->m_pSimuWidget;
                    double u_wind_speed=pSimuWidget->m_pctrlWindspeed->getValue();

                    pGrid->addEdit(P::OriginalVelocity, NumberEditType, new NumberEdit(),
                                  "Original flow velocity (U) []:", u_wind_speed, SPEED);
                    pGrid->addEdit(P::OriginalChordLength, NumberEditType, new NumberEdit(),
                                  "Original airfoil Chord length (C) []:", 1, LENGTH);
                    pGrid->addEdit(P::OriginalMach, NumberEditType, new NumberEdit(),
                                  "Original flow Mach Number (M):", Mach_initial);//0.21 Mach_initial
                    //Alexandre MOD
                    pGrid->addEdit(P::IntegralLengthScale, NumberEditType, new NumberEdit(),
                                  "Turbulence integral length scale (l) []:", 1, LENGTH);
                    pGrid->addEdit(P::TurbulenceIntensity, NumberEditType, new NumberEdit(),
                                  "Turbulence intensity []:", 0.125, PERCENT);
                    //End Alexandre MOD
                    pGrid->addEdit(P::DStarChordStation, NumberEditType, new NumberEdit(),
                                  "δ* at chord station:", 0.98);//Sara
                    pGrid->addEdit(P::DStarScalingFactor, NumberEditType, new NumberEdit(),
                                  "δ* scaling factor:", 1);//Sara
                    pGrid->addEdit(P::EddyConvectionMach, NumberEditType, new NumberEdit(),
                                  "Eddy Convection Mach number []:", 0.8, PERCENT);
                    pGrid->addEdit(P::DirectivityTheta, NumberEditType, new NumberEdit(),
                                  "Directivity angle θe [deg]:", 90);//Sara
                    pGrid->addEdit(P::DirectivityPhi, NumberEditType, new NumberEdit(),
                                  "Directivity angle ψe [deg]:", 90);//Sara


            QVBoxLayout *vBox = new QVBoxLayout;
            hBox->addLayout(vBox);
                QLabel *imageLabel = new QLabel;
                imageLabel->setPixmap(QPixmap(":/images/noise_3d_plate.png"));
                vBox->addWidget(imageLabel, 0, Qt::AlignHCenter);
                groupBox = new QGroupBox ("TE noise source contributions");
                vBox->addWidget(groupBox);
                pGrid = new ParameterGrid<P>(this);
                groupBox->setLayout(pGrid);
                pGrid->addEdit(P::SeparatedFlow, CheckBox, new QCheckBox ("enable"),"Separated flow on the suction side (high Reynolds flow):", true);
                pGrid->addEdit(P::SuctionSide, CheckBox, new QCheckBox ("enable"),"Suction side of airfoil (attached flow):", true);
                pGrid->addEdit(P::PressureSide, CheckBox, new QCheckBox ("enable"),"Pressure side of airfoil (attached flow):", true);

                //Sara
                groupBox = new QGroupBox ("LE noise source contributions");
                vBox->addWidget(groupBox);
                pGrid = new ParameterGrid<P>(this);
                groupBox->setLayout(pGrid);
QComboBox *Lowson_type_combobox = new QComboBox;
pGrid->addEdit(P::Lowson_type,ComboBox, Lowson_type_combobox,"Lowson's Model:","");
Lowson_type_combobox->insertItem(0,"None");
Lowson_type_combobox->insertItem(1,"Von Kármán");
Lowson_type_combobox->insertItem(2,"Rapid Distortion");

groupBox = new QGroupBox ("Quasi 3D Simulation");
vBox->addWidget(groupBox);
pGrid = new ParameterGrid<P>(this);
groupBox->setLayout(pGrid);
QComboBox *qs3DSim_combobox = new QComboBox;
pGrid->addEdit(P::qs3DSim,ComboBox, qs3DSim_combobox,"Quasi 3D:","");
qs3DSim_combobox->insertItem(0,"Rotor");
qs3DSim_combobox->insertItem(1,"Blade");
qs3DSim_combobox->insertItem(2,"Disable");

connect(qs3DSim_combobox, QOverload<int>::of(&QComboBox::currentIndexChanged),
    [=](int index){

    NoiseCalculation *pNoiseCalculation = (NoiseCalculation *) g_mainFrame->m_pBEM;
    pNoiseCalculation->user_sel=index;

if (index == 0){
    tabWidget->setTabEnabled(2, true);
    tabWidget->setTabEnabled(3, true);
}
if (index == 1){
    tabWidget->setTabEnabled(2, true);
    tabWidget->setTabEnabled(3, false);
}
if (index == 2){
    tabWidget->setTabEnabled(2, false);
    tabWidget->setTabEnabled(3, false);
}

if(index!=2){check=true;}

});
                //Sara
                vBox->addStretch();

    widget = new QWidget;
    tabWidget->addTab(widget, "Op. Points");
        hBox = new QHBoxLayout;
        widget->setLayout(hBox);
            groupBox = new QGroupBox ("Operational Points");
            groupBox->setMinimumHeight(300);
            hBox->addWidget(groupBox);
                QGridLayout *grid = new QGridLayout;
                grid->setColumnStretch(5, 1);
                groupBox->setLayout(grid);
                    QLabel *label = new QLabel ("Select operational points from");
                    grid->addWidget(label, 0, 0, 1, 1);
                    m_selectFromButtons = new QButtonGroup;
                    connect(m_selectFromButtons, SIGNAL(buttonClicked(int)), this, SLOT(onSelectButtonsClicked(int)));
                    QRadioButton *radioButton = new QRadioButton ("this polar:");
                    m_selectFromButtons->addButton(radioButton, NoiseParameter::OnePolar);
                    grid->addWidget(radioButton, 0, 1, 1, 1);
                    m_airfoilComboBox = new FoilComboBox (&g_foilStore);
                    grid->addWidget(m_airfoilComboBox, 0, 2, 1, 1);
                    m_polarComboBox = new PolarComboBox (&g_polarStore);
                    m_polarComboBox->setParentBox(m_airfoilComboBox);
                    connect(m_polarComboBox, SIGNAL(valueChanged(int)), this, SLOT(onPolarBoxChange()));
                    grid->addWidget(m_polarComboBox, 0, 3, 1, 1);
                    radioButton = new QRadioButton ("all polars");
                    m_selectFromButtons->addButton(radioButton, NoiseParameter::MultiplePolars);
                    grid->addWidget(radioButton, 1, 1, 1, 3);
                    radioButton = new QRadioButton ("original BPM δ* correlations");
                    m_selectFromButtons->addButton(radioButton, NoiseParameter::OriginalBpm);
                    grid->addWidget(radioButton, 2, 1, 1, 3);

                    m_opPointScrollArea = new QScrollArea;
                    grid->addWidget(m_opPointScrollArea, 3, 0, 1, 4);
                        // scroll area is filled in NoiseCreatorDialog::fillOpPointView

                    m_originalBpmWidget = new QWidget;
                    grid->addWidget(m_originalBpmWidget, 3,0,1,4, Qt::AlignLeft | Qt::AlignTop);
                    pGrid = new ParameterGrid<P>(this);
                    m_originalBpmWidget->setLayout(pGrid);
                    pGrid->addEdit(P::Aoa, NumberEditType, new NumberEdit, "AOA (α) [deg]:", 0);//Sara
                    pGrid->addEdit(P::ChordBasedReynolds, NumberEditType, new NumberEdit,
                                           "Chord based Reynolds number (Rc):", 100000);
                    pGrid->addComboBox(P::Transition, "Type of Transition:", NoiseParameter::TransitionFlow,
                                               QStringList()<<"Fully turbulent"<<"Transition flow");

                            //Sara
                            widget = new QWidget;
                            tabWidget->addTab(widget, "quasi 3D Blade");
                            vBox = new QVBoxLayout;
                            hBox = new QHBoxLayout;
                            widget->setLayout(hBox);
                            groupBox = new QGroupBox ("quasi 3D Blade Simulation Parameters");
                            hBox->addWidget(groupBox);
                            pGrid = new ParameterGrid<P>(this);
                            groupBox->setLayout(pGrid);

m_rot_speed_check = new QCheckBox("rot. speed set:");
pGrid->addEdit(P::rot_speed_check, CheckBox, m_rot_speed_check,"", rot_speed_in);
connect(m_rot_speed_check,SIGNAL(toggled(bool)),this,SLOT(OnRotSpeedCheck(bool)));
m_rot_speed_check->setToolTip("Select 2 of 3");
connect(m_rot_speed_check,SIGNAL(clicked()),this,SLOT(OnWarningSet3()));

m_rot_speed_numberedit = new NumberEdit ();
m_rot_speed_numberedit->setEnabled(m_rot_speed_check->isChecked());
pGrid->addEdit(P::rot_speed, NumberEditType, m_rot_speed_numberedit,"Rotational Speed [rpm]:",6);
m_rot_speed_numberedit->setAutomaticPrecision(3);

m_u_wind_speed_check = new QCheckBox("wind speed set:");
pGrid->addEdit(P::u_wind_speed_check, CheckBox, m_u_wind_speed_check,"", u_wind_speed_in);
connect(m_u_wind_speed_check,SIGNAL(toggled(bool)),this,SLOT(OnWindSpeedCheck(bool)));
m_u_wind_speed_check->setToolTip("Select 2 of 3");
connect(m_u_wind_speed_check,SIGNAL(clicked()),this,SLOT(OnWarningSet3()));

m_u_wind_speed_numberedit = new NumberEdit ();
m_u_wind_speed_numberedit->setEnabled(m_u_wind_speed_check->isChecked());
m_u_wind_speed_numberedit->setAutomaticPrecision(3);
pGrid->addEdit(P::u_wind_speed, NumberEditType, m_u_wind_speed_numberedit,"Uniform Wind Speed []:",u_wind_speed,SPEED);

m_TSR_check = new QCheckBox ("TSR set");
pGrid->addEdit(P::TSR_check, CheckBox, m_TSR_check,"", TSR_in);
connect(m_TSR_check,SIGNAL(toggled(bool)),this,SLOT(OnTSRCheck(bool)));
m_TSR_check->setToolTip("Select 2 of 3");
connect(m_TSR_check,SIGNAL(clicked()),this,SLOT(OnWarningSet3()));

m_TSR_spinbox = new QDoubleSpinBox;
m_TSR_spinbox->setEnabled(m_TSR_check->isChecked());
m_TSR_spinbox->setLocale(QLocale("en_us"));
pGrid->addEdit(P::TSRtd,DoubleSpinBox, m_TSR_spinbox,"TSR:", 7);
double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
double lend  =   pSimuWidget->m_pctrlLELineEdit->getValue();

m_TSR_spinbox->setRange(lstart, lend);
m_TSR_spinbox->setSingleStep(ldelta);
m_TSR_spinbox->setDecimals(1);

dstar_combobox = new QComboBox;
pGrid->addEdit(P::dstar_type,ComboBox, dstar_combobox,"δ* source:","");
dstar_combobox->insertItem(0,"XFoil");
dstar_combobox->insertItem(1,"BPM");
dstar_combobox->insertItem(2,"User");
connect(dstar_combobox,SIGNAL(currentIndexChanged(int)),this,SLOT(OnSetDstarButton(int)));

pGrid->setSizeConstraint(QLayout::SetMinimumSize);
buttonle = new QPushButton ("δ* User Input");
buttonle->setMinimumWidth(QFontMetrics(QFont()).width("δ* User Input") * 1.8);
                            buttonle->setEnabled(delta_user_in);
                            pGrid->addWidget(buttonle,9,1);//,8,2
                            connect(buttonle,SIGNAL(clicked()),this,SLOT(OnImportStarredD()));

                            hBox->addLayout(vBox);
                                QLabel *imageLabela = new QLabel;
                                imageLabela->setPixmap(QPixmap(":/images/noise_3d_position.png"));
                                vBox->addWidget(imageLabela, 0, Qt::AlignHCenter);

                                groupBox = new QGroupBox ("Observer Position");
                                vBox->addWidget(groupBox);
                                pGrid = new ParameterGrid<P>(this);
                                groupBox->setLayout(pGrid);

                                pGrid->addEdit(P::obs_x_pos, NumberEditType, new NumberEdit(),"XB:", 10);
                                pGrid->addEdit(P::obs_y_pos, NumberEditType, new NumberEdit(),"YB:", 10);

                                QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
                                double hub_radius=pbem->m_pBlade->m_HubRadius;
                                double outer_radius=pbem->m_pTData->OuterRadius;
                                double blade_radius=(outer_radius-hub_radius);
                                double z_pos=blade_radius/2.;

                                pGrid->addEdit(P::obs_z_pos, NumberEditType, new NumberEdit(),"ZB:", z_pos);

                                widget = new QWidget;
                                tabWidget->addTab(widget, "quasi 3D Rotor");
                                vBox = new QVBoxLayout;
                                hBox = new QHBoxLayout;
                                widget->setLayout(hBox);
                                groupBox = new QGroupBox ("quasi 3D Rotor Simulation Parameters");
                                hBox->addWidget(groupBox);
                                pGrid = new ParameterGrid<P>(this);
                                groupBox->setLayout(pGrid);
                                mode_combobox = new QComboBox;
                                pGrid->addEdit(P::state_ss_us,ComboBox, mode_combobox,"State:","");
                                mode_combobox->insertItem(0,"Steady");
                                mode_combobox->insertItem(1,"Unsteady");
                                connect(mode_combobox,SIGNAL(currentIndexChanged(int)),this,SLOT(OnModeDefine(int)));

                                double tower_height_in;
                            if(g_windFieldStore.size() == 0){tower_height_in=100-hub_radius;}else{tower_height_in=g_windFieldModule->getShownWindField()->getHubheight()-hub_radius;}
                                m_tower_height_numberedit = new NumberEdit ();
                                pGrid->addEdit(P::tower_height, NumberEditType, m_tower_height_numberedit,"Tower Height []:",tower_height_in,LENGTH);

                                m_initial_azimuth_spinbox = new QDoubleSpinBox;
                                m_initial_azimuth_spinbox->setLocale(QLocale("en_us"));
                                pGrid->addEdit(P::initial_azimuth,DoubleSpinBox, m_initial_azimuth_spinbox,"Initial Azimuth [deg]:", 0);
                                m_initial_azimuth_spinbox->setRange(0, 360);
                                m_initial_azimuth_spinbox->setSingleStep(1);
                                m_initial_azimuth_spinbox->setDecimals(0);

                                m_yaw_angle_spinbox = new QDoubleSpinBox;
                                m_yaw_angle_spinbox->setLocale(QLocale("en_us"));
                                pGrid->addEdit(P::yaw_angle,DoubleSpinBox, m_yaw_angle_spinbox,"Yaw Angle [deg]:", 0);
                                m_yaw_angle_spinbox->setRange(0, 360);
                                m_yaw_angle_spinbox->setSingleStep(1);
                                m_yaw_angle_spinbox->setDecimals(0);

                                m_tower_to_hub_distance_numberedit = new NumberEdit ();
                                pGrid->addEdit(P::tower_to_hub_distance, NumberEditType, m_tower_to_hub_distance_numberedit,"Tower To Hub Distance []:",4,LENGTH);

                                rotation_combobox = new QComboBox;
                                pGrid->addEdit(P::rotation_type,ComboBox, rotation_combobox,"Reference:","");
                                rotation_combobox->insertItem(0,"Angular");
                                rotation_combobox->insertItem(1,"Time");
                                connect(rotation_combobox,SIGNAL(currentIndexChanged(int)),this,SLOT(OnRotationDefine(int)));
                                rotation_combobox->setEnabled(rot_in);

                                m_number_loops_numberedit = new NumberEdit ();
                                pGrid->addEdit(P::number_loops, NumberEditType, m_number_loops_numberedit,"Number of Loops:",1);
                                m_number_loops_numberedit->setEnabled(loops_in);

                                m_anglesteps_numberedit = new NumberEdit ();
                                pGrid->addEdit(P::anglesteps, NumberEditType, m_anglesteps_numberedit,"Anglesteps [deg]:",45);
                                m_anglesteps_numberedit->setMinimum(1);
                                m_anglesteps_numberedit->setMaximum(45);
                                m_anglesteps_numberedit->setEnabled(anglesteps_in);

                                m_time_numberedit = new NumberEdit ();
                                pGrid->addEdit(P::time, NumberEditType, m_time_numberedit,"Time [s]:",100);
                                m_time_numberedit->setEnabled(time_in);

                                m_timesteps_numberedit = new NumberEdit ();
                                pGrid->addEdit(P::timesteps, NumberEditType, m_timesteps_numberedit,"Timesteps [ms]:",5);
                                m_timesteps_numberedit->setEnabled(timesteps_in);

                                groupBox = new QGroupBox ("Observer Position");
                                vBox->addWidget(groupBox);
                                pGrid = new ParameterGrid<P>(this);
                                groupBox->setLayout(pGrid);
                                double distx = 1.5*outer_radius;
                                pGrid->addEdit(P::obs_x_pos_rotor, NumberEditType, new NumberEdit(),"XF:", distx);
                                pGrid->addEdit(P::obs_y_pos_rotor, NumberEditType, new NumberEdit(),"YF:", 0);
                                pGrid->addEdit(P::obs_z_pos_rotor, NumberEditType, new NumberEdit(),"ZF:", 0);

    groupBox = new QGroupBox ("Shear layer effect");
    vBox->addWidget(groupBox);
    pGrid = new ParameterGrid<P>(this);
    groupBox->setLayout(pGrid);
    m_shear_check = new QCheckBox("Shear layer:");
    pGrid->addEdit(P::shear_check, CheckBox, m_shear_check,"", shear_in);
    connect(m_shear_check,SIGNAL(toggled(bool)),this,SLOT(OnShearLayerCheck(bool)));

    m_shear_roughness_numberedit = new NumberEdit();
    pGrid->addEdit(P::shear_roughness, NumberEditType, m_shear_roughness_numberedit,"Roughness length []:", 0.01,LENGTH);
    m_shear_roughness_numberedit->setEnabled(m_shear_check->isChecked());

    m_shear_height_numberedit = new NumberEdit();
    pGrid->addEdit(P::shear_height, NumberEditType, m_shear_height_numberedit,"Measurement height []:", 10,LENGTH);
    m_shear_height_numberedit->setEnabled(m_shear_check->isChecked());
    m_shear_height_numberedit->setMinimum(0.001);

    m_shear_speed_numberedit = new NumberEdit();
    pGrid->addEdit(P::shear_speed, NumberEditType, m_shear_speed_numberedit,"Mean Wind Speed []:", 13, SPEED);
    m_shear_speed_numberedit->setEnabled(m_shear_check->isChecked());

    hBox->addLayout(vBox);
    QLabel *imageLabelb = new QLabel;
    imageLabelb->setPixmap(QPixmap(":/images/noise_3d_position_rotor.png"));
    hBox->addWidget(imageLabelb, 0, Qt::AlignHCenter);
//Sara

    setUnitContainingLabels();
    initView();
}

void NoiseCreatorDialog::initView() {
    loadObject(m_editedSimulation);

    if (m_editedSimulation) {
        m_selectFromButtons->button(m_editedSimulation->getSelectFrom())->setChecked(true);
        onSelectButtonsClicked(m_editedSimulation->getSelectFrom()); // needed because no signal if button is yet selected

        if (m_editedSimulation->getSelectFrom() == NoiseParameter::OnePolar) {
            if (!m_editedSimulation->getAnalyzedOpPoints().isEmpty()) {
                m_airfoilComboBox->setCurrentObject(
                        static_cast<CFoil*>(m_editedSimulation->getAnalyzedOpPoints()[0]->getParent()->getParent()));
                m_polarComboBox->setCurrentObject(
                        static_cast<CPolar*>(m_editedSimulation->getAnalyzedOpPoints()[0]->getParent()));
            }
        }

        QVector<OpPoint*> analyzedOpPoints = m_editedSimulation->getAnalyzedOpPoints();
        for (int i = 0; i < m_opPointRecords.size(); ++i) {
            if (analyzedOpPoints.contains(m_opPointRecords[i].opPoint)) {
                m_opPointRecords[i].checkBox->setChecked(true);
            }
        }
    } else {
        get<QLineEdit>(P::Name)->setText(g_noiseSimulationStore.getNextName("Noise Simulation"));
        m_selectFromButtons->button(NoiseParameter::OnePolar)->setChecked(true);
        onSelectButtonsClicked(NoiseParameter::OnePolar);
    }
}

void NoiseCreatorDialog::prepareOpPointRecords(bool allPolars) {
    m_opPointRecords.clear();

    for (int i = 0; i < g_oppointStore.size(); ++i) {
        OpPoint *opPoint = g_oppointStore.at(i);
        if (allPolars || static_cast<CPolar*>(opPoint->getParent()) == m_polarComboBox->currentObject()) {
            m_opPointRecords.append(OpPointRecord(opPoint->getParent()->getParent()->getName(),
                                                  opPoint->getParent()->getName(),
                                                  opPoint->getName(),
                                                  opPoint,
                                                  opPoint->m_readyForNoise));
        }
    }

    std::sort(m_opPointRecords.begin(), m_opPointRecords.end(), OpPointRecord::sort);
}

void NoiseCreatorDialog::fillOpPointView() {
    if (m_opPointViewWidget) {
        m_opPointViewWidget->deleteLater();
        m_opPointViewWidget = NULL;
    }

    if (!m_opPointRecords.isEmpty()) {
        m_opPointViewWidget = new QWidget;
        m_opPointViewWidget->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        m_opPointScrollArea->setWidget(m_opPointViewWidget);
            QGridLayout *grid = new QGridLayout;
            grid->setHorizontalSpacing(20);
            grid->setSizeConstraint(QLayout::SetMinimumSize);
            m_opPointViewWidget->setLayout(grid);
                QPushButton *button = new QPushButton ("All");
                button->setMinimumWidth(QFontMetrics(QFont()).width("All") * 1.8);
                button->setCheckable(true);
                connect(button, &QPushButton::toggled, this, &NoiseCreatorDialog::onAllButtonToggled);
                grid->addWidget(button, 0, 0, 1, 1);
                QLabel *label = new QLabel("<center>Op. Point</center>");
                grid->addWidget(label, 0, 1, 1, 1);
                label = new QLabel("<center>Polar</center>");
                grid->addWidget(label, 0, 2, 1, 1);
                label = new QLabel("<center>Airfoil</center>");
                grid->addWidget(label, 0, 3, 1, 1);
                QFont font ("Monospace");
                for (int i = 0; i < m_opPointRecords.size(); ++i) {
                    m_opPointRecords[i].checkBox = new QCheckBox;
                    if (!m_opPointRecords[i].readyForNoise) {
                        m_opPointRecords[i].checkBox->setEnabled(false);
                        m_opPointRecords[i].checkBox->setToolTip("This operational point lacks some data needed for "
                                                                 "noise simulations. Repeat the analysis in the Direct "
                                                                 "Analysis module");
                    }
                    grid->addWidget(m_opPointRecords[i].checkBox, i+1, 0, 1, 1, Qt::AlignHCenter);
                    label = new QLabel (QString("%1").arg(m_opPointRecords[i].name.toDouble(), 6, 'f', 2));
                    label->setFont(font);
                    grid->addWidget(label, i+1, 1);
                    label = new QLabel (m_opPointRecords[i].polar);
                    label->setFont(font);
                    grid->addWidget(label, i+1, 2);
                    label = new QLabel (m_opPointRecords[i].airfoil);
                    label->setFont(font);
                    grid->addWidget(label, i+1, 3);
                }
    }
}

void NoiseCreatorDialog::onSelectButtonsClicked(int id) {
    m_airfoilComboBox->setEnabled(id == NoiseParameter::OnePolar);
    m_polarComboBox->setEnabled(id == NoiseParameter::OnePolar);
    m_opPointScrollArea->setHidden(id == NoiseParameter::OriginalBpm);
    m_originalBpmWidget->setVisible(id == NoiseParameter::OriginalBpm);

    prepareOpPointRecords(id == NoiseParameter::MultiplePolars);
    fillOpPointView();
}

void NoiseCreatorDialog::onPolarBoxChange() {
    prepareOpPointRecords(false);
    fillOpPointView();
}

void NoiseCreatorDialog::onAllButtonToggled(bool pressed) {
    for (const OpPointRecord &record : m_opPointRecords) {
        if (record.checkBox->isEnabled()){//Sara
        record.checkBox->setChecked(pressed);
        }//Sara
    }
}

void NoiseCreatorDialog::onCreateButtonClicked() {
    /* check for parameter sanity */
    bool hasOpPoints = (m_selectFromButtons->checkedId() == NoiseParameter::OriginalBpm ? true : false);
    for (const OpPointRecord &record : m_opPointRecords) {
        if (record.checkBox->isChecked()) {
            hasOpPoints = true;
            break;
        }
    }
//Sara
    QString message ("");
    if (!hasOpPoints) {
//        QMessageBox::critical(this, "Create Noise Simulation",
//        "The following error(s) occured:\n\n - Simulation has no Op. Points", QMessageBox::Ok);
//        return;

        message.prepend("\n\n- Simulation has no Op. Points");
    }
//is quasi 3d?
if(sum!=2){
//        QMessageBox::information(this, "Create Noise Simulation","Select just two options to input values:\n -Rotational Speed; \n -Uniform Wind Speed; \n -TSR.",QMessageBox::Ok);
//        return;
message.prepend("\n Select just two options to input values:\n    -Rotational Speed; \n    -Uniform Wind Speed;    \n    -TSR.");

if (message != NULL){
message.prepend("The following error(s) occured:\n");
QMessageBox::critical(this, "Create Noise Simulation",message, QMessageBox::Ok);
return;}
        } else {
if (message != NULL){message.prepend("The following error(s) occured:\n");
    QMessageBox::critical(this, "Create Noise Simulation",message, QMessageBox::Ok);
return;}
//Sara

    /* create new simulation */
    NoiseSimulation *newSimulation = new NoiseSimulation (this);

    newSimulation->setSelectFrom(static_cast<NoiseParameter::OpPointSource> (m_selectFromButtons->checkedId()));

    QList<OpPoint*> analyzedOpPoints;

    for (int i = 0; i < m_opPointRecords.size(); ++i) {
        if (m_opPointRecords[i].checkBox->isChecked()) {
            analyzedOpPoints.append(m_opPointRecords[i].opPoint);
        }
    }
    newSimulation->setAnalyzedOpPoints(analyzedOpPoints.toVector());

    try {
        newSimulation->simulate();
        if (g_noiseSimulationStore.add(newSimulation)) {
            m_module->setShownSimulation(newSimulation);
            accept();  // leave dialog only if adding was successful
        }
    } catch (NoiseException &e) {
        delete newSimulation;
        QMessageBox::critical(g_mainFrame, "Simulation Error", e.what());
    }}
    onVerifyDeltaFor3D();//Sara
//Sara
}

//Sara begin
void NoiseCreatorDialog::onVerifyDeltaFor3D(){
if (check){

    QXDirect *pXDirect = (QXDirect *) g_mainFrame->m_pXDirect;

    if(pXDirect->AlphaDeltaNoise!=0){
        QMessageBox::information(this, "Step Angle Resolution!",
                              "The following error occured:\n\n - Use maximum 0.25º step angle resolution for noise prediction models in XDirect.", QMessageBox::Ok);
        return;
    }}

}

void NoiseCreatorDialog::onVerifyWindfield(){
//    QXDirect *pXDirect = (QXDirect *) g_mainFrame->m_pXDirect;

    if(g_windFieldStore.size() == 0){
        QMessageBox::critical(this, "Wind Field Error!",
                              "The following error(s) occured:\n\n - Define Windfield.", QMessageBox::Ok);
        return;
    }
}
//Sara end

void NoiseCreatorDialog::OnImportStarredD(){
    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    int number_of_elements = pBEM->dlg_elements;

QMessageBox::information (this, "δ* User Input Instructions",QString("For user input the δ* you must to follow the instructions:\n\n - Select a csv file with three columns and no header;\n\n- The first column must be filled with the index from 1 to %1;\n\n - The second column must be filled with the values of δ* on the suction side;\n\n- The third column must be filled with the values of δ* on the pressure side.").arg(number_of_elements),QMessageBox::Ok);

NoiseParameter *pNoiseParameter = (NoiseParameter *) g_mainFrame->m_pSimuWidget;

        QString PathName, strong, header;
        a_D_starred_index_user.clear();
        a_D_starred_S_user.clear();
        a_D_starred_P_user.clear();

PathName = QFileDialog::getOpenFileName(g_mainFrame, tr("Open csv File"),g_mainFrame->m_LastDirName,tr("δ* Input (*.csv)"));

        if(!PathName.length())		return ;
        int pos = PathName.lastIndexOf("/");
        if(pos>0) g_mainFrame->m_LastDirName = PathName.left(pos);

        QFile XFile(PathName);
        if (!XFile.open(QIODevice::ReadOnly))
        {
            QString strange = tr("Could not read the file\n")+PathName;
            QMessageBox::warning(g_mainFrame, tr("Warning"), strange);
            return;
        }

        QTextStream in(&XFile);

QString line;

QStringList qList;
QFile File(PathName);

if (!File.open(QIODevice::ReadOnly | QIODevice::Text))
     return;

int w=0;
//File.readLine(); // read first line and ignore
while (!File.atEnd()) {
    QString line = File.readLine(); // read wavelength line and store it
    const QStringList fields { line.split(',') };
    const QString D_starred_index { fields[0] };
    const QString D_starred_S_aux { fields[1] };
    const QString D_starred_P_aux { fields[2] };
    a_D_starred_index_user.append(D_starred_index.toInt());
    a_D_starred_S_user.append(D_starred_S_aux.toDouble());
    a_D_starred_P_user.append(D_starred_P_aux.toDouble());
    pNoiseParameter->D_starred_index_user[w]=D_starred_index.toInt();
    int u = pNoiseParameter->D_starred_index_user[w]-1.;
    pNoiseParameter->D_starred_S_user[u]=D_starred_S_aux.toDouble();
    pNoiseParameter->D_starred_P_user[u]=D_starred_P_aux.toDouble();
    ++w;
}
File.close();

if((a_D_starred_index_user.size()==pBEM->dlg_elements) & (a_D_starred_S_user.size()==pBEM->dlg_elements) & (a_D_starred_P_user.size()==pBEM->dlg_elements)){
    QMessageBox::information(g_mainFrame, tr("Import Sucessfull!"), tr("The δ* was sucessfully imported."),QMessageBox::Ok);
    return;
}
}

void NoiseCreatorDialog::OnSetDstarButton(int index){
    if (index==2){
   buttonle->setEnabled(true);}
}

//Sara
void NoiseCreatorDialog::OnModeDefine(int index){
//when unsteady state
    if (index==1){
loops_in=false;
anglesteps_in=false;
time_in=false;
shear_in=false;
timesteps_in=false;

    rotation_combobox->setCurrentIndex(0);
    m_number_loops_numberedit->setEnabled(loops_in);
    m_anglesteps_numberedit->setEnabled(anglesteps_in);
    m_time_numberedit->setEnabled(time_in);
    m_timesteps_numberedit->setEnabled(timesteps_in);
    rotation_combobox->setEnabled(false);
    m_shear_check->setEnabled(shear_in);
    onVerifyWindfield();
    }
    else {
        loops_in=true;
        anglesteps_in=true;
        time_in=false;
        shear_in=true;
        timesteps_in=false;

    rotation_combobox->setEnabled(true);
    m_number_loops_numberedit->setEnabled(loops_in);
    m_anglesteps_numberedit->setEnabled(anglesteps_in);
    m_shear_check->setEnabled(shear_in);
    m_time_numberedit->setEnabled(time_in);
    m_timesteps_numberedit->setEnabled(timesteps_in);
    }
}

void NoiseCreatorDialog::OnRotationDefine(int index){
    if (index==0){
        loops_in=true;
        anglesteps_in=true;
        time_in=false;
        timesteps_in=false;
    }
    if (index==1){
        loops_in=false;
        anglesteps_in=false;
        time_in=true;
        timesteps_in=true;
    }
    m_number_loops_numberedit->setEnabled(loops_in);
    m_anglesteps_numberedit->setEnabled(anglesteps_in);
    m_time_numberedit->setEnabled(time_in);
    m_timesteps_numberedit->setEnabled(timesteps_in);
}
//Sara

void NoiseCreatorDialog::OnRotSpeedCheck(bool index){
rot_speed_in=index;
m_rot_speed_numberedit->setEnabled(index);
}

void NoiseCreatorDialog::OnWindSpeedCheck(bool index){
u_wind_speed_in=index;
m_u_wind_speed_numberedit->setEnabled(index);
}

void NoiseCreatorDialog::OnTSRCheck(bool index){
TSR_in=index;
m_TSR_spinbox->setEnabled(index);
}

void NoiseCreatorDialog::OnShearLayerCheck(bool index){
shear_roughness_in=index;
shear_height_in=index;
shear_speed_in=index;

m_shear_roughness_numberedit->setEnabled(index);
m_shear_height_numberedit->setEnabled(index);
m_shear_speed_numberedit->setEnabled(index);
}

void NoiseCreatorDialog::OnWarningSet3(){
    sum=0;
    if(m_rot_speed_check->isChecked()){++sum;}
    if(m_u_wind_speed_check->isChecked()){++sum;}
    if(m_TSR_check->isChecked()){++sum;}
//    if(sum!=2){
//    QMessageBox::information(this, "Create Noise Simulation","Select just two options to input values:\n\n -Rotational Speed; \n -Uniform Wind Speed; \n -TSR.",QMessageBox::Ok);
//    return;
//    }
}
//Sara
