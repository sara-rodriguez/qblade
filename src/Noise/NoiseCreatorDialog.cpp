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
#include "../XDirect/FoilPolarDlg.h"
//Sara

typedef Parameter::NoiseSimulation P;

NoiseCreatorDialog::NoiseCreatorDialog(NoiseSimulation *presetSimulation, NoiseModule *module)
	: m_module(module),
	  m_editedSimulation(presetSimulation),
	  m_opPointViewWidget(NULL)
{
    setWindowTitle("Noise Simulation");	//Sara
	
	QTabWidget *tabWidget = new QTabWidget;
	m_contentVBox->insertWidget(0, tabWidget);
	
	QWidget *widget = new QWidget;
	tabWidget->addTab(widget, "Parameters");
		QHBoxLayout *hBox = new QHBoxLayout;
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
								  "Original flow Mach Number (M):", 0.21);
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
								  "Directivity angle θe [deg]:", 90);
                    pGrid->addEdit(P::DirectivityPhi, NumberEditType, new NumberEdit(),
								  "Directivity angle ψe [deg]:", 90);

					
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
                pGrid = new ParameterGrid<P>(this);//Sara urgente
                groupBox->setLayout(pGrid);
QComboBox *Lowson_type_combobox = new QComboBox;
pGrid->addEdit(P::Lowson_type,ComboBox, Lowson_type_combobox,"Lowson's Model:","");
Lowson_type_combobox->insertItem(0,"None");
Lowson_type_combobox->insertItem(1,"Von Kármán");
Lowson_type_combobox->insertItem(2,"Rapid Distortion");
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
							pGrid->addEdit(P::Aoa, NumberEditType, new NumberEdit, "AOA (α) [deg]:", 0);
							pGrid->addEdit(P::ChordBasedReynolds, NumberEditType, new NumberEdit,
										   "Chord based Reynolds number (Rc):", 100000);
							pGrid->addComboBox(P::Transition, "Type of Transition:", NoiseParameter::TransitionFlow,
											   QStringList()<<"Fully turbulent"<<"Transition flow");

                            //Sara
                            widget = new QWidget;
                            tabWidget->addTab(widget, "3D Analysis");
                            hBox = new QHBoxLayout;
                            widget->setLayout(hBox);
                            groupBox = new QGroupBox ("3D Simulation Parameters");
                            hBox->addWidget(groupBox);
                            pGrid = new ParameterGrid<P>(this);
                            groupBox->setLayout(pGrid);
                            NumberEdit *numEdit = new NumberEdit();
                            QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
                            double max_sects = pBEM->dlg_elements;
                            pGrid->addEdit(P::sects, NumberEditType, numEdit,"Number of Segments:", 40);
                            numEdit->setRange(13,max_sects);
                            numEdit->setValue(30,true);
                            numEdit->setAutomaticPrecision(0);

m_rot_speed_check = new QCheckBox("rot. speed set:");
pGrid->addEdit(P::rot_speed_check, CheckBox, m_rot_speed_check,"", 0);
    connect(m_rot_speed_check,SIGNAL(clicked()),this,SLOT(OnRotSpeedCheck()));//int

m_rot_speed_numberedit = new NumberEdit ();
m_rot_speed_numberedit->setAutomaticPrecision(3);
pGrid->addEdit(P::rot_speed, NumberEditType, m_rot_speed_numberedit,"Rotational Speed [rpm]:",1);
if(m_rot_speed_check->isChecked()){m_rot_speed_numberedit->setEnabled(true);}else{m_rot_speed_numberedit->setEnabled(false);}

m_u_wind_speed_check = new QCheckBox("wind speed set:");
pGrid->addEdit(P::u_wind_speed_check, CheckBox, m_u_wind_speed_check,"", 1);
connect(m_u_wind_speed_check,SIGNAL(clicked()),this,SLOT(OnWindSpeedCheck()));//int stateChanged(int)


m_u_wind_speed_numberedit = new NumberEdit ();
m_u_wind_speed_numberedit->setAutomaticPrecision(3);
pGrid->addEdit(P::u_wind_speed, NumberEditType, m_u_wind_speed_numberedit,"Uniform Wind Speed []::",u_wind_speed, SPEED);
if(m_u_wind_speed_check->isChecked()){m_u_wind_speed_numberedit->setEnabled(true);}else{m_u_wind_speed_numberedit->setEnabled(false);}

m_TSR_spinbox = new QDoubleSpinBox;
m_TSR_spinbox->setLocale(QLocale("en_us"));

m_TSR_check = new QCheckBox ("TSR set");
//const bool blocked = m_TSR_check->signalsBlocked();
//m_TSR_check->blockSignals(true);
pGrid->addEdit(P::TSR_check, CheckBox, m_TSR_check,"", 1);
//m_TSR_check->blockSignals(blocked);

connect(m_TSR_check,SIGNAL(clicked()),this,SLOT(OnTSRCheck()));

pGrid->addEdit(P::TSRtd,DoubleSpinBox, m_TSR_spinbox,"TSR:", 7);

double lstart  =   pSimuWidget->m_pctrlLSLineEdit->getValue();
double ldelta  =   pSimuWidget->m_pctrlLDLineEdit->getValue();
double lend  =   pSimuWidget->m_pctrlLELineEdit->getValue();

     m_TSR_spinbox->setRange(lstart, lend);
     m_TSR_spinbox->setSingleStep(ldelta);
     m_TSR_spinbox->setDecimals(1);
     m_TSR_spinbox->valueChanged(change_TSR);

QLabel *labeltd = new QLabel ("Select Blade Type from Database:");
pGrid->addWidget(labeltd);
m_airfoilComboBoxtd = new FoilComboBox (&g_foilStore);
pGrid->addWidget(m_airfoilComboBoxtd);

dstar_combobox = new QComboBox;
pGrid->addEdit(P::dstar_type,ComboBox, dstar_combobox,"δ* type:","");
dstar_combobox->insertItem(0,"BPM");
dstar_combobox->insertItem(1,"XFoil");
dstar_combobox->insertItem(2,"User");
connect(dstar_combobox,SIGNAL(currentIndexChanged(int)),this,SLOT(OnSetDstarButton(int)));

pGrid->setSizeConstraint(QLayout::SetMinimumSize);
buttonle = new QPushButton ("δ* User Input");
buttonle->setMinimumWidth(QFontMetrics(QFont()).width("δ* User Input") * 1.8);
buttonle->setEnabled(false);
pGrid->addWidget(buttonle,9,2);//, 0, 0, 1, 1 , 7, 5, 1, 1
connect(buttonle,SIGNAL(clicked()),this,SLOT(OnImportStarredD()));

QComboBox *phi_combobox = new QComboBox;
pGrid->addEdit(P::phi_type,ComboBox, phi_combobox,"Φ type:","");
phi_combobox->insertItem(0,"90º");
phi_combobox->insertItem(1,"free");

//groupBox = new QGroupBox ("Observer Position:");
//hBox->addWidget(groupBox);
//pGrid = new ParameterGrid<P>(this);
//groupBox->setLayout(pGrid);

QLabel *labelte = new QLabel ("Observer Position:");
pGrid->addWidget(labelte,13,0);
pGrid->addEdit(P::obs_x_pos, NumberEditType, new NumberEdit(),"X:", 10);
pGrid->addEdit(P::obs_y_pos, NumberEditType, new NumberEdit(),"Y:", 10);

QBEM *pbem = (QBEM *) g_mainFrame->m_pBEM;
double hub_radius=pbem->m_pBlade->m_HubRadius;
double outer_radius=pbem->m_pTData->OuterRadius;
double blade_radius=(outer_radius-hub_radius);
double z_pos=blade_radius/2.;

pGrid->addEdit(P::obs_z_pos, NumberEditType, new NumberEdit(),"Z:", z_pos);
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
		record.checkBox->setChecked(pressed);
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
    if (!hasOpPoints) {
        QMessageBox::critical(this, "Create Noise Simulation",
                              "The following error(s) occured:\n\n - Simulation has no Op. Points", QMessageBox::Ok);
        return;
    }
    //Sara todo criar alertas para usuario
//    if (!hasOpPoints) {
//		QMessageBox::critical(this, "Create Noise Simulation",
//							  "The following error(s) occured:\n\n - Reynolds Lower 600,000", QMessageBox::Ok);
//		return;
//	}
//    if (!hasOpPoints) {
//		QMessageBox::critical(this, "Create Noise Simulation",
//							  "Alert:\n\n - Reynolds Upper 2,400,000", QMessageBox::Ok);
//		return;
//	}
//    if (!hasOpPoints) {
//		QMessageBox::critical(this, "Create Noise Simulation",
//							  "Alert:\n\n - Mach Upper 0.21", QMessageBox::Ok);
//		return;
//	}
//    if (!hasOpPoints) {
//		QMessageBox::critical(this, "Create Noise Simulation",
//							  "Alert:\n\n - AOA(abs) Upper 19.8''", QMessageBox::Ok);
//		return;
//	}
    //Sara
	
	/* create new simlation */
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
	}
}

void NoiseCreatorDialog::OnImportStarredD(){
    QBEM *pBEM = (QBEM *) g_mainFrame->m_pBEM;
    int number_of_elements = pBEM->dlg_elements;

QMessageBox::information (this, "δ* User Input Instructions",QString("For user input the δ* you must to follow the instructions:\n\n - Select a csv file with tho columns and no header;\n\n- The first column must be filled with δ* on  suction side;\n\n- The second column must be filled with δ* on pressure side;\n\n -The csv file must have %1 rows filled.").arg(number_of_elements),QMessageBox::Ok);

NoiseParameter *pNoiseParameter = (NoiseParameter *) g_mainFrame->m_pSimuWidget;

        QString PathName, strong, header;
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
    const QString D_starred_S_aux { fields[0] };
    const QString D_starred_P_aux { fields[1] };
    a_D_starred_S_user.append(D_starred_S_aux.toDouble());
    a_D_starred_P_user.append(D_starred_P_aux.toDouble());
    pNoiseParameter->D_starred_S_user[w]=D_starred_S_aux.toDouble();
    pNoiseParameter->D_starred_P_user[w]=D_starred_P_aux.toDouble();
    ++w;
}
File.close();

if((a_D_starred_S_user.size()==pBEM->dlg_elements) & (a_D_starred_P_user.size()==pBEM->dlg_elements)){
    QMessageBox::information(g_mainFrame, tr("Import Sucessfull!"), tr("The δ* was sucessfully imported."),QMessageBox::Ok);
    return;
}
}

void NoiseCreatorDialog::OnSetDstarButton(int index){
    if (index==2){
   buttonle->setEnabled(true);}
}

void NoiseCreatorDialog::OnRotSpeedCheck(){
    if (!m_rot_speed_check->isChecked()){m_rot_speed_numberedit->setEnabled(false);}
    else{m_rot_speed_numberedit->setEnabled(true);}
    OnWarningSet3();
}

void NoiseCreatorDialog::OnWindSpeedCheck(){
    if (!m_u_wind_speed_check->isChecked()){m_u_wind_speed_numberedit->setEnabled(false);}
    else{m_u_wind_speed_numberedit->setEnabled(true);}
    OnWarningSet3();
}

void NoiseCreatorDialog::OnTSRCheck(){
    if (!m_TSR_check->isChecked()){m_TSR_spinbox->setEnabled(false);}
    else{m_TSR_spinbox->setEnabled(true);}
        OnWarningSet3();
}

void NoiseCreatorDialog::OnWarningSet3(){
    int sum=0;
    if(m_rot_speed_check->isChecked()){++sum;}
    if(m_u_wind_speed_check->isChecked()){++sum;}
    if(m_TSR_check->isChecked()){++sum;}
    if(sum!=2){
    QMessageBox::information(this, "Create Noise Simulation","Select 2 options to input values!",QMessageBox::Ok);
    return;
    }
}
//Sara
