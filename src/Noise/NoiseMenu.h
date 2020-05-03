#ifndef NOISEMENU_H
#define NOISEMENU_H

#include <QMenu>
class QMainWindow;

class NoiseModule;


class NoiseMenu : public QMenu
{
	Q_OBJECT
	
public:
	NoiseMenu (QMainWindow *parent, NoiseModule *module);
	
private:
	NoiseModule *m_module;
    QAction *m_exportNoise, *m_modelValidityHint, *m_exportqs3DNoiseLog, *m_exportqs3DNoise_blade, *m_exportqs3DNoise_rotor, *m_exportqs3DNoise_rotor_loops; //Sara
	
private slots:
	void onAboutToShow ();
	void onExportNoise ();
    void onExportqs3DNoise_blade ();//Sara
    void onExportqs3DNoise_rotor ();//Sara
    void onExportqs3DNoise_rotor_loops ();//Sara
    void onExportqs3DNoiseLog  (); //Sara
	void onModelValidityHint ();
};

#endif // NOISEMENU_H
