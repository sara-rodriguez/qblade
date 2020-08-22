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
    QAction *m_exportNoise, *m_modelValidityHint, *m_exportqs3DNoiseLog, *m_exportqs3DNoise_blade, *m_exportqs3DNoise_rotor, *m_exportqs3DNoise_rotor_loops, *m_loopsReMaalpha; //Sara
	
private slots:
	void onAboutToShow ();
	void onExportNoise ();
    //Sara
    void onExportqs3DNoise_blade ();
    void onExportqs3DNoise_rotor ();
    void onExportqs3DNoiseLog  ();
    void onloopsReMaalpha ();
    //Sara
	void onModelValidityHint ();
};

#endif // NOISEMENU_H
