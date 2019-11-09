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
    QAction *m_exportNoise, *m_modelValidityHint, *m_exportqs3DNoiseLog, *m_exportqs3DNoise_final; //Sara *m_exportqs3DNoise, *m_exportqs3DNoiseComplete,
	
private slots:
	void onAboutToShow ();
	void onExportNoise ();
    void onExportqs3DNoise_final ();//Sara
    void onExportqs3DNoiseLog  (); //Sara
    void onExportqs3DNoiseComplete (); //Sara
    void onExportqs3DNoise (); //Sara
	void onModelValidityHint ();
};

#endif // NOISEMENU_H
