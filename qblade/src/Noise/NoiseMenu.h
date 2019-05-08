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
    QAction *m_exportNoise, *m_modelValidityHint, *m_export3DNoise; //Sara
	
private slots:
	void onAboutToShow ();
	void onExportNoise ();
    void onExport3DNoise (); //Sara
	void onModelValidityHint ();
};

#endif // NOISEMENU_H
