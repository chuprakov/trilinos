#include <QApplication>
#include <QtGui>
#include "Optika_GUI.hpp"
#include <iostream>
namespace Optika{


void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters){
	{
		using namespace Qt;
		QApplication a(0,0);
		MetaWindow *theWindow = new MetaWindow(validParameters);
		theWindow->show();
		a.exec();
	}
}

void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet){	
	{
		using namespace Qt;
		QApplication a(0,0);
		MetaWindow *theWindow = new MetaWindow(validParameters, dependencySheet);
		theWindow->show();
		a.exec();
	}
}

}

