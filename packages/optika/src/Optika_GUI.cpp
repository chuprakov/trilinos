#include <QApplication>
#include <QtGui>
#include "Optika_GUI.hpp"
#include <iostream>
#include <pthread.h>
namespace Optika{

void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters){
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow *theWindow = new MetaWindow(validParameters);
		theWindow->show();
		a.exec();
	}
}

void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>)){
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow *theWindow = new MetaWindow(validParameters, customFunc);
		theWindow->show();
		a.exec();
	}
}

void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet){	
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow *theWindow = new MetaWindow(validParameters, dependencySheet);
		theWindow->show();
		a.exec();
	}
}

void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet, void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>)){	
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow *theWindow = new MetaWindow(validParameters, dependencySheet, customFunc);
		theWindow->show();
		a.exec();
	}
}


}

