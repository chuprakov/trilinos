#include <QApplication>
#include <QtGui>
#include "Optika_GUI.hpp"
#include <iostream>
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

OptikaGUI::OptikaGUI(Teuchos::RCP<Teuchos::ParameterList> validParameters){
	this->validParameters = validParameters;
}


OptikaGUI::OptikaGUI(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet){
	this->validParameters = validParameters;
	this->dependencySheet = dependencySheet;
}

void OptikaGUI::exec(){
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow *theWindow;
		if(Teuchos::is_null(dependencySheet)){
			theWindow = new MetaWindow(validParameters, customFunc);
		}
		else{
			theWindow = new MetaWindow(validParameters, dependencySheet, customFunc);
		}
		if(title != ""){
			theWindow->setWindowTitle(QString::fromStdString(title));
		}
		if(iconFilePath != ""){
			QIcon windowIcon(QString::fromStdString(iconFilePath));
			theWindow->setWindowIcon(windowIcon);
		}
		if(styleSheetFilePath != ""){
			QString str;
			QFile file(QString::fromStdString(styleSheetFilePath));
			if (file.open(QIODevice::ReadOnly | QIODevice::Text)){
				QTextStream in(&file);
				while (!in.atEnd()) {
					str += in.readLine();
				}
				a.setStyleSheet(str);
			}
		}
		theWindow->show();
		a.exec();
	}
	
}

void OptikaGUI::setWindowTitle(std::string title){
	this->title = title;
}

void OptikaGUI::setWindowIcon(std::string filePath){
	this->iconFilePath = filePath;
}

void OptikaGUI::setStyleSheet(std::string filePath){
	this->styleSheetFilePath = filePath;
} 
void OptikaGUI::setCustomFunction(void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>)){
	this->customFunc = customFunc;
}

std::string OptikaGUI::getWindowTitle(){
	return title;
}

std::string OptikaGUI::getWindowIcon(){
	return iconFilePath;
}

std::string OptikaGUI::getStyleSheet(){
	return styleSheetFilePath;
}

}

