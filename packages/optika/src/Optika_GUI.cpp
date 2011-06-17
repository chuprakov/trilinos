// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include <QApplication>
#include <QtGui>
#include <QString>
#include "Optika_GUI.hpp"
namespace Optika{

void getInput(RCP<ParameterList> validParameters){
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		//MetaWindow *theWindow = new MetaWindow(validParameters);
		MetaWindow theWindow(validParameters);
		theWindow.show();
		a.exec();
	}
}

void getInput(RCP<ParameterList> validParameters, void (*customFunc)(RCP<const ParameterList>)){
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		//MetaWindow *theWindow = new MetaWindow(validParameters, customFunc);
		MetaWindow theWindow(validParameters, customFunc);
		theWindow.show();
		a.exec();
	}
}

void getInput(RCP<ParameterList> validParameters, RCP<DependencySheet> dependencySheet){	
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow theWindow(validParameters, dependencySheet);
		theWindow.show();
		a.exec();
	}
}

void getInput(RCP<ParameterList> validParameters, RCP<DependencySheet> dependencySheet, void (*customFunc)(RCP<const ParameterList>)){	
	{
		using namespace Qt;
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow theWindow(validParameters, dependencySheet, customFunc);
		theWindow.show();
		a.exec();
	}
}

void getInput(
  const std::string& nameOfXmlFile,
  RCP<ParameterList> userInput,
  void (*customFunc)(RCP<const ParameterList>))
{
  {
		using namespace Qt;
    RCP<DependencySheet> depSheet = rcp(new DependencySheet());
    userInput = getParametersFromXmlFile(nameOfXmlFile, depSheet);
		int argNum=1;
		char* args[1];
		std::string appName ="Optika";
		args[0] = &appName[0];
		QApplication a(argNum,args);
		MetaWindow theWindow(userInput, depSheet, customFunc);
		theWindow.show();
		a.exec();
	}
}

OptikaGUI::OptikaGUI(RCP<ParameterList> validParameters):
	validParameters(validParameters),
  customFunc(0)
  {}


OptikaGUI::OptikaGUI(RCP<ParameterList> validParameters, RCP<DependencySheet> dependencySheet):
	validParameters(validParameters),
	dependencySheet(dependencySheet),
  customFunc(0)
  {}

OptikaGUI::OptikaGUI(const std::string& xmlFileName){
  dependencySheet = rcp(new DependencySheet);
  validParameters = getParametersFromXmlFile(xmlFileName, dependencySheet);
  customFunc =0;
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
		if(is_null(dependencySheet)){
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
			QApplication::setWindowIcon(windowIcon);
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
		if(aboutInfo != ""){
			theWindow->setAboutInfo(QString::fromStdString(aboutInfo));
		}
    theWindow->setActionButtonText(QString::fromStdString(actionButtonText));
		theWindow->show();
		theWindow->activateWindow();
		a.exec();
    delete theWindow;
	}
	
}

void OptikaGUI::setWindowTitle(const std::string& title){
	this->title = title;
}

void OptikaGUI::setWindowIcon(const std::string& filePath){
	this->iconFilePath = filePath;
}

void OptikaGUI::setStyleSheet(const std::string& filePath){
	this->styleSheetFilePath = filePath;
} 

void OptikaGUI::setCustomFunction(void (*customFunc)(RCP<const ParameterList>)){
	this->customFunc = customFunc;
}

void OptikaGUI::setAboutInfo(const std::string& aboutInfo){
	this->aboutInfo = aboutInfo;
}

void OptikaGUI::setActionButtonText(const std::string& text){
  actionButtonText = text;
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

std::string  OptikaGUI::getAboutInfo(){
	return aboutInfo;
}

}

