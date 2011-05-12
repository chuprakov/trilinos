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
#include <Optika_GUI.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <QFileDialog>
#include <QApplication>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <Teuchos_XMLParameterListHelpers.hpp>

  
static const std::string aboutInfo = 
"This is a generic ParameterList configurator for XML files "
"that conform to the standards of Teuchos' ParameterList serialization."
"\n\n"
"Developed by Kurtis Nusbaum and licensed under the LGPL\n\n";

const QString lastSaveDirTag("LastSaveDir");

QString getSettingsFileName(){
  QDir pldir = QDir::homePath() + "/.plconfigurator";
  if(!pldir.exists()){
    pldir.mkdir(pldir.absolutePath());
  }
  return QDir::homePath() + "/.plconfigurator/PLSettings.xml";
}


QString getLastSaveLocation(){
  QFile file(getSettingsFileName());
  QString toReturn = QDir::homePath();
  if(file.open(QIODevice::ReadOnly)){
    QXmlStreamReader xmlReader(&file);
		while(!xmlReader.isEndDocument()){
			if(xmlReader.isStartElement() && xmlReader.name() == lastSaveDirTag){
        toReturn = xmlReader.readElementText();
        break;
      }
			xmlReader.readNext();
    } 
    file.close();
  }
  return toReturn;
}

void saveLastSettings(QString lastSaveDir){
  QFile settings(getSettingsFileName());
  settings.open(QIODevice::WriteOnly);
  QXmlStreamWriter xmlWriter(&settings);
  xmlWriter.setAutoFormatting(true);
  xmlWriter.writeStartDocument();
  xmlWriter.writeStartElement("settings");
 	  xmlWriter.writeStartElement(lastSaveDirTag);
    xmlWriter.writeCharacters(lastSaveDir);
    xmlWriter.writeEndElement();
  xmlWriter.writeEndElement();
  xmlWriter.writeEndDocument();

}

void retreiveFileName(std::string& fileName){
		using namespace Qt;
    QString currentPath = QDir::currentPath();
		int argNum=1;
		char* args[1];
		std::string appName ="PLConfigurator";
		args[0] = &appName[0];
		QApplication a(argNum,args);
    QString caption("XML File to Configure");
    QString dirToStart = getLastSaveLocation();; 
    QString filter = "XML Files (*.xml)";
    QString retrieved = 
      QFileDialog::getOpenFileName(0, caption, dirToStart, filter);
    QFileInfo info(retrieved);
    saveLastSettings(info.dir().absolutePath());
    fileName = retrieved.toStdString();
}


int main(int argc, char* argv[]){
  using Teuchos::CommandLineProcessor;
  using Optika::OptikaGUI;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::DependencySheet;
  using Teuchos::getParametersFromXmlFile;

  std::string xmlFileName = "";
  CommandLineProcessor clp;
  clp.setOption(
    "xml-filename", 
    &xmlFileName,
    "The name of the xmlfile you wish to configure");

  if(xmlFileName == ""){
    retreiveFileName(xmlFileName);
    if(xmlFileName == ""){
      std::cout << "No file name provided. Exiting..." << std::endl;
      return 0;
    }
  }

  if(!QFile::exists(QString::fromStdString(xmlFileName))){
    std::cerr << "File does not exists!" << std::endl;
    std::cerr << "Bad file name: " << xmlFileName << "." << std::endl;
    return -1;
  }


  RCP<DependencySheet> deps = rcp(new DependencySheet);

  RCP<ParameterList> parameters = getParametersFromXmlFile(xmlFileName, deps);

  OptikaGUI og(parameters, deps);
  og.setAboutInfo(aboutInfo);
  og.setWindowTitle("PLEditor");
  og.setActionButtonText("Quit");
  og.setWindowIcon("plicon.svg");
  og.exec();

}
