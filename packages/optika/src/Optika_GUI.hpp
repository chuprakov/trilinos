// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
#ifndef OPTIKA_GUI_HPP_
#define OPTIKA_GUI_HPP_
#include "Optika_DependencySheet.hpp"
#include "Optika_metawindow.hpp"
namespace Optika{
	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters);

	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it. Also runs the function specified whenever the user hits the submit
	 * button.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 * @param customFunc Custom function to run whenever the user clicks the submit button.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>));

	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet);

	/**
	 * Retreives the input for a Teuchos Parameter List using a GUI. Note the Parameter List will be edited.
	 * All user input will be stored in it. Also runs the function specified whenever the user hits the submit
	 * button.
	 *
	 * @param validParameters A list of parameters from which the users may specify values.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 * @param customFunc Custom function to run whenever the user clicks the submit button.
	 */
	void getInput(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet, void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>));

class OptikaGUI{
public:
	OptikaGUI(Teuchos::RCP<Teuchos::ParameterList> validParameters);
	OptikaGUI(Teuchos::RCP<Teuchos::ParameterList> validParameters, Teuchos::RCP<DependencySheet> dependencySheet);

	void exec();

	void setWindowTitle(std::string title);
	void setWindowIcon(std::string filePath);
	void setStyleSheet(std::string filePath);
	void setCustomFunction(void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>));

	std::string getWindowTitle();
	std::string getWindowIcon();
	std::string getStyleSheet();

private:
	Teuchos::RCP<Teuchos::ParameterList> validParameters;
	Teuchos::RCP<DependencySheet> dependencySheet;
	std::string title, iconFilePath, styleSheetFilePath;
	void (*customFunc)(Teuchos::RCP<const Teuchos::ParameterList>);
};

}


#endif //OPTIKA_GUI_HPP_
