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
