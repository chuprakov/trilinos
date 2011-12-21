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
#ifndef OPTIKA_MODEL_HPP_
#define OPTIKA_MODEL_HPP_
#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>
#include <QDir>
#include "Optika_treeitem.hpp"

/*! \file Optika_TreeModel.hpp
    \brief The model used by Optika
    in its implementation of the MVC
    framework.
*/

class QDomElement;

namespace Optika{

class TreeItem;

/**
 * \brief TreeModel is a type of QAbstractItemModel that has a Tree like structure.
 *
 * Note: For all undocumented functions, please refer to the Qt api. They will have a good desciption.
 */
class TreeModel : public QAbstractItemModel{
	Q_OBJECT
public:

  /** \name Constructors/Destructor */
  //@{

	/**
	 * \brief Constructs the TreeModel.
	 * 
	 * @param validParameters A list of parameters for which the users must enter values.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 * @param saveFileName Name of a save file used in a previous attempt to get values for the validParameters ParameterList.
	 * @param parent The parent object.
	 */
	TreeModel(
    RCP<ParameterList> validParameters,
    RCP<DependencySheet> dependencySheet=null,
	  QString saveFileName=QString(), 
    QObject *parent=0);

	/**
	 *
	 * \brief Deconstructor for the TreeModel.
	 */
	~TreeModel();

  /** \name Overridden from QAbstractItemModel */
  //@{

  /** * \brief .  */
	QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
  /** * \brief .  */
	Qt::ItemFlags flags(const QModelIndex &index) const;
  /** * \brief .  */
	QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;
  /** * \brief .  */
	QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const;
  /** * \brief .  */
	QModelIndex parent(const QModelIndex &index) const;
  /** * \brief .  */
	bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole);
  /** * \brief .  */
	int rowCount(const QModelIndex &parent = QModelIndex()) const;
  /** * \brief .  */
	int columnCount(const QModelIndex &parent = QModelIndex()) const;

  //! @name Helper Functions
  //@{

	/**
   * \brief Issues any signals that need to be emitted right away.
   *
	 * If this TreeModel has a dependent Parameter List, then all the depndencies need to be evaluated before the Parameter List may be displayed.
	 * Certain items might need to be hidden before the user even starts entering data. 
	 * This function goes through all of the depndees in the Dependent Parameter List and issues a signal saying they've changed.
	 * They really haven't changed yet, but this allows all the depndencies to be evaluated and any initial visual settings to be
	 * displayed correctly.
	 */
	void issueInitilizationSignals();

  //@}

  //! @name Debug Functions
  //@{
  
	/**
	 * \brief Prints out the model.
	 */
	void printOut() const;

  //@}

  //! @name Input/Output Functions
  //@{

	/**
	 * \brief Writes out the state of the current parameters in xml format.
	 *
	 * @param fileName The name of the file to which the TreeModel should write the XML output.
	 */
	bool writeOutput(QString fileName);

	/**
	 * \brief Reads an xml file that describes the state of current parameters in xml format.
	 *
	 * @param fileName The name of the file from which the TreeModel should read parameter values.
	 */
	void readInput(QString fileName);

  //@}

  //! @name Getters and Setters
  //@{
  
	/**
	 * \brief Gets the name of the save file with which the TreeModel is associated.
   *
	 * If the TreeModel has yet to be saved and thus has no save file associated with it, the funtion will return an empty string.
	 *
	 * @return The name of the save file with which the TreeModel is associated.
	 */
	QString getSaveFileName();

	/**
	 * \brief Determines wether or not the current state of TreeModel has been saved.
	 *
	 * @return True if the current state of the TreeModel has been saved. False otherwise.
	 */
	bool isSaved();

	/**
	 * \brief Resets all the inputs to their default values.
	 */
	void reset();

	/**
	 * \brief Returns the type of item located at the specified QModelIndex.
	 *
	 * @param index The index of the TreeItem.
	 * @return The type of the item at the index.
	 */
	QString itemType(const QModelIndex &index) const;

	/**
	 * \brief Determines whether or not a Dependent Parameter List is being used in the TreeModel.
	 *
	 * @return True if the TreeModel has dependencies, false otherwise.
	 */
	bool hasDependencies();

	/**
	 * \brief Determines whether or not the value at the valueToCheck 
	 * is valid.
	 *
	 * @param valueToCheck The index of the item whose valididty
	 * is in questions.
	 * @return True if the value at index is valid, false otherwise.
	 */
	bool hasValidValue(QModelIndex valueToCheck) const;

	/**
	 * \brief Gets the validator for a particular TreeItem.
	 *
	 * @param index The index of the TreeItem whose validators is sought.
	 * @return The validator at the given index.
	 */
	RCP<const ParameterEntryValidator> getValidator(const QModelIndex &index) const;

	/**
	 * \brief Gets the array for a particular TreeItem.
	 *
	 * @param index The index of the TreeItem whose array is sought.
	 * @return The array at the given index.
	 */
	template <class S>
	Array<S> getArray(const QModelIndex& index){
		return any_cast<Array<S> >(itemEntry(index)->getAny()); 
	}

	/**
	 * \brief Gets the TwoDArray for a particular TreeItem.
	 *
	 * @param index The index of the TreeItem whose TwoDArray is sought.
	 * @return The TwoDArray at the given index.
	 */
	template <class S>
	TwoDArray<S> getTwoDArray(const QModelIndex& index){
		return any_cast<TwoDArray<S> >(itemEntry(index)->getAny()); 
	}

	/**
	 * \brief Get a ParameterList containing all of the parameters at their current settings.
	 *
	 * @return A ParameterList containing all of the parameters at their current settings.
	 */
	RCP<const ParameterList> getCurrentParameters();

	/**
	 * \brief Finds the index of a particular parameter entry.
	 *
	 * @param parameterEntry The ParameterEntry whose index is being sought.
	 */
	QModelIndex findParameterEntryIndex(RCP<const ParameterEntry> parameterEntry);

  //@}

  //! @name Constant Getting Functions.
  //@{

  /**
   * \brief Returns constant representing the RawDataRole
   *
   * @return The constant representing the RawDataRole.
   */
  static const int getRawDataRole(){
    static const int RawDataRole = Qt::UserRole;
    return RawDataRole;
  }

  //@}

signals:

  //! @name Public Signals
  //@{

	/**
	 * \brief Emitted when a row should be hidden.
	 *
	 * @param row The row of the item that should be hidden.
	 * @param parent The parent of the item that should be hidden.
	 */
	void hideData(int row, const QModelIndex& parent);

	/**
	 * \brief Emitted when a row should be shown.
	 *
	 * @param row The row of the item that should be shown.
	 * @param parent The parent of the item that should be shown.
	 */
	void showData(int row, const QModelIndex& parent);

	/**
	 * \brief Emitted when it has been determined that a TreeItem no longer has a valid value.
	 *
	 * @param badItem The index of the item that now has a bad value.
	 * @param message A message describing what happened to cause the
	 * item to obtain an invalid value.
	 */
	void badValue(QModelIndex badItem,  QString message);

  //@}

private:
  /** \name Private Members */
  //@{
  
	/**
	 * \brief Whether or not the model has been saved since it was last modified.
	 */
	bool saved;

	/**
	 * \brief Whether or not the model has any dependencies.
	 */
	bool dependencies;

	/**
	 * \brief The name of the savefile associated with the model.
	 */
	QString saveFileName;

	/**
	 * \brief The root item of the model.
	 */
	TreeItem *rootItem;

	/**
	 * \brief The list of valid parameters.
	 */
	RCP<ParameterList> validParameters;

	/**
	 * \brief A canonical list of what the validParameters were when they were first
	 * passed to the treemodel. Used by the reset function.
	 */
	RCP<const ParameterList> canonicalList;

	/**
	 * \brief The dependency sheet being used to determine any
	 * depdendencies between parameters.
	 */
	RCP<DependencySheet> dependencySheet;

  //@}

  /** \name Private Functions */
  //@{

	/**
	 * \brief Gets the ParameterEntry object given a QModelIndex.
	 *
	 * @param index Index of the TreeItem for which the ParameterEntry is desired.
	 * @return The ParameterEntry associated with the QModelIndex.
	 */
	RCP<const ParameterEntry> 
    itemEntry(const QModelIndex &index) const;

	/**
	 * \brief Reads in the parameter list to be represented by the model.
	 * 
	 * @param validParameters The list to be read.
	 * @param parentItem The initial parent tree item to be used.
	 */
	void readInParameterList(RCP<ParameterList> validParameters, TreeItem *parentItem);

	/**
	 * \brief Inserts a new parameter list into the model.
	 *
	 * @param parameterList The ParameterList to be inserted.
	 * @param listEntry The ParameterEntry of the ParameterList to be inserted.
	 * @param plname The name of the ParameterList.
	 * @param The parent TreeItem.
	 */
	void insertParameterList(RCP<ParameterList> parameterList, RCP<ParameterEntry> listEntry, std::string plname, TreeItem *parent);

	/**
	 * \brief Inserts a new parameter into the model.
	 *
	 * @param listEntry The ParameterEntry of the Parameter to be inserted.
	 * @param name The name of the Parameter.
	 * @param The parent TreeItem.
	 */
	void insertParameter(RCP<ParameterEntry> parameter, std::string name, TreeItem *parent);

	/**
	 * \brief Basic setup shared by each of the constructors
	 *
	 * @param saveFileName The saveFileName parameter passed to the constructors.
	 */
	void basicSetup(QString saveFileName);

	/**
	 * \brief Checks the state of a dependent after it's dependency has been evaluated. 
   *
   * Takes appropriate action if any more modifications to the model need to be made or if
	 * the view needs to know anything as a result of the change.
	 */
	void checkDependentState(const QModelIndex dependee, RCP<Dependency> dependency);

  /**
   *
   * \brief Finds the QModelIndex associated with a parameter entry.
   *
   * @param start The index where we should start looking.
   * @param parameterEntry The parameter entry we're looking for.
   * @return The index associated with the parameter entry.
   */
  QModelIndex parameterEntryMatch(const QModelIndex &start,
    const RCP<const ParameterEntry> &parameterEntry) const;


  /**
   * \brief Given a Dom element, searches for the corresponding parameter
   * in the model, updates it's value with the value from the Dom element,
   * and then recusively does the same for all children.
   *
   * @param element The element for which the corresponding parameter
   * in the model and it's children should be updated.
   */
  void processInputElement(const QDomElement& element);

  /**
   * \brief Determines whether or not a model index corresponds to the
   * parameter represented by the DomElement.
   *
   * This funciton determines whether or not the Dom Element and the 
   * model index actually represent the same parameter by verifying they have
   * the same set of parent nodes.
   */
  bool isRealMatch(
    const QDomElement& element, 
    const QModelIndex& potentialMatch) const;

  /**
   * Determines whether or not the given index is the root index.
   *
   * @param The index in question
   * @return True is the index is the root index, false otherwise.
   */
  bool isRootIndex(const QModelIndex& index) const;
  

  //@}

private slots:

  //! @name Private Slots
  //@{

	/**
	 * \brief When the state of any of the MainTree's items is changed, this slot should be called
	 */
	void currentFileNowModified();

	/**
	 * \brief Listens to see if any data has changed. 
   *
   * If so and the item has dependencies, this function will make sure all appropriate signals are emitted,
	 * and any changes that need to be made to the model are made.
	 * 
	 * @param index1 The start index of the data that changed.
	 * @param index2 The end index of the data that changed.
	 */
	void dataChangedListener(const QModelIndex& index1, const QModelIndex& index2);

  //@}
};


} //end namespace 
#endif /* OPTIKA_MODEL_HPP_ */
