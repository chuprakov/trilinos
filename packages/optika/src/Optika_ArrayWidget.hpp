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
#ifndef OPTIKA_ARRAYWIDGET_HPP_
#define OPTIKA_ARRAYWIDGET_HPP_

#include <QDialog>
#include <QModelIndex>
#include <QPushButton>
#include <QGridLayout>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QLineEdit>
#include <QGridLayout>
#include <QScrollArea>
#include <QLabel>
#include <vector>
#include "Optika_treemodel.hpp"
#include "Optika_FileNameWidget.hpp"
#include "Optika_SpinBoxApplier.hpp"

/*! \file Optika_ArrayWidget.hpp
    \brief A collection of Widgets used to edit
    Arrays in  Optika.
*/

namespace Optika {

template<class S>
class GenericArrayWidget : public QDialog{

public:

	GenericArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0);

	/**
	 * Gets the type of array being edited.
	 *
	 * @return The type of array being edited.
	 */
	const QString getType() const {
		return type;
	}

  const QString getName() const{
    return name;
  }

  const RCP<const ParameterEntryValidator> getEntryValidator() const{
    return entryValidator;
  }

  /**
   * \brief called when the user is done entering data
   * into the widget. MUST BE IMPLEMENTED AS A SLOT IN
   * SUBCLASSES!
   */
  virtual void accept() =0;


protected: 

	/**
	 * \brief Sets up the layout for the arrayContainer, including adding what ever editing
	 * widget should be used for the particual type of array.
	 */
	virtual void setupArrayLayout(){
    if(arrayContainer->layout() == NULL){
      arrayContainer->setLayout(getArrayLayout());
    }
  }

  /** @name Protected Members */
  //@{

	/**
	 * \brief The widget containing all of the editing widgets (e.g.
	 * QLineEdits, and QSpinBoxes) that comprise the array editor.
	 */
	QWidget *arrayContainer;


  virtual QLayout* getArrayLayout() =0;  

	/**
   * \brief Gathers all the user inputed data and closes the dialog.
	 */
	virtual void doAcceptWork() =0;


  //@}
	

private:


	/**
	 * \brief The type of array.
	 */
	QString type;

  /**
   * The name of the Parameter being edited.
   */
  QString name;

	/**
	 * \brief The validator being used on the array.
	 */
	RCP<const ParameterEntryValidator> entryValidator;	

};

template<class S>
GenericArrayWidget<S>::GenericArrayWidget(
  QString name, 
  QString type, 
  const RCP<const ParameterEntryValidator> validator,
  QWidget *parent):
  QDialog(parent),
  type(type),
  name(name),
  entryValidator(validator)
{
	setModal(true);
	setSizeGripEnabled(true);
	arrayContainer = new QWidget(this);

	QScrollArea *scrollArea = new QScrollArea(this);
	scrollArea->setWidget(arrayContainer);
	scrollArea->setWidgetResizable(true);

	QPushButton *doneButton = new QPushButton(tr("Done"));
	QPushButton *cancelButton = new QPushButton(tr("Cancel"));
	connect(doneButton, SIGNAL(clicked(bool)), this, SLOT(accept()));
	connect(cancelButton, SIGNAL(clicked(bool)), this, SLOT(reject()));
	QGridLayout *layout = new QGridLayout(this);
	layout->addWidget(scrollArea,0,0,1,3);
	layout->addWidget(doneButton,1,2);
	layout->addWidget(cancelButton,1,1);

	this->setLayout(layout);

	setWindowTitle(name);
}


template <class S>
class Generic2DArrayWidget : public GenericArrayWidget<S>{

	Generic2DArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0);

  void initData(TwoDArray<S> array){
    baseArray = array;
    this->setupArrayLayout();
  }

  void doAcceptWork(){
    baseArray.clear();
    baseArray = getArrayFromWidgets();
    this->done(QDialog::Accepted);
  }


protected:

  
  virtual QWidget* getEditorWidget(int row, int col) =0;

  virtual TwoDArray<S> getArrayFromWidgets()=0;

  TwoDArray<QWidget*> widgetArray;
private:
  TwoDArray<S> baseArray;
	QLayout* getArrayLayout();

};

template<class S>
Generic2DArrayWidget<S>::Generic2DArrayWidget(
  QString name, 
  QString type, 
  const RCP<const ParameterEntryValidator> validator,
  QWidget *parent):
  GenericArrayWidget<S>(name, type, validator, parent)
{}


template<class S>
QLayout* Generic2DArrayWidget<S>::getArrayLayout(){
 QGridLayout *widgetLayout = new QGridLayout;
  for(int i =0; i < baseArray.getNumColumns(); ++i){
		widgetLayout->addWidget(new QLabel("Column: " +QString::number(i)),0,i+1,Qt::AlignLeft);
  }
  for(int i =0; i < baseArray.getNumRows(); ++i){
		widgetLayout->addWidget(new QLabel("Row: " +QString::number(i)),i+1,0,Qt::AlignLeft);
  }
  for(int i =0; i < baseArray.getNumRows(); ++i){
    for(int j =0; j < baseArray.getNumRows(); ++j){
		  QWidget* editorWidget = getEditorWidget(i,j);
		  widgetLayout->addWidget(editorWidget,i+1,j+1,Qt::AlignLeft);
		  widgetArray(i,j) = editorWidget;
    }
  }
  return widgetLayout;
}


/**
 * \brief A templated abstract base class for all other array editing widgets. 
 *
 * Note the absence of the Q_OBJECT
 * macro. This is becuase classes using the Q_OBJECT macro can't be templated (bummer). The macro is therfore
 * present in the subclasses.
 *
 * \reference Array
 */ 
template <class S>
class Generic1DArrayWidget : public GenericArrayWidget<S>{
public:

  /** \name Constructor */
  //@{

  /**
   * \brief Constructs a Generic1DArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	Generic1DArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0);

  //@}
	
  //! @name Attribute/Query Methods 
  //@{


  /**
   * Return the array backing this widget.
   *
   * @param The array backing this widget.
   */
  const Array<S> getData() const{
    return baseArray;
  }

  //@}


  /** @name Miscellaneous */
  //@{

  /**
   * \brief Initialize the data in the widget with the currecnt values of the baseArray as
   * well as setting up the widgets layout.
   */
  void initData(Array<S> array){
    baseArray = array;
    this->setupArrayLayout();
  }

	/**
	 * \brief Gets the widget to be used as an editor for each entry in the array.
	 */
	virtual QWidget* getEditorWidget(int index) = 0;


  //@}

protected:

  /** \name Protected types */
  //@{

	/**
	 * \brief Convienece typedef. Represents an array of QWidgets.
	 */
	typedef std::vector<QWidget*> WVector;

  //@}

  /** \name Protected members */
  //@{

	/**
	 * \brief Conatins the editing widgets (e.g. QLineEdits and QSpinBoxes) comprising the array editor.
	 */
	WVector widgetVector;

	/**
	 * The array to be edited.
	 */
	Array<S> baseArray;

  void doAcceptWork();
  //@}


private:

  /** \name Private Functions */
  //@{

  /**
   * \brief Get a new array reflecting the current values entered in the widgets.
   *
   * @return A new array reflecting the currecnt values entered in the widgets.
   */
  virtual Array<S> getArrayFromWidgets() = 0;

  QLayout* getArrayLayout();

  //@}
};

template<class S>
Generic1DArrayWidget<S>::Generic1DArrayWidget(
  QString name, 
  QString type, 
  const RCP<const ParameterEntryValidator> validator,
  QWidget *parent):
  GenericArrayWidget<S>(name, type, validator, parent)
{}


template<class S>
QLayout* Generic1DArrayWidget<S>::getArrayLayout(){
  QGridLayout *widgetLayout = new QGridLayout;
  for(int i=0; i<baseArray.size(); ++i){
	  widgetLayout->addWidget(new QLabel("Item: " +QString::number(i)),0,i,Qt::AlignLeft);
	  QWidget* editorWidget = getEditorWidget(i);
	  widgetLayout->addWidget(editorWidget,1,i,Qt::AlignLeft);
	  widgetVector.push_back(editorWidget);
  }
  return widgetLayout;
}

template<class S>
void Generic1DArrayWidget<S>::doAcceptWork(){
  baseArray.clear();
  baseArray = getArrayFromWidgets();
  this->done(QDialog::Accepted);
}

/**
 * \brief A widget for editing Arrays of type int.
 */
class IntArrayWidget: public Generic1DArrayWidget<int>{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a IntArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	IntArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<int>(name, type, validator,parent){}

  //@}

public slots:

  /** \name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
  void accept(){
    doAcceptWork();
  }


private:

  /**
   * \brief .
   */
	QWidget* getEditorWidget(int index){
		QSpinBox *newSpin = new QSpinBox(this);
		RCP<const EnhancedNumberValidator<int> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<int>, int> >(getEntryValidator(),true)->getPrototype();
		}
		SpinBoxApplier<int>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}

  /**
   * \brief .
   */
	Array<int> getArrayFromWidgets(){
    Array<int> toReturn(widgetVector.size(), 0);
    for(size_t i=0; i < widgetVector.size(); ++i){
      toReturn[i]= ((QSpinBox*)widgetVector[i])->value();
    }
    return toReturn;
	}

  //@}
};

/**
 * \brief A widget for editing Arrays of type short.
 */
class ShortArrayWidget: public Generic1DArrayWidget<short>
{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a ShortArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	ShortArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<short>(name, type, validator,parent){}

  //@}

public slots:
  /** \name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
	void accept(){
    doAcceptWork();
	}

  //@}

private:
  /** * \brief .  */
	QWidget* getEditorWidget(int index){
		QSpinBox *newSpin = new QSpinBox(this);
		RCP<const EnhancedNumberValidator<short> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<short>, short> >(getEntryValidator(),true)->getPrototype();
		}
		SpinBoxApplier<short>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}

  /** * \brief .  */
	Array<short> getArrayFromWidgets(){
    Array<short> toReturn(widgetVector.size(), 0);
    for(size_t i=0; i < widgetVector.size(); ++i){
      toReturn[i]= ((QSpinBox*)widgetVector[i])->value();
    }
    return toReturn;
	}

  //@}
};

/**
 * \brief A widget for editing Arrays of type double.
 */
class DoubleArrayWidget: public Generic1DArrayWidget<double>
{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a DoubleArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	DoubleArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<double>(name, type, validator,parent){}

  //@}
  
public slots:
  /** \name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
	void accept(){
    doAcceptWork();
	}

private:

  /** * \brief .  */
	QWidget* getEditorWidget(int index){
		QDoubleSpinBox *newSpin = new QDoubleSpinBox(this);
		RCP<const EnhancedNumberValidator<double> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<double>, double> >(getEntryValidator(),true)->getPrototype();
		}
		SpinBoxApplier<double>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}

  /** * \brief .  */
	Array<double> getArrayFromWidgets(){
    Array<double> toReturn(widgetVector.size(), 0.0);
    for(size_t i=0; i < widgetVector.size(); ++i){
      toReturn[i]= ((QSpinBox*)widgetVector[i])->value();
    }
    return toReturn;
  }

  //@}
};

/**
 * \brief A widget for editing Arrays of type short.
 */
class FloatArrayWidget: public Generic1DArrayWidget<float>
{
	Q_OBJECT
public:
  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a FloatArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
	FloatArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<float>(name, type, validator,parent){}

  //@}
  
public slots:

  /** \name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
	void accept(){
    doAcceptWork();
	}

private:

  /** * \brief .  */
	QWidget* getEditorWidget(int index){
		QDoubleSpinBox *newSpin = new QDoubleSpinBox(this);
		RCP<const EnhancedNumberValidator<float> > validator = null;
		if(!is_null(getEntryValidator())){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<float>, float> >(getEntryValidator(),true)->getPrototype();
		}
		SpinBoxApplier<float>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}

  /** * \brief .  */
	Array<float> getArrayFromWidgets(){
    Array<float> toReturn(widgetVector.size(), 0.0);
    for(size_t i=0; i < widgetVector.size(); ++i){
      toReturn[i]= ((QSpinBox*)widgetVector[i])->value();
    }
    return toReturn;
  }

  //@}
};

/**
 * \brief A widget for editing an array of strings
 */
class StringArrayWidget : public Generic1DArrayWidget<std::string>{
	Q_OBJECT
public:

  /** \name Constructors */
  //@{

  /**
   * \brief Constructs a StringArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
  StringArrayWidget(
    QString name, 
    QString type, 
    const RCP<const ParameterEntryValidator> validator,
    QWidget *parent=0):
    Generic1DArrayWidget<std::string>(name, type, validator,parent)
  {}

  //@}

  /** \name Overridden from Generic1DArrayWidget */
  //@{

  /** * \brief .  */
	void accept(){
    doAcceptWork();
	}


private:

  /** * \brief .  */
	QWidget* getEditorWidget(int index){
    QString currentData = QString::fromStdString(baseArray[index]);
		if(is_null(getEntryValidator())){
			return new QLineEdit(currentData,this);
		}
		else if(!is_null(rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(getEntryValidator()))){
			return new FileNameWidget(
        currentData, 
        rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(getEntryValidator())->getPrototype()->fileMustExist(), this);
		}
		else if(getEntryValidator()->validStringValues()->size() != 0){
			RCP<const Array<std::string> > options = getEntryValidator()->validStringValues();
			QComboBox *newCombo = new QComboBox(this);
			for(Array<std::string>::const_iterator itr = options->begin(); itr != options->end(); ++itr){
				newCombo->addItem(QString::fromStdString(*itr));
			}
      int selectedItem = newCombo->findText(currentData);
      newCombo->setCurrentIndex(selectedItem);
			return newCombo;
		}
		else{
			return new QLineEdit(currentData,this);
		}
	}

  /** * \brief .  */
	Array<std::string> getArrayFromWidgets(){
    Array<std::string> toReturn(widgetVector.size(), "");
    for(size_t i=0; i < widgetVector.size(); ++i){
			if(is_null(getEntryValidator())){
        toReturn[i] = ((QLineEdit*)widgetVector[i])->text().toStdString();
			}
			else if(!is_null(rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(getEntryValidator()))){
        toReturn[i] = ((FileNameWidget*)widgetVector[i])->getCurrentFileName().toStdString();
			}
			else if(getEntryValidator()->validStringValues()->size() !=0){
        toReturn[i] = ((QComboBox*)widgetVector[i])->currentText().toStdString();
			}
			else{
        toReturn[i] = ((QLineEdit*)widgetVector[i])->text().toStdString();
			}
		}
    return toReturn;
  }

  //@}
};


} //end namespace
#endif //OPTIKA_ARRAYWIDGET_HPP_
