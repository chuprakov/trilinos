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

namespace Optika {


/**
 * A templated abstract base class for all other array editing widgets. Note the absence of the Q_OBJECT
 * macro. This is becuase classes using the Q_OBJECT macro can't be templated (bummer). The macro is therfore
 * present in the subclasses.
 */ 
template <class S>
class GenericArrayWidget : public QDialog{
public:

  /**
   * Constructs a GenericArrayWidget.
   *
   * @param name The name of the parmaeter being edited.
   * @param type The array's template type.
   * @param validator  The validator on the array (null if there is none).
   * @param parent The parent widget.
   */
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
	QString getType(){
		return type;
	}

	/**
   * Gathers all the user inputed data and closes the dialog.
	 */
	void doAcceptWork();

  /**
   * Return the array backing this widget.
   *
   * @param The array backing this widget.
   */
  const Array<S> getData() const{
    return baseArray;
  }

  /**
   * Initialize the data in the widget with the currecnt values of the baseArray as
   * well as setting up the widgets layout.
   */
  void initData(Array<S> array){
    baseArray = array;
    setupArrayLayout();
  }

  /**
   * Must be implemented as a slot in base classes.
   */
  virtual void accept() =0;

protected:
	/**
	 * Convienece typedef. Represents an array of QWidgets.
	 */
	typedef std::vector<QWidget*> WVector;

	/**
	 * Conatins the editing widgets (e.g. QLineEdits and QSpinBoxes) comprising the array editor.
	 */
	WVector widgetVector;

	/**
	 * The validator being used on the array.
	 */
	RCP<const ParameterEntryValidator> entryValidator;	


	/**
	 * The array to be edited.
	 */
	Array<S> baseArray;

  /**
   * The name of the Parameter being edited.
   */
  QString name;

private:
	/**
	 * The widget containing all of the editing widgets (e.g.
	 * QLineEdits, and QSpinBoxes) that comprise the array editor.
	 */
	QWidget *arrayContainer;

	/**
	 * The type of array.
	 */
	QString type;
	
	/**
	 * Gets the widget to be used as an editor for each entry in the array.
	 */
	virtual QWidget* getEditorWidget(int index) = 0;

	/**
	 * Sets up the layout for the arrayContainer, including adding what ever editing
	 * widget should be used for the particual type of array.
	 */
	void setupArrayLayout();

  /**
   * Get a new array reflecting the current values entered in the widgets.
   *
   * @return A new array reflecting the currecnt values entered in the widgets.
   */
  virtual Array<S> getArrayFromWidgets() = 0;
};

template<class S>
GenericArrayWidget<S>::GenericArrayWidget(
  QString name, 
  QString type, 
  const RCP<const ParameterEntryValidator> validator,
  QWidget *parent):
  QDialog(parent),
  entryValidator(validator),
  name(name),
  type(type)
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


template<class S>
void GenericArrayWidget<S>::setupArrayLayout(){
  if(arrayContainer->layout() == NULL){
	  QGridLayout *widgetLayout = new QGridLayout;
	  for(int i=0; i<baseArray.size(); ++i){
		  widgetLayout->addWidget(new QLabel("Item: " +QString::number(i)),0,i,Qt::AlignLeft);
		  QWidget* editorWidget = getEditorWidget(i);
		  widgetLayout->addWidget(editorWidget,1,i,Qt::AlignLeft);
		  widgetVector.push_back(editorWidget);
	  }
	  arrayContainer->setLayout(widgetLayout);
  }
}

template<class S>
void GenericArrayWidget<S>::doAcceptWork(){
  baseArray.clear();
  baseArray = getArrayFromWidgets();
  done(QDialog::Accepted);
}

/**
 * A widget for editing Arrays of type int.
 */
class IntArrayWidget: public GenericArrayWidget<int>{
	Q_OBJECT
public:

  /**
   * Constructs a IntArrayWidget.
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
    GenericArrayWidget<int>(name, type, validator,parent){}

public slots:
  void accept(){
    doAcceptWork();
  }

private:
	QWidget* getEditorWidget(int index){
		QSpinBox *newSpin = new QSpinBox(this);
		RCP<const EnhancedNumberValidator<int> > validator = null;
		if(!is_null(entryValidator)){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<int>, int> >(entryValidator,true)->getPrototype();
		}
		SpinBoxApplier<int>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}

	Array<int> getArrayFromWidgets(){
    Array<int> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			baseArray.push_back(((QSpinBox*)(*it))->value());
		}
    return toReturn;
	}
};

/**
 * A widget for editing Arrays of type short.
 */
class ShortArrayWidget: public GenericArrayWidget<short>
{
	Q_OBJECT
public:

  /**
   * Constructs a ShortArrayWidget.
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
    GenericArrayWidget<short>(name, type, validator,parent){}

public slots:
	void accept(){
    doAcceptWork();
	}

private:
	QWidget* getEditorWidget(int index){
		QSpinBox *newSpin = new QSpinBox(this);
		RCP<const EnhancedNumberValidator<short> > validator = null;
		if(!is_null(entryValidator)){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<short>, short> >(entryValidator,true)->getPrototype();
		}
		SpinBoxApplier<short>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}

	Array<short> getArrayFromWidgets(){
    Array<short> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			toReturn.push_back(((QSpinBox*)(*it))->value());
		}
    return toReturn;
  }
};
/**
 * A widget for editing Arrays of type double.
 */
class DoubleArrayWidget: public GenericArrayWidget<double>
{
	Q_OBJECT
public:

  /**
   * Constructs a DoubleArrayWidget.
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
    GenericArrayWidget<double>(name, type, validator,parent){}
public slots:
	void accept(){
    doAcceptWork();
	}

private:
	QWidget* getEditorWidget(int index){
		QDoubleSpinBox *newSpin = new QDoubleSpinBox(this);
		RCP<const EnhancedNumberValidator<double> > validator = null;
		if(!is_null(entryValidator)){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<double>, double> >(entryValidator,true)->getPrototype();
		}
		SpinBoxApplier<double>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}
	Array<double> getArrayFromWidgets(){
    Array<double> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			baseArray.push_back(((QDoubleSpinBox*)(*it))->value());
		}
    return toReturn;
  }
};

/**
 * A widget for editing Arrays of type short.
 */
class FloatArrayWidget: public GenericArrayWidget<float>
{
	Q_OBJECT
public:
  /**
   * Constructs a FloatArrayWidget.
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
    GenericArrayWidget<float>(name, type, validator,parent){}
public slots:
	void accept(){
    doAcceptWork();
	}

private:
	QWidget* getEditorWidget(int index){
		QDoubleSpinBox *newSpin = new QDoubleSpinBox(this);
		RCP<const EnhancedNumberValidator<float> > validator = null;
		if(!is_null(entryValidator)){
			validator = rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<float>, float> >(entryValidator,true)->getPrototype();
		}
		SpinBoxApplier<float>::applyToSpinBox(validator, newSpin);
    newSpin->setValue(baseArray[index]);
		return newSpin;
	}
	Array<float> getArrayFromWidgets(){
    Array<float> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			toReturn.push_back(((QDoubleSpinBox*)(*it))->value());
		}
    return toReturn;
  }
};

/**
 * A widget for editing an array of strings
 */
class StringArrayWidget : public GenericArrayWidget<std::string>{
	Q_OBJECT
public:
  /**
   * Constructs a StringArrayWidget.
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
    GenericArrayWidget<std::string>(name, type, validator,parent){}

	void accept(){
    doAcceptWork();
	}


private:
	QWidget* getEditorWidget(int index){
    QString currentData = QString::fromStdString(baseArray[index]);
		if(is_null(entryValidator)){
			return new QLineEdit(currentData,this);
		}
		else if(!is_null(rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(entryValidator))){
			return new FileNameWidget(
        currentData, 
        rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(entryValidator)->getPrototype()->fileMustExist(), this);
		}
		else if(entryValidator->validStringValues()->size() != 0){
			RCP<const Array<std::string> > options = entryValidator->validStringValues();
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

	Array<std::string> getArrayFromWidgets(){
    Array<std::string> toReturn;
		for(WVector::iterator it = widgetVector.begin(); it != widgetVector.end(); ++it){
			if(is_null(entryValidator)){
				toReturn.push_back(((QLineEdit*)(*it))->text().toStdString());
			}
			else if(!is_null(rcp_dynamic_cast<const ArrayValidator<FileNameValidator, std::string> >(entryValidator))){
				toReturn.push_back(((FileNameWidget*)(*it))->getCurrentFileName().toStdString());
			}
			else if(entryValidator->validStringValues()->size() !=0){
				toReturn.push_back(((QComboBox*)(*it))->currentText().toStdString());
			}
			else{
				toReturn.push_back(((QLineEdit*)(*it))->text().toStdString());
			}
		}
    return toReturn;
  }
};


} //end namespace
#endif //OPTIKA_ARRAYWIDGET_HPP_
