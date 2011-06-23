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
#include <QtGui>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QComboBox>
#include <QFileDialog>
#include "Optika_delegate.hpp"
#include "float.h"
#include <iostream>
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace Optika{

Delegate::Delegate(QObject *parent):QItemDelegate(parent){}

QWidget* Delegate::createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/ , const QModelIndex &index ) const{
	QWidget *editor = 0;
	if(index.column() != 1)
		return editor;

	RCP<const ParameterEntryValidator> paramValidator = ((TreeModel*)(index.model()))->getValidator(index);
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);

	if(itemType == intId){
		editor = new QSpinBox(parent);
		RCP<const EnhancedNumberValidator <int> > intValidator;
		if(!is_null(paramValidator)){
			intValidator = rcp_dynamic_cast<const EnhancedNumberValidator<int> >(paramValidator);
		}
		SpinBoxApplier<int>::applyToSpinBox(intValidator, (QSpinBox*)editor);
	}
	else if(itemType == shortId){
		editor = new QSpinBox(parent);
		RCP<const EnhancedNumberValidator<short> > shortValidator;
		if(!is_null(paramValidator)){
			shortValidator = rcp_dynamic_cast<const EnhancedNumberValidator<short> >(paramValidator);
		}
		SpinBoxApplier<short>::applyToSpinBox(shortValidator, (QSpinBox*)editor);
	}
/*	else if(itemType == longlongId){
		editor = new QwwLongSpinBox(parent);
		RCP<const EnhancedNumberValidator<long long> > longlongValidator;
		if(!is_null(paramValidator)){
			longlongValidator = rcp_dynamic_cast<const EnhancedNumberValidator<long long> >(paramValidator);
		}
		EnhancedNumberValidator<long long>::applyToSpinBox(longlongValidator, (QDoubleSpinBox*)editor);
	}*/
	else if(itemType == doubleId){
		editor = new QDoubleSpinBox(parent);
		RCP<const EnhancedNumberValidator<double> > doubleValidator;
		if(!is_null(paramValidator)){
			doubleValidator = rcp_dynamic_cast<const EnhancedNumberValidator<double> >(paramValidator);
		}
		SpinBoxApplier<double>::applyToSpinBox(doubleValidator, (QDoubleSpinBox*)editor);
	}
	else if(itemType == floatId){
		editor = new QDoubleSpinBox(parent);
		RCP<const EnhancedNumberValidator<float> > floatValidator; 
		if(!is_null(paramValidator)){
			floatValidator = rcp_dynamic_cast<const EnhancedNumberValidator<float> >(paramValidator);
		}
		SpinBoxApplier<float>::applyToSpinBox(floatValidator, (QDoubleSpinBox*)editor);
	}
	else if(itemType == boolId){
		editor = new QComboBox(parent);
		static_cast<QComboBox*>(editor)->addItem(getBoolEditorTrue());
		static_cast<QComboBox*>(editor)->addItem(getBoolEditorFalse());
	}
	else if(itemType == stringId){
		if(is_null(paramValidator)){
			editor = new QLineEdit(parent);
		}
		else if(!is_null(rcp_dynamic_cast<const FileNameValidator>(paramValidator))){
			QString paramName = 
				((TreeModel*)(index.model()))->data(index.sibling(index.row(), 0),Qt::DisplayRole).toString();
			QString currentPath = ((TreeModel*)(index.model()))->data(index,Qt::DisplayRole).toString();
			if(currentPath.size() == 0){
				currentPath = QDir::homePath();
			}
			QString filename;
			if(rcp_dynamic_cast<const FileNameValidator>(paramValidator)->fileMustExist()){
				filename = QFileDialog::getOpenFileName(parent, paramName, currentPath);
			}
			else{
				filename = QFileDialog::getSaveFileName(parent, paramName, currentPath);
			}
			if(filename != ""){
				((TreeModel*)(index.model()))->setData(index, filename);
			}
		}
		else if(paramValidator->validStringValues()->size() != 0){
			RCP<const Array<std::string> > options = paramValidator->validStringValues();
			editor = new QComboBox(parent);
			for(Array<std::string>::const_iterator itr = options->begin(); itr != options->end(); ++itr){
				static_cast<QComboBox*>(editor)->addItem(QString::fromStdString(*itr));
			}
		}
		else{
			editor = new QLineEdit(parent);
		}
	}
	else if(itemType.contains(arrayId)){
		editor = getArrayEditor(index, getArrayType(itemType), parent);
	}

	return editor;
}

void Delegate::setEditorData(QWidget *editor, const QModelIndex &index) const{
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);
	if(itemType == intId){
		int value = index.model()->data(index).toInt();
		static_cast<QSpinBox*>(editor)->setValue(value);
	}
	else if(itemType == shortId){
		short value = index.model()->data(index).toInt();
		static_cast<QSpinBox*>(editor)->setValue(value);
	}
	else if(itemType == doubleId){
		double value = index.model()->data(index).toDouble();
		static_cast<QDoubleSpinBox*>(editor)->setValue(value);
	}
	else if(itemType == floatId){
		float value = index.model()->data(index).toDouble();
		static_cast<QDoubleSpinBox*>(editor)->setValue(value);
	}
	else if(itemType == boolId){
		QString value = index.model()->data(index).toString();
		static_cast<QComboBox*>(editor)->setEditText(value);
	}
	else if(itemType == stringId){
		QString value = index.model()->data(index).toString();
		RCP<const ParameterEntryValidator> validator = ((TreeModel*)(index.model()))->getValidator(index);
		if(is_null(validator) || validator->validStringValues()->size()==0)
			static_cast<QLineEdit*>(editor)->setText(value);
	 	else
			static_cast<QComboBox*>(editor)->setEditText(value);
	}
  else if(itemType.contains(arrayId)){
    setArrayWidgetData(editor, itemType.section(" ", -1), index);
  }
}

void Delegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const{
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);
	if(itemType == intId){
		QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, spinBox->value(), Qt::EditRole);
	}
	else if(itemType == shortId){
		QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, (short)spinBox->value(), Qt::EditRole);
	}
	else if(itemType == doubleId){
		QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, spinBox->value(), Qt::EditRole);
	}
	else if(itemType == floatId){
		QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, (float)spinBox->value(), Qt::EditRole);
	}
	else if(itemType == boolId){
		bool value = static_cast<QComboBox*>(editor)->currentText() 
      == getBoolEditorTrue(); 
		model->setData(index, value, Qt::EditRole);
	}
	else if(itemType == stringId){
		RCP<const ParameterEntryValidator> validator = 
			static_cast<const TreeModel*>(index.model())->getValidator(index);
		QString value;
		if(is_null(validator)){
			value = static_cast<QLineEdit*>(editor)->text();
		}
		else{
			value = static_cast<QComboBox*>(editor)->currentText(); 
		}
		model->setData(index, value, Qt::EditRole);
	}
  else if(itemType.contains(arrayId)){
    QVariant value = extractValueFromArray(editor, itemType.section(" ", -1));
    model->setData(index, value, Qt::EditRole);
  }
}

QVariant Delegate::extractValueFromArray(QWidget* editor, QString itemType) const
{
  if(itemType == intId){
    IntArrayWidget* intEditor = (IntArrayWidget*)editor;
    return QVariant::fromValue(intEditor->getData());
  }
  else if(itemType == shortId){
    ShortArrayWidget* shortEditor = (ShortArrayWidget*)editor;
    return QVariant::fromValue(shortEditor->getData());
  }
  else if(itemType == floatId){
    FloatArrayWidget* floatEditor = (FloatArrayWidget*)editor;
    return QVariant::fromValue(floatEditor->getData());
  }
  else if(itemType == doubleId){
    DoubleArrayWidget* doubleEditor = (DoubleArrayWidget*)editor;
    return QVariant::fromValue(doubleEditor->getData());
  }
  else if(itemType == stringId){
    StringArrayWidget* stringEditor = (StringArrayWidget*)editor;
    return QVariant::fromValue(stringEditor->getData());
  }
  else{
    return QVariant();
  }
}
 

void Delegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &/*index*/) const{
	editor->setGeometry(option.rect);
}

QWidget* Delegate::getArrayEditor(const QModelIndex& index, QString type, QWidget *parent) const{
  TreeModel* model = (TreeModel*)index.model();
  QString name = model->data(
    index.sibling(index.row(),0),Qt::DisplayRole).toString();
  RCP<const ParameterEntryValidator> validator = 
    model->getValidator(index);
    
	if(type == intId){
		return new IntArrayWidget(name, type, validator, parent);
	}
	else if(type == shortId){
		return new ShortArrayWidget(name, type, validator, parent);
	}
	/*else if(type == longlongId){
		return new IntArrayWidget(index, type, parent);
	}*/
	else if(type == doubleId){
		return new DoubleArrayWidget(name, type, validator, parent);
  }
	else if(type == floatId){
		return new FloatArrayWidget(name, type, validator, parent);
	}
	else if(type == stringId){
		return new StringArrayWidget(name, type, validator, parent);
	}
  else{
    return 0;
  }
}

void Delegate::setArrayWidgetData(QWidget* editor, QString type, const QModelIndex& index) const{
  QVariant newData = index.model()->data(index, TreeModel::getRawDataRole());
	if(type == intId){
    ((IntArrayWidget*)editor)->initData(newData.value<Array<int> >());
	}
	else if(type == shortId){
    ((ShortArrayWidget*)editor)->initData(newData.value<Array<short> >());
	}
	/*else if(type == longlongId){
	}*/
	else if(type == doubleId){
    ((DoubleArrayWidget*)editor)->initData(newData.value<Array<double> >());
  }
	else if(type == floatId){
    ((FloatArrayWidget*)editor)->initData(newData.value<Array<float> >());
	}
	else if(type == stringId){
    ((StringArrayWidget*)editor)->initData(
      newData.value<Array<std::string> >());
	}
}


} //End namespace

