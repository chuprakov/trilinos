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
#ifndef OPTIKA_SPINBOXAPPLIER_HPP_
#define OPTIKA_SPINBOXAPPLIER_HPP_
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <float.h>
#include <limits>
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace Optika{

template<class S> class SpinBoxApplier{
public:

	static void applyToSpinBox(Teuchos::RCP<const Teuchos::EnhancedNumberValidator<S> > validator, QSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(Teuchos::EnhancedNumberTraits<S>::min());
			spinBox->setMaximum(Teuchos::EnhancedNumberTraits<S>::max());
			spinBox->setSingleStep(Teuchos::EnhancedNumberTraits<S>::defaultStep());
		}
	}
};

template <>
class SpinBoxApplier<int>{
public:

	/**
	 * Applies an EnhancedNumberValidator of type int to a QSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 */
	static void applyToSpinBox(Teuchos::RCP<const Teuchos::EnhancedNumberValidator<int> > validator, QSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(Teuchos::EnhancedNumberTraits<int>::min());
			spinBox->setMaximum(Teuchos::EnhancedNumberTraits<int>::max());
			spinBox->setSingleStep(Teuchos::EnhancedNumberTraits<int>::defaultStep());
		}
	}
};

/**
 * A specific validator used to validate values of type short.
 */
template<>
class SpinBoxApplier<short>{
public:

	/**
	 * Applies an EnhancedNumberValidator of type short to a QSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 */
	static void applyToSpinBox(Teuchos::RCP<const Teuchos::EnhancedNumberValidator<short> > validator, QSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(Teuchos::EnhancedNumberTraits<short>::min());
			spinBox->setMaximum(Teuchos::EnhancedNumberTraits<short>::max());
			spinBox->setSingleStep(Teuchos::EnhancedNumberTraits<short>::defaultStep());
		}
	}
};

/**
 * A specific validator used to validate entry's of type long long int.
 */
/*template<long long int>
class EnhancedNumberValidator<long long int> : public Teuchos::EnhancedNumberValidator<long long int>{
public:
	**
	 * Construcsts an EnhancedNumberValidator of type long long with no
	 * minimum or maximum.
	 *
	EnhancedNumberValidator():Teuchos::EnhancedNumberValidator(longlongId, longlongDefaultStep){}
	**
	 * Constructs an Enhanced number validator for type long long int.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 *
	EnhancedNumberValidator(long long int min, long long int max, long long int step=longlongDefaultStep):
	:Teuchos::EnhancedNumberValidator<long long int>(longlongId, min, max, step){}

	**
	 * Applies an EnhancedNumberValidator of type long long int to a QwwLongSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 *
	static void applyToSpinBox(Teuchos::RCP<const EnhancedNumberValidator<long long int> > validator, QwwLongSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(tempMin);
			spinBox->setMaximum(validator->max());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(LLONG_MIN);
			spinBox->setMaximum(LLONG_MAX);
			spinBox->setSingleStep(longlongDefaultStep);
		}
	}
};*/

template<>
class SpinBoxApplier<double>{
public:

	/**
	 * Applies an EnhancedNumberValidator of type double to a QDoubleSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The SpinBox on which to apply the validator.
	 */
	static void applyToSpinBox(Teuchos::RCP<const Teuchos::EnhancedNumberValidator<double> > validator, QDoubleSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
			spinBox->setDecimals(validator->getPrecision());
		}
		else{
			spinBox->setMinimum(Teuchos::EnhancedNumberTraits<double>::min());
			spinBox->setMaximum(Teuchos::EnhancedNumberTraits<double>::max());
			spinBox->setSingleStep(Teuchos::EnhancedNumberTraits<double>::defaultStep());
			spinBox->setDecimals(Teuchos::EnhancedNumberTraits<float>::defaultPrecision());
		}
	}
};

template<>
class SpinBoxApplier<float>{
public:

	/**
	 * Applies an EnhancedNumberValidator of type float to a QDoubleSpinBox.
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The SpinBox on which to apply the validator.
	 */
	static void applyToSpinBox(Teuchos::RCP<const Teuchos::EnhancedNumberValidator<float> > validator, QDoubleSpinBox *spinBox){
		if(!Teuchos::is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
			spinBox->setDecimals(validator->getPrecision());
		}
		else{
			spinBox->setMinimum(Teuchos::EnhancedNumberTraits<float>::min());
			spinBox->setMaximum(Teuchos::EnhancedNumberTraits<float>::max());
			spinBox->setSingleStep(Teuchos::EnhancedNumberTraits<float>::defaultStep());
			spinBox->setDecimals(Teuchos::EnhancedNumberTraits<float>::defaultPrecision());
		}
	}
}; 

}
#endif // OPTIKA_SPINBOXAPPLIER_HPP_
