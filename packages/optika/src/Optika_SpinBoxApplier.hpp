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
#include "Optika_ConfigDefs.hpp"

/*! \file Optika_SpinBoxApplier.hpp
    \brief A collection of objects which
    apply the restriction of a EnhancedNumberValidator
    to a spinbox.
*/

namespace Optika{

/**
 * \brief A templated class that applies the minimum,
 * maximum and step specified in an EnhancedNumberValidator
 * onto a QSpingBox.
 */
template<class S> class SpinBoxApplier{
public:
  //! @name Applier Functions
  //@{

  /**
   * \brief Applied attributes of the validator to the spin box.
   *
   * @param validator Validator whose attributes are to be applied to
   * the spinbox
   * @param spinBox A point to the spinbox upon which the validator's
   * attributes should be applied.
   */
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<S> > validator, QSpinBox *spinBox){
		if(!is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(EnhancedNumberTraits<S>::min());
			spinBox->setMaximum(EnhancedNumberTraits<S>::max());
			spinBox->setSingleStep(EnhancedNumberTraits<S>::defaultStep());
		}
	}

  //@}

};

/**
 * \brief Template specialzation on type in of the SpinBoxApplier
 */
template <>
class SpinBoxApplier<int>{
public:
  
  //! @name Applier Functions
  //@{

	/**
	 * \brief Applies an EnhancedNumberValidator of type int to a QSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 */
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<int> > validator, QSpinBox *spinBox){
		if(!is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(EnhancedNumberTraits<int>::min());
			spinBox->setMaximum(EnhancedNumberTraits<int>::max());
			spinBox->setSingleStep(EnhancedNumberTraits<int>::defaultStep());
		}
	}

  //@}
};

/**
 * \brief A specific validator used to validate values of type short.
 */
template<>
class SpinBoxApplier<short>{
public:

  //! @name Applier Functions
  //@{

	/**
	 * \brief Applies an EnhancedNumberValidator of type short to a QSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 */
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<short> > validator, QSpinBox *spinBox){
		if(!is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(EnhancedNumberTraits<short>::min());
			spinBox->setMaximum(EnhancedNumberTraits<short>::max());
			spinBox->setSingleStep(EnhancedNumberTraits<short>::defaultStep());
		}
	}

  //@}
};

/**
 * \brief A specific validator used to validate entry's of type long long int.
 */
/*template<long long int>
class EnhancedNumberValidator<long long int> : public EnhancedNumberValidator<long long int>{
public:
	**
	 * Construcsts an EnhancedNumberValidator of type long long with no
	 * minimum or maximum.
	 *
	EnhancedNumberValidator():EnhancedNumberValidator(longlongId, longlongDefaultStep){}
	**
	 * Constructs an Enhanced number validator for type long long int.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the value should be changed. This is mostly used for 
	 * the QSpinBox that is used in the Optika GUI. If you're not using the GUI, you may ignore this parameter.
	 *
	EnhancedNumberValidator(long long int min, long long int max, long long int step=longlongDefaultStep):
	:EnhancedNumberValidator<long long int>(longlongId, min, max, step){}

	**
	 * Applies an EnhancedNumberValidator of type long long int to a QwwLongSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The spinBox on which to apply the validator.
	 *
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<long long int> > validator, QwwLongSpinBox *spinBox){
		if(!is_null(validator)){
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

/**
 * \brief A template spcialization of the SpinBoxApplier
 * class on the type double.
 */
template<>
class SpinBoxApplier<double>{
public:

  //! @name Applier Functions
  //@{

	/**
	 * \brief Applies an EnhancedNumberValidator of type double to a QDoubleSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The SpinBox on which to apply the validator.
	 */
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<double> > validator, QDoubleSpinBox *spinBox){
		if(!is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
			spinBox->setDecimals(validator->getPrecision());
		}
		else{
			spinBox->setMinimum(EnhancedNumberTraits<double>::min());
			spinBox->setMaximum(EnhancedNumberTraits<double>::max());
			spinBox->setSingleStep(EnhancedNumberTraits<double>::defaultStep());
			spinBox->setDecimals(EnhancedNumberTraits<float>::defaultPrecision());
		}
	}

  //@}
};

/**
 * \brief A template specialzation of the SpinBoxApplier
 * class on the type float.
 */
template<>
class SpinBoxApplier<float>{
public:

  //! @name Applier Functions
  //@{

	/**
	 * \brief Applies an EnhancedNumberValidator of type float to a QDoubleSpinBox.
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The SpinBox on which to apply the validator.
	 */
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<float> > validator, QDoubleSpinBox *spinBox){
		if(!is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
			spinBox->setDecimals(validator->getPrecision());
		}
		else{
			spinBox->setMinimum(EnhancedNumberTraits<float>::min());
			spinBox->setMaximum(EnhancedNumberTraits<float>::max());
			spinBox->setSingleStep(EnhancedNumberTraits<float>::defaultStep());
			spinBox->setDecimals(EnhancedNumberTraits<float>::defaultPrecision());
		}
	}

  //@}
}; 

}
#endif // OPTIKA_SPINBOXAPPLIER_HPP_
