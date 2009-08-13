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
#ifndef OPTIKA_VALIDATORFACTORY_HPP_
#define OPTIKA_VALIDATORFACTORY_HPP_
#include "Optika_SpecificParameterEntryValidators.hpp"


namespace Optika{

class ValidatorFactory{
public:
	enum ValidatorType{Int, Short, Double, Float, FileName,
	IntArray, ShortArray, DoubleArray, FloatArray, FileNameArray};

	static Teuchos::RCP<Teuchos::ParameterEntryValidator> createValidator(ValidatorType valiType);
	static Teuchos::RCP<EnhancedNumberValidator<int> > getIntValidator();
	static Teuchos::RCP<EnhancedNumberValidator<short> > getShortValidator();
	static Teuchos::RCP<EnhancedNumberValidator<double> > getDoubleValidator();
	static Teuchos::RCP<EnhancedNumberValidator<float> > getFloatValidator();
	static Teuchos::RCP<FileNameValidator> getFileNameValidator();
	static Teuchos::RCP<ArrayNumberValidator<int> > getArrayIntValidator();
	static Teuchos::RCP<ArrayNumberValidator<short> > getArrayShortValidator();
	static Teuchos::RCP<ArrayNumberValidator<double> > getArrayDoubleValidator();
	static Teuchos::RCP<ArrayNumberValidator<float> > getArrayFloatValidator();
	static Teuchos::RCP<ArrayFileNameValidator> getArrayFileNameValidator();
};

}

#endif /* OPTIKA_VALIDATORFACTORY_HPP_ */
