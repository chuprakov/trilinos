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
