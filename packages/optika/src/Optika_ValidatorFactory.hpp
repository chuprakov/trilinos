#ifndef OPTIKA_VALIDATORFACTORY_HPP_
#define OPTIKA_VALIDATORFACTORY_HPP_
#include "Optika_SpecificParameterEntryValidators.hpp"


namespace Optika{

class ValidatorFactory{
public:
	enum ValidatorType{IntValidator, ShortValidator, DoubleValidator, FloatValidator, FilenameValidator,
	IntArrayValidator, ShortArrayValidator, DoubleArrayValidator, FloatArrayValidator, ArrayFilenameValidator};

	static Teuchos::RCP<Teuchos::ParameterEntryValidator> createValidator(ValidatorType valiType);
};

}

#endif /* OPTIKA_VALIDATORFACTORY_HPP_ */
