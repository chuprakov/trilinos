#include "Optika_ValidatorFactory.hpp"
#include "Teuchos_ENull.hpp"

namespace Optika{

Teuchos::RCP<Teuchos::ParameterEntryValidator> ValidatorFactory::createValidator(ValidatorType valiType){
	switch(valiType){
		case IntValidator:
			return Teuchos::RCP<EnhancedNumberValidator<int> >(new EnhancedNumberValidator<int>());
			break;
		case ShortValidator:
			return Teuchos::RCP<EnhancedNumberValidator<short> >(new EnhancedNumberValidator<short>());
			break;
		case DoubleValidator:
			return Teuchos::RCP<EnhancedNumberValidator<double> >(new EnhancedNumberValidator<double>());
			break;
		case FloatValidator:
			return Teuchos::RCP<EnhancedNumberValidator<float> >(new EnhancedNumberValidator<float>());
			break;
		case IntArrayValidator:
			return Teuchos::RCP<ArrayNumberValidator<int> >( new ArrayNumberValidator<int>( Teuchos::RCP<EnhancedNumberValidator<int> >( new EnhancedNumberValidator<int>())));
			break;
		case ShortArrayValidator:
			return Teuchos::RCP<ArrayNumberValidator<short> >( new ArrayNumberValidator<short>( Teuchos::RCP<EnhancedNumberValidator<short> >( new EnhancedNumberValidator<short>())));
			break;
		case DoubleArrayValidator:
			return Teuchos::RCP<ArrayNumberValidator<double> >( new ArrayNumberValidator<double>( Teuchos::RCP<EnhancedNumberValidator<double> >( new EnhancedNumberValidator<double>())));
			break;
		case FloatArrayValidator:
			return Teuchos::RCP<ArrayNumberValidator<float> >( new ArrayNumberValidator<float>( Teuchos::RCP<EnhancedNumberValidator<float> >( new EnhancedNumberValidator<float>())));
			break;
		case FilenameValidator:
			return Teuchos::RCP<FileNameValidator>(new FileNameValidator());
			break;
		case ArrayFilenameValidator:
			return Teuchos::RCP<ArrayFileNameValidator>(new ArrayFileNameValidator(Teuchos::RCP<FileNameValidator>(new FileNameValidator())));
			break;
		default:
			Teuchos::RCP<Teuchos::ParameterEntryValidator> toReturn;
			return toReturn;
			break;
	}
}

}
