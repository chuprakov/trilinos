#include "Optika_ValidatorFactory.hpp"
#include "Teuchos_ENull.hpp"
namespace Optika{

Teuchos::RCP<Teuchos::ParameterEntryValidator> ValidatorFactory::createValidator(ValidatorType valiType){
	switch(valiType){
		case Int:
			return Teuchos::RCP<EnhancedNumberValidator<int> >(new EnhancedNumberValidator<int>());
			break;
		case Short:
			return Teuchos::RCP<EnhancedNumberValidator<short> >(new EnhancedNumberValidator<short>());
			break;
		case Double:
			return Teuchos::RCP<EnhancedNumberValidator<double> >(new EnhancedNumberValidator<double>());
			break;
		case Float:
			return Teuchos::RCP<EnhancedNumberValidator<float> >(new EnhancedNumberValidator<float>());
			break;
		case IntArray:
			return Teuchos::RCP<ArrayNumberValidator<int> >( new ArrayNumberValidator<int>( Teuchos::RCP<EnhancedNumberValidator<int> >( new EnhancedNumberValidator<int>())));
			break;
		case ShortArray:
			return Teuchos::RCP<ArrayNumberValidator<short> >( new ArrayNumberValidator<short>( Teuchos::RCP<EnhancedNumberValidator<short> >( new EnhancedNumberValidator<short>())));
			break;
		case DoubleArray:
			return Teuchos::RCP<ArrayNumberValidator<double> >( new ArrayNumberValidator<double>( Teuchos::RCP<EnhancedNumberValidator<double> >( new EnhancedNumberValidator<double>())));
			break;
		case FloatArray:
			return Teuchos::RCP<ArrayNumberValidator<float> >( new ArrayNumberValidator<float>(Teuchos::RCP<EnhancedNumberValidator<float> >( new EnhancedNumberValidator<float>())));
			break;
		case FileName:
			return Teuchos::RCP<FileNameValidator>(new FileNameValidator());
			break;
		case FileNameArray:
			return Teuchos::RCP<ArrayFileNameValidator>(new ArrayFileNameValidator(Teuchos::RCP<FileNameValidator>(new FileNameValidator())));
			break;
		default:
			Teuchos::RCP<Teuchos::ParameterEntryValidator> toReturn;
			return toReturn;
			break;
	}
}

Teuchos::RCP<EnhancedNumberValidator<int> > ValidatorFactory::getIntValidator(){
	return Teuchos::RCP<EnhancedNumberValidator<int> >(new EnhancedNumberValidator<int>());
}

Teuchos::RCP<EnhancedNumberValidator<short> > ValidatorFactory::getShortValidator(){
	return Teuchos::RCP<EnhancedNumberValidator<short> >(new EnhancedNumberValidator<short>());
}

Teuchos::RCP<EnhancedNumberValidator<double> > ValidatorFactory::getDoubleValidator(){
	return Teuchos::RCP<EnhancedNumberValidator<double> >(new EnhancedNumberValidator<double>());
}

Teuchos::RCP<EnhancedNumberValidator<float> > ValidatorFactory::getFloatValidator(){
	return Teuchos::RCP<EnhancedNumberValidator<float> >(new EnhancedNumberValidator<float>());
}

Teuchos::RCP<FileNameValidator> ValidatorFactory::getFileNameValidator(){
	return Teuchos::RCP<FileNameValidator>(new FileNameValidator());
}

Teuchos::RCP<ArrayNumberValidator<int> > ValidatorFactory::getArrayIntValidator(){
	return Teuchos::RCP<ArrayNumberValidator<int> >( new ArrayNumberValidator<int>( Teuchos::RCP<EnhancedNumberValidator<int> >( new EnhancedNumberValidator<int>())));
}

Teuchos::RCP<ArrayNumberValidator<short> > ValidatorFactory::getArrayShortValidator(){
	return Teuchos::RCP<ArrayNumberValidator<short> >( new ArrayNumberValidator<short>( Teuchos::RCP<EnhancedNumberValidator<short> >( new EnhancedNumberValidator<short>())));
}

Teuchos::RCP<ArrayNumberValidator<double> > ValidatorFactory::getArrayDoubleValidator(){
	return Teuchos::RCP<ArrayNumberValidator<double> >( new ArrayNumberValidator<double>( Teuchos::RCP<EnhancedNumberValidator<double> >( new EnhancedNumberValidator<double>())));
}

Teuchos::RCP<ArrayNumberValidator<float> > ValidatorFactory::getArrayFloatValidator(){
	return Teuchos::RCP<ArrayNumberValidator<float> >( new ArrayNumberValidator<float>( Teuchos::RCP<EnhancedNumberValidator<float> >( new EnhancedNumberValidator<float>())));
}

Teuchos::RCP<ArrayFileNameValidator> ValidatorFactory::getArrayFileNameValidator(){
	return Teuchos::RCP<ArrayFileNameValidator>(new ArrayFileNameValidator(Teuchos::RCP<FileNameValidator>(new FileNameValidator())));
}

}
