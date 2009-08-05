#include "Optika_SpecificParameterEntryValidators.hpp"

namespace Optika{

FileNameValidator::FileNameValidator():ParameterEntryValidator(){}

Teuchos::RCP<const Teuchos::Array<std::string> > FileNameValidator::validStringValues() const{
	return Teuchos::null;
}

void FileNameValidator::validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
	Teuchos::any anyValue = entry.getAny(true);
	if(!(anyValue.type() == typeid(std::string) )){
		const std::string &entryName = entry.getAny(false).typeName();
		std::stringstream oss;
		std::string msg;
		oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
		" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
		"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
		"can help you figure out what went wrong.\n\n"
		"Error: The value that you entered was the wrong type." <<
		"Parameter: " << paramName << "\n" << 
		"Type specified: " << entryName << "\n" <<
		"Type accepted: " << typeid(std::string).name() << "\n";
		msg = oss.str();
		throw Teuchos::Exceptions::InvalidParameterType(msg);
	}
}

void FileNameValidator::printDoc(std::string const &docString, std::ostream &out) const{
	Teuchos::StrUtils::printLines(out,"# ",docString);
	out << "#  Validator Used: \n";
	out << "#	FileName Validator\n";
}

StringValidator::StringValidator(Teuchos::Array<std::string> validStrings):ParameterEntryValidator(){
	this->validStrings = validStrings;	
}

const Teuchos::Array<std::string> StringValidator::setValidStrings(Teuchos::Array<std::string> validStrings){
	this->validStrings = validStrings;
	return this->validStrings;
}

Teuchos::RCP<const Teuchos::Array<std::string> > StringValidator::validStringValues() const{
	Teuchos::RCP<const Teuchos::Array<std::string> > toReturn = Teuchos::RCP<const Teuchos::Array<std::string> >(new Teuchos::Array<std::string>(validStrings));
	return toReturn;
}

void StringValidator::validate(Teuchos::ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
	Teuchos::any anyValue = entry.getAny(true);
	if(!(anyValue.type() == typeid(std::string) )){
		const std::string &entryName = entry.getAny(false).typeName();
		std::stringstream oss;
		std::string msg;
		oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
		" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
		"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
		"can help you figure out what went wrong.\n\n"
		"Error: The value that you entered was the wrong type." <<
		"Parameter: " << paramName << "\n" << 
		"Type specified: " << entryName << "\n" <<
		"Type accepted: " << typeid(std::string).name() << "\n";
		msg = oss.str();
		throw Teuchos::Exceptions::InvalidParameterType(msg);
	}
	else{
		Teuchos::Array<std::string>::const_iterator it = std::find(validStrings.begin(), validStrings.end(), Teuchos::getValue<std::string>(entry));
		if(it == validStrings.end()){
			std::stringstream oss;
			std::string msg;
			oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
			" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
			"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
			"can help you figure out what went wrong.\n\n"
			"Error: The value that was entered doesn't fall with in " <<
			"the range set by the validator." <<
			"Parameter: " << paramName << "\n" <<
			"Acceptable Values: " << validStrings << "\n" <<
			"Value entered: " << Teuchos::getValue<std::string>(entry) << "\n";
			msg = oss.str();
			throw Teuchos::Exceptions::InvalidParameterValue(msg);
		}

	}
		

}

void StringValidator::printDoc(std::string const &docString, std::ostream &out) const{
	Teuchos::StrUtils::printLines(out,"# ",docString);
	out << "#  Validator Used: \n";
	out << "#	FileName Validator\n";
}



}

