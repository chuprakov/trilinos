#ifndef OPTIKA_INVALIDDEPENDENCYEXCEPTION_HPP_
#define OPTIKA_INVALIDDEPENDENCYEXCEPTION_HPP_
#include <stdexcept>

namespace Optika {


/**
 * Thrown when some aspect of a Dependency has been determined to be invalid.
 */
class InvalidDependencyException : public std::logic_error{
public: 
	/**
	 * Constructs an InvalidDependencyException
	 *
	 * @param what_arg The error message to be associated with this error.
	 */
	InvalidDependencyException(const std::string& what_arg):std::logic_error(what_arg){}
};



}
#endif //OPTIKA_INVALIDDEPENDENCYEXCEPTION_HPP_
