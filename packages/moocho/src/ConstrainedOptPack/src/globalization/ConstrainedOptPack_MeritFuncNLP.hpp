// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLP.h

#ifndef MERIT_FUNC_NLP_H
#define MERIT_FUNC_NLP_H

#include <iosfwd>

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** Base class for all merit functions for NonLinear Programs (NLP) {abstract}.
  */
class MeritFuncNLP {
public:

	///
	class InvalidInitialization : public std::logic_error
	{public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	virtual ~MeritFuncNLP() {}

	///
	/** Return the value of the merit function at f(x), c(x).
	  * This interface requires the client to compute f(x)
	  * and c(x) and pass it to this function to have
	  * the value of phi(f,c) calculated.
	  *
	  * If the merit function has not been initialized properly
	  * then a #InvalidInitialization# exception will be thrown.
	  */
	virtual value_type value(value_type f, const VectorSlice& c) const = 0;

	///
	/** Return the value of the directional derivative of the 
	  * merit function w.r.t. alpha at alpha = 0.  In other words
	  * compute return d( phi( f(x) , c(x) ) ) / d(alpha_k) at alpha_k = 0
	  * where x = x_k + alpha_k * d_k.
	  *
	  * If the merit function has not been initialized properly
	  * then a #InvalidInitialization# exception will be thrown.
	  */
	virtual value_type deriv() const = 0;

	///
	/** Print the merit funciton
	  */
	virtual void print_merit_func(std::ostream& out
		, const std::string& leading_str) const = 0;

};	// end class MeritFuncNLP

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_NLP_H
