// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLE.h

#ifndef MERIT_FUNC_NLE_H
#define MERIT_FUNC_NLE_H

#include <iosfwd>

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** Base class for all merit functions for systems of NonLinear Equations (NLE) {abstract}.
  */
class MeritFuncNLE {
public:

	///
	class InvalidInitialization : public std::logic_error
	{public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	virtual ~MeritFuncNLE() {}

	///
	/** Return the value of the merit function at c(x).
	  * This interface requires the client to compute
	  * c(x) and pass it to this function to have
	  * the value of phi(c) calculated.
	  *
	  * If the merit function has not been initialized properly
	  * then a #InvalidInitialization# exception will be thrown.
	  */
	virtual value_type value(const VectorSlice& c) const = 0;

	///
	/** Return the value of the directional derivative of the 
	  * merit function w.r.t. alpha at alpha = 0.  In other words
	  * compute return d( phi(c(x)) ) / d(alpha_k) at alpha_k = 0
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

};	// end class MeritFuncNLE

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_NLE_H