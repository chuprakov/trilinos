// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncCalc1D.h

#ifndef MERIT_FUNC_CALC_1D_H
#define MERIT_FUNC_CALC_1D_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** Abstracts a 1D merit function {abstract}.
  *
  * This is the interface that line search algorithms use to compute
  * the value of the merit function at alpha (phi(alpha)) and
  * to retrieve the initial descent derivative of the merit function
  * (deriv()).
  */
class MeritFuncCalc1D {
public:

	///
	virtual ~MeritFuncCalc1D() {}

	/// Return the value of the merit function at alpha.
	virtual value_type operator()(value_type alpha) const = 0;

	/// Return the derivative of the merit function at alpha = 0
	virtual value_type deriv() const = 0;

	/// Print the particular merit function
	virtual void print_merit_func(std::ostream& out
		, const std::string& leading_str) const = 0;

};	// end class MeritFuncCalc1D

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_CALC_1D_H