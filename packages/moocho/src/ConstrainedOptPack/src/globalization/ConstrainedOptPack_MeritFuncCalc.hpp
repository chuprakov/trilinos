// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncCalc.h

#ifndef MERIT_FUNC_CALC_H
#define MERIT_FUNC_CALC_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** Abstract iterface for n-D merit functions {abstract}.
  *
  * Used to compute the value of the merit at a point x (phi(x)) and
  * to retrieve the derivative (phi.deriv()) along some direction d from
  * some base point xo.
  */
class MeritFuncCalc  {
public:

	///
	virtual ~MeritFuncCalc() {}

	///
	/** Return the value of the merit function at x.
	  */
	virtual value_type operator()(const VectorSlice& x) const= 0;

	/// Calls value(d_k) on aggregate merit_func.
	virtual value_type deriv() const = 0;

	/// Print what this merit function is
	virtual void print_merit_func(std::ostream& out
		, const std::string& leading_str) const = 0;

};	// end class MeritFuncCalc

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_CALC_H
