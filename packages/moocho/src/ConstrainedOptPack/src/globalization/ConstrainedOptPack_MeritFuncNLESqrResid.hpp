// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLESqrResid.h

#ifndef MERIT_FUNC_NLE_SQR_RESID_H
#define MERIT_FUNC_NLE_SQR_RESID_H

#include "MeritFuncNLE.h"

namespace ConstrainedOptimizationPack {

///
/** A merit function for the square of the constriant values.
  *
  * phi(x) = 1/2 * c(x)'*c(x)
  *
  * Dphi(x_k,d_k) = - c(x)'*c(x)
  *
  * Note that the definition of Dphi(x_k,d_k) assumes
  * that Gc_k'*d_k + c_k = 0.  In otherwords, d_k must
  * satisfiy the linearized equality constraints at
  * at x_k.
  *
  * Implicit copy constructor and assignment operators
  * are allowed.
  */
class MeritFuncNLESqrResid : public MeritFuncNLE {
public:

	/// Initializes deriv() = 0
	MeritFuncNLESqrResid();

	///
	value_type calc_deriv( const VectorSlice& c_k );

	// ////////////////////////////////
	// Overridden from MeritFuncNLE

	///
	value_type value(const VectorSlice& c) const;

	///
	value_type deriv() const;

	///
	void print_merit_func(std::ostream& out
		, const std::string& leading_str) const;

private:
	value_type deriv_;

};	// end class MeritFuncNLESqrResid

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_NLE_SQR_RESID_H