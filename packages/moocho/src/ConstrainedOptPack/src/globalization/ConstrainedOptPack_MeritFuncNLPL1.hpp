// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLPL1.h

#ifndef MERIT_FUNC_NLP_L1_H
#define MERIT_FUNC_NLP_L1_H

#include "MeritFuncNLP.h"
#include "MeritFuncNLPDirecDeriv.h"
#include "MeritFuncPenaltyParam.h"

namespace ConstrainedOptimizationPack {

///
/** The L1 merit function.
  *
  * phi(x) = f(x) + mu * norm(c(x),1)
  *
  * Dphi(x_k,d_k) = Gf_k' * d_k - mu * norm(c_k,1)
  *
  * Note that the definition of Dphi(x_k,d_k) assumes
  * that Gc_k'*d_k + c_k = 0.  In otherwords, d_k must
  * satisfiy the linearized equality constraints at
  * at x_k.
  *
  * Implicit copy constructor and assignment operators
  * are allowed.
  */
class MeritFuncNLPL1
	: public MeritFuncNLP
	, public MeritFuncNLPDirecDeriv
	, public MeritFuncPenaltyParam
{
public:

	/// Initializes deriv() = 0 and mu() = 0
	MeritFuncNLPL1();

	// ////////////////////////////////
	// Overridden from MeritFuncNLP

	///
	value_type value(value_type f, const VectorSlice& c) const;

	///
	value_type deriv() const;

	///
	void print_merit_func(std::ostream& out
		, const std::string& leading_str) const;

	// ////////////////////////////////
	// Overridden from MeritFuncNLPDirecDeriv

	///
	value_type calc_deriv( const VectorSlice& Gf_k, const VectorSlice& c_k
		, const VectorSlice& d_k );

	// ////////////////////////////////
	// Overridden from MeritFuncPenaltyParam

	///
	void mu(value_type mu);

	///
	value_type mu() const;

private:
	value_type deriv_;
	value_type mu_;

};	// end class MeritFuncNLPL1

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_NLP_L1_H
