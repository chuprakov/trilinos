// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLPModL1.h

#ifndef MERIT_FUNC_NLP_MOD_L1_H
#define MERIT_FUNC_NLP_MOD_L1_H

#include "MeritFuncNLP.h"
#include "MeritFuncNLPDirecDeriv.h"
#include "MeritFuncPenaltyParams.h"
#include "LinAlgPack/include/VectorClass.h"

namespace ConstrainedOptimizationPack {

///
/** The modified L1 merit function using different penatly parameters for each constriant.
  *
  * phi(x) = f) + sum( mu(j) * abs(c(j)), j = 1,...,m )
  *
  * Dphi(x_k,d_k) = Gf_k' * d_k - sum( mu(j) * abs(c(j)), j = 1,...,m )
  *
  * Note that the definition of Dphi(x_k,d_k) assumes
  * that Gc_k'*d_k + c_k = 0.  In otherwords, d_k must
  * satisfiy the linearized equality constraints at
  * at x_k.
  *
  * Implicit copy constructor and assignment operators
  * are allowed.
  */
class MeritFuncNLPModL1
	: public MeritFuncNLP
	, public MeritFuncNLPDirecDeriv
	, public MeritFuncPenaltyParams
{
public:

	/// Initializes deriv() = 0 and mu() = 0
	MeritFuncNLPModL1();

	// ////////////////////////////////
	// Overridden from MeritFuncNLP

	///
	value_type value(value_type f, const VectorSlice& c) const;

	///
	value_type deriv() const;

	///
	void print_merit_func(std::ostream& out
		, const std::string& leading_str) const;

	// /////////////////////////////////////////////
	// Overridden from MeritFuncNLPDirecDeriv

	///
	/** If the value n passed to resize(n) does not
	  * equal the size of the vector parameters then
	  * an exception #MeritFuncNLP::InvalidInitialization#
	  * will be thrown.
	  */
	value_type calc_deriv( const VectorSlice& Gf_k, const VectorSlice& c_k
		, const VectorSlice& d_k );

	// //////////////////////////////////////////
	// Overridden from MeritFuncPenaltyParams

	///
	void resize( size_type n );

	///
	VectorSlice mu();

	///
	const VectorSlice mu() const;

private:
	value_type	deriv_;
	Vector		mu_;

};	// end class MeritFuncNLPModL1

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_NLP_MOD_L1_H