// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLPDirecDeriv.h

#ifndef MERIT_FUNC_NLP_DIREC_DERIV_H
#define MERIT_FUNC_NLP_DIREC_DERIV_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** This class provides the interface for allowing subclass merit
  * functions to compute the directional 1D derivative at a base point
  * {abstract}.
  *
  * The quantities Gf(xo) (gradient of f(xo))
  * c(xo) and d are used by several
  * types of merit functions to calculate the derivative of:\\
  * d(phi(x_k + alpha_k*d_k))/d(alpha_k) at alpha_k = 0.
  *
  * It is generally assumed that d satisfies Gc_k'*d_k + c_k = 0 otherwise the
  * merit function would need Gc_k to compute this directional derivative
  * properly.
  */
class MeritFuncNLPDirecDeriv {
public:

	///
	virtual ~MeritFuncNLPDirecDeriv() {}

	///
	/** Calculate d(phi(x_k + alpha_k*d_k))/d(alpha_k) at alpha_k = 0.
	  *
	  * The value is stored internally by the subclass are returned by its
	  * deriv() member usually.  The value is also returned from this
	  * function.
	  *
	  * If the sizes of the vectors input do not aggree then
	  * #std::length_error# exception will be thrown.
	  */
	virtual value_type calc_deriv( const VectorSlice& Gf_k, const VectorSlice& c_k
		, const VectorSlice& d_k ) = 0;

};	// end class MeritFuncBase

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_NLP_DIREC_DERIV_H