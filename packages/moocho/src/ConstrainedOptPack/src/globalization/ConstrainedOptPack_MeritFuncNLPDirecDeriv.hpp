// /////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef MERIT_FUNC_NLP_DIREC_DERIV_H
#define MERIT_FUNC_NLP_DIREC_DERIV_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

///
/** This class provides a mix-in interface for allowing subclass merit
 * functions to compute the directional 1D derivative at a base point.
 *
 * The quantities Gf(xo) (gradient of f(xo))
 * c(xo), h(xo) and d are used by several
 * types of merit functions to calculate the derivative of:<br>
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

	/** @name To be overridden by subclasses */
	//@{

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
	virtual value_type calc_deriv(
		const Vector    &Gf_k
		,const Vector   *c_k
		,const Vector   *h_k
		,const Vector   *hl
		,const Vector   *hu
		,const Vector   &d_k
		) = 0;

	//@}

};	// end class MeritFuncNLPDirecDeriv

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLP_DIREC_DERIV_H
