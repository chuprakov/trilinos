// /////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncNLPL1.hpp
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

#ifndef MERIT_FUNC_NLP_L1_H
#define MERIT_FUNC_NLP_L1_H

#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp"
#include "ConstrainedOptPack_MeritFuncPenaltyParam.hpp"

namespace ConstrainedOptPack {

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
 * ToDo: Add a term for general inequalities hl <= h <= hu
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

	/** @name Overridden from MeritFuncNLP */
	//@{

	///
	MeritFuncNLP& operator=(const MeritFuncNLP&);
	///
	value_type value(
		value_type             f
		,const Vector    *c
		,const Vector    *h
		,const Vector    *hl
		,const Vector    *hu
		) const;
	///
	value_type deriv() const;
	///
	void print_merit_func(
		std::ostream& out, const std::string& leading_str ) const;

	//@}

	/** @name Overridden from MeritFuncNLPDirecDeriv */
	//@{

	///
	value_type calc_deriv(
		const Vector    &Gf_k
		,const Vector   *c_k
		,const Vector   *h_k
		,const Vector   *hl
		,const Vector   *hu
		,const Vector   &d_k
		);

	//@}

	/** @name Overridden from MeritFuncPenaltyParam */
	//@{

	///
	void mu(value_type mu);

	///
	value_type mu() const;

	//@}

private:
	value_type deriv_;
	value_type mu_;

};	// end class MeritFuncNLPL1

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLP_L1_H
