// //////////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncCalcNLE.hpp
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

#ifndef MERIT_FUNC_CALC_NLE_H
#define MERIT_FUNC_CALC_NLE_H

#include "ConstrainedOptPack_MeritFuncCalc.hpp"
#include "ConstrainedOptPack_MeritFuncNLE.hpp"
#include "StandardAggregationMacros.hpp"

namespace ConstrainedOptPack {

///
/** Adds the ability to compute phi(c(x)) at x
  * directly instead of having to compute c first.
  * This class uses an aggregate NLP to perform the computations of c(x).
  */
class MeritFuncCalcNLE : public MeritFuncCalc {
public:

	/// <<std aggr>> stereotype members for phi.
	STANDARD_CONST_AGGREGATION_MEMBERS( MeritFuncNLE, phi )

	/// <<std aggr>> stereotype members for nlp.
	STANDARD_CONST_AGGREGATION_MEMBERS( NLP, nlp )

	///
	MeritFuncCalcNLE( const MeritFuncNLE* phi = 0, const NLP* nlp = 0 );

	// ////////////////////////////////////////////
	// Overridden from MeritFuncCalc

	///
	/** Return the value of the merit function at x.
	  * Here phi(x) is calculated directly using the nlp.
	  */
	value_type operator()(const Vector& x) const;

	/// Calls phi().deriv() on phi.
	value_type deriv() const;

	/// Calls phi().print_merit_func(....).
	void print_merit_func(std::ostream& out
		, const std::string& leading_str) const;

};	// end class MeritFuncCalcNLE

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_NLE_H
