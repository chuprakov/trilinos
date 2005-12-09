// /////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncCalc.hpp
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

#ifndef MERIT_FUNC_CALC_H
#define MERIT_FUNC_CALC_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

///
/** Abstract iterface for n-D merit functions {abstract}.
  *
  * Used to compute the value of the merit at a point \a x (\c phi(x)) and
  * to retrieve the derivative (\c phi.deriv()) along some direction \a d from
  * some base point \a xo.
  */
class MeritFuncCalc  {
public:

	///
	virtual ~MeritFuncCalc() {}

	///
	/** Return the value of the merit function at \a x.
	  */
	virtual value_type operator()(const Vector& x) const= 0;

	/// Calls value(d_k) on aggregate merit_func.
	virtual value_type deriv() const = 0;

	/// Print what this merit function is
	virtual void print_merit_func(
		std::ostream& out, const std::string& leading_str ) const = 0;

};	// end class MeritFuncCalc

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_H
