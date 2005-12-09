// //////////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncCalc1D.hpp
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

#ifndef MERIT_FUNC_CALC_1D_H
#define MERIT_FUNC_CALC_1D_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

///
/** Abstracts a 1D merit function {abstract}.
  *
  * This is the interface that line search algorithms use to compute
  * the value of the merit function at alpha (phi(alpha)) and
  * to retrieve the initial descent derivative of the merit function
  * (using \c deriv()).
  */
class MeritFuncCalc1D {
public:

	///
	virtual ~MeritFuncCalc1D() {}

	/// Return the value of the merit function at alpha.
	virtual value_type operator()( value_type alpha ) const = 0;

	/// Return the derivative of the merit function at alpha = 0
	virtual value_type deriv() const = 0;

	/// Print the particular merit function
	virtual void print_merit_func(
		std::ostream& out, const std::string& leading_str ) const = 0;

};	// end class MeritFuncCalc1D

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_1D_H
