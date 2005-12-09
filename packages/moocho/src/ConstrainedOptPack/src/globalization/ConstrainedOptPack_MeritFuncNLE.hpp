// //////////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncNLE.hpp
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

#ifndef MERIT_FUNC_NLE_H
#define MERIT_FUNC_NLE_H

#include <iosfwd>

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

///
/** Base class for all merit functions for systems of NonLinear Equations (NLE) {abstract}.
  */
class MeritFuncNLE {
public:

	///
	class InvalidInitialization : public std::logic_error
	{public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	virtual ~MeritFuncNLE() {}

	///
	/** Return the value of the merit function at <tt>c(x)</tt>.
	 * This interface requires the client to compute
	 * <tt>c(x)</tt> and pass it to this function to have
	 * the value of phi(c) calculated.
	 *
	 * If the merit function has not been initialized properly
	 * then a <tt>InvalidInitialization</tt> exception will be thrown.
	 */
	virtual value_type value(const Vector& c) const = 0;

	///
	/** Return the value of the directional derivative of the 
	 * merit function w.r.t. <tt>alpha</tt> at <tt>alpha = 0</tt>.  In other words
	 * compute return <tt>d(phi(c(x)))/d(alpha_k)</tt> at <tt>alpha_k = 0</tt>
	 * where <tt>x = x_k + alpha_k * d_k</tt>.
	 *
	 * If the merit function has not been initialized properly
	 * then a <tt>InvalidInitialization</tt> exception will be thrown.
	 */
	virtual value_type deriv() const = 0;

	///
	/** Print the merit funciton
	  */
	virtual void print_merit_func(std::ostream& out
		, const std::string& leading_str) const = 0;

};	// end class MeritFuncNLE

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLE_H
