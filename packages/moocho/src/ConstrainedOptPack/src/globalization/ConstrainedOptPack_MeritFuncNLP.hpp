// /////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncNLP.hpp
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

#ifndef MERIT_FUNC_NLP_H
#define MERIT_FUNC_NLP_H

#include <iosfwd>

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

///
/** Base class for all merit functions for NonLinear Programs (NLP) {abstract}.
  */
class MeritFuncNLP {
public:

	///
	class InvalidInitialization : public std::logic_error
	{public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	virtual ~MeritFuncNLP() {}

	/** @name To be overridden by subclasses */
	//@{

	///
	/** Assign the state of one Merit functions.
	 *
	 * The default implementation throws an <tt>std::logic_error</tt> exception
	 * unless it is assignment to self.
	 */
	virtual MeritFuncNLP& operator=(const MeritFuncNLP&);

	///
	/** Return the value of the merit function at f(x), c(x), h(x).
	 * This interface requires the client to compute f(x)
	 * c(x) and h(x) and pass it to this function to have
	 * the value of phi(f,c,h,hl,hu) calculated.
	 *
	 * If the merit function has not been initialized properly
	 * then a <tt>InvalidInitialization</tt> exception will be thrown.
	 */
	virtual value_type value(
		value_type             f
		,const Vector    *c
		,const Vector    *h
		,const Vector    *hl
		,const Vector    *hu
		) const = 0;

	///
	/** Return the value of the directional derivative of the 
	 * merit function w.r.t. alpha at alpha = 0.  In other words
	 * compute return d( phi(f(x),c(x),h(x)) ) / d(alpha_k) at alpha_k = 0
	 * where x = x_k + alpha_k * d_k.
	 *
	 * If the merit function has not been initialized properly
	 * then a <tt>InvalidInitialization</tt> exception will be thrown.
	 */
	virtual value_type deriv() const = 0;

	///
	/** Print the merit funciton
	 */
	virtual void print_merit_func(
		std::ostream& out, const std::string& leading_str ) const = 0;

	//@}

};	// end class MeritFuncNLP

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLP_H
