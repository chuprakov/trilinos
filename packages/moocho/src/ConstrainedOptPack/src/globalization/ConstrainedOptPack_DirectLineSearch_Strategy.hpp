// //////////////////////////////////////////////////////////////////////////////////
// DirectLineSearch_Strategy.h
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

#ifndef DIRECT_LINE_SEARCH_STRATEGY_H
#define DIRECT_LINE_SEARCH_STRATEGY_H

#include <stdexcept>
#include <iosfwd>

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** This is the interface for strategy objects that
  * perform a line search from an initial point along
  * a search direction given a merit function {abstract}.
  */
class DirectLineSearch_Strategy {
public:

	/** @name Exceptions */
	//@{

	/// Thrown if the direction vector d_k is not a descent direction for the merit funciton
	class NotDescentDirection : public std::logic_error
	{public: NotDescentDirection(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	///
	virtual ~DirectLineSearch_Strategy() {}

	/// Set the maximum number of iterations
	virtual void set_max_iter(int max_iter) = 0;

	/// Get the maximum number of iterations
	virtual int max_iter() const = 0;

	/// Get the number of iterations performed
	virtual int num_iterations() const = 0;

	///
	/** Called to perform the linesearch.
	  *
	  * This operaion  computes the
	  * approximate minimum to a merit function along a search direcation.
	  * More specifically the following problem is approximatly solved:\\
	  *
	  * min  phi(alpha)  s.t. alpha = [0, alpha_upper]\\
	  *
	  * Actually, if the initial alpha satisfys an internal descent requirement, then
	  * it will be choosen over smaller values of alpha that may result in a 
	  * greater reduction in the given merit funciton.
	  *
	  * If the maximum number of iterations is exceeded then the subclass will return
	  * false and will return the values of alpha_k, x_kp1, and phi_kp1 for the 
	  * lowest value of phi_kp1 found, and the last call to phi.value(x) will
	  * be this best x_kp1.
	  *
	  * Preconditions: \begin{itemize}
	  * \item #phi.deriv(d_k) < 0# (throw NotDescentDirection)
	  * \end{itemize}
	  *
	  * @param	phi		[I]		The merit function object that will compute #phi.value(alpha)#
	  *							and the descent derivative.
	  * @param	phi_k	[I]		The value of #phi.value(0)#.  Not computed internally
	  *							for the sake of efficency.
	  * @param	alpha_k	[I/O]	The initial #alpha_k# to try on input (usually 1).
	  *							On output #alpha_k# is the accepted value for a successful
	  *							line search, or it will be the alpha_k for the minimum phi
	  *							found for a line search failure.
	  * @param	phi_kp1	[I/O]	Merit function at new point.
	  *							On input it must be the value of #phi.value(alpha_k)#
	  *							and on output is set to #phi.value(alpha_k)#.
	  *	@parm	out		[O]		If != 0 then output is sent to this stream to record
	  *							the progress of the linesearch iterations.  The default
	  *							is zero.
	  *
	  * @return					true: Successful line search; false: Line search failure.
	  */
	virtual bool do_line_search( const MeritFuncCalc1D& phi, value_type phi_k
		, value_type* alpha_k, value_type* phi_kp1
		, std::ostream* out = 0 ) = 0;

	///
	/** Print the direct line search algorithm.
	  *
	  * The default does nothing.
	  */
	virtual void print_algorithm(std::ostream& out, const std::string& leading_str) const
	{}

};	// end class DirectLineSearch_Strategy

}	// end namespace ConstrainedOptimizationPack

#endif	// DIRECT_LINE_SEARCH_STRATEGY_H
