// //////////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MeritFuncCalc1DQuadratic.hpp
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

#ifndef MERIT_FUNC_CALC_1D_QUADRATIC_H
#define MERIT_FUNC_CALC_1D_QUADRATIC_H

#include "ConstrainedOptPack_MeritFuncCalc1D.hpp"
#include "ConstrainedOptPack_MeritFuncCalc.hpp"

namespace ConstrainedOptPack {

///
/** Adds the ability to compute phi(alpha) at alpha of a given set of vectors.
  *
  * Computes <tt>phi( x = sum( alpha^k * d[k], k = 0...p-1 ) )</tt> where 
  * <tt>1 <= p <= 2</tt>.
  */
class MeritFuncCalc1DQuadratic : public MeritFuncCalc1D {
public:

	///
	typedef const Vector*   const_VectorWithOp_ptr;

	/** @name Constructors */
	//@{

	///
	/** The only constructor.
	  *
	  * Note that \c *x and \c *d gets updated as <tt>operator()(alpha)</tt> is called.
	  *
	  * The client must ensure that the memory pointed to by the vectors in d must not
	  * be desturbed while this object is in use.  To do so may have bad side effects.
	  *
	  * @param  phi  [in] The merit function to use.
	  * @param  p    [in] The number of vectors in \c d[].
	  * @param  d    [in] Array (length \c p) of pointers to the rhs \c d[] vectors.
	  * @param  x    [out] The vector that gets updated.
	  *
	  * Preconditions:<ul>
	  * <li> <tt>1 <= p <= 3<tt>
	  * <li> <tt>d[k]->space().is_compatible(x->space()), for k = 0...p-1</tt>
	  * </ul>
	  */
	MeritFuncCalc1DQuadratic(
		const MeritFuncCalc&      phi
		,size_type                p
		,const_VectorWithOp_ptr   d[]
		,VectorMutable*     x
		);

	//@}

	/** @name Overridden from MeritFuncCalc1D */
	//@{

	/// Returns <tt>phi( x = sum( alpha^k * d[k], k = 0...p-1 ) )</tt>.
	value_type operator()(value_type alpha) const;

	/// Returns phi.deriv()
	value_type deriv() const;

	/// Calls <tt>phi->print_merit_func()</tt>.
	void print_merit_func(
		std::ostream& out, const std::string& leading_str ) const;

	//@}

private:
	const MeritFuncCalc&        phi_;
	size_type                   p_;
	const_VectorWithOp_ptr      d_[3];
	VectorMutable         *x_;

	// not defined and not to be called
	MeritFuncCalc1DQuadratic();
	MeritFuncCalc1DQuadratic& operator=( const MeritFuncCalc1DQuadratic& );

};	// end class MeritFuncCalc1DQuadratic

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_1D_QUADRATIC_H
