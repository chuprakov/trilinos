// //////////////////////////////////////////
// ExampleNLPObjGradient.h
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

#ifndef EXAMPLE_NLP_OBJ_GRADIENT_H
#define EXAMPLE_NLP_OBJ_GRADIENT_H

#include "NLPInterfacePack/src/NLPObjGradient.h"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/src/VectorSpace.h"
#include "AbstractLinAlgPack/src/VectorSpaceCompositeStd.h"

namespace NLPInterfacePack {

///
/** Simple example %NLP subclass to illustrate how to implement the
 * \c NLPObjGradient interface for a specialized \c NLP.
 *
 * The example %NLP we will use is a scalable problem where
 * the basis of the jacobian of the constraints is a diagonal
 * matrix (however it is not computed here).
 \verbatim

    min    f(x) = (1/2) * sum( x(i)^2, for i = 1..n )
    s.t.   c(x)(j) = x(j) * (x(m+j) -1) - 10 * x(m+j) = 0, for j = 1..m
          0.01 < x(i) < 20, for i = p...p+m

    where:
        m = n/2
        p = 1 if dep_bounded == true or m+1 if dep_bounded = false
 \endverbatim
 * This is not really a fully functional NLP in the sense that there is
 * no derivative information for the constraints.
 */
class ExampleNLPObjGradient
	: virtual public NLPObjGradient
{
public:

	///
	/** Constructor.
	 *
	 * @param  vec_space  [in] Smart pointer to a vector space object that will
	 *                    be used to define the spaces of dependent and independent
	 *                    variables.
	 * @param  xo         [in] The initial starting guess for \a x.
	 * @param  has_bounds [in] If \c true, then the NLP will have bounds.  If \c false
	 *                    then it will not have bounds.
	 * @param  dep_bouned [in] If \c true, then the bounds will be on the dependent
	 *                    variables.  If \c false, then the bounds will be on the
	 *                    independent variable.  This argument is ignored if
	 *                    <tt>has_bounds == false</tt>.
	 */
	ExampleNLPObjGradient(
		const VectorSpace::space_ptr_t&  vec_space
		,value_type                      xo
		,bool                            has_bounds
		,bool                            dep_bounded
		);

	/** @name Helper methods to be used by subclasses. */
	//@{

	///
	virtual Range1D var_dep() const;
	///
	virtual Range1D var_indep() const;

	//@}

	/** @name Overridden public members from NLP */
	//@{

	///
	void initialize(bool test_setup);
	///
	bool is_initialized() const;
	///
	size_type n() const;
	///
	size_type m() const;
	///
	vec_space_ptr_t space_x() const;
	///
	vec_space_ptr_t space_c() const;
	/// Returns <tt>return.get() == NULL</tt>.
	vec_space_ptr_t space_h() const;
	///
    size_type num_bounded_x() const;
	///
	void force_xinit_in_bounds(bool force_xinit_in_bounds);
	///
	bool force_xinit_in_bounds() const;
	///
	const VectorWithOp& xinit() const;
	///
	const VectorWithOp& xl() const;
	///
	const VectorWithOp& xu() const;
	///
	value_type max_var_bounds_viol() const;
	/// Throws exception.
	const VectorWithOp& hl() const;
	/// Throws exception.
	const VectorWithOp& hu() const;
	///
	void scale_f( value_type scale_f );
	///
	value_type scale_f() const;
	///
	void report_final_solution(
		const VectorWithOp&    x
		,const VectorWithOp*   lambda
		,const VectorWithOp*   lambdaI
		,const VectorWithOp*   nu
		,bool                  optimal
		) const;

	//@}

protected:

	/** @name Overridden protected members from NLP */
	//@{

	///
	void imp_calc_f(
		const VectorWithOp& x, bool newx
		,const ZeroOrderInfo& zero_order_info) const;
	///
	void imp_calc_c(
		const VectorWithOp& x, bool newx
		,const ZeroOrderInfo& zero_order_info) const;
	/// This implementation does nothing (should never be called though).
	void imp_calc_h(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const;

	//@}

	/** @name Overridden protected members from NLPObjGradient */
	//@{

	///
	void imp_calc_Gf(
		const VectorWithOp& x, bool newx
		,const ObjGradInfo& obj_grad_info) const;

	//@}

private:

	// /////////////////////////////////////////
	// Private data members

	VectorSpace::space_ptr_t    vec_space_;       // The vector space for dependent and indepenent variables and c(x).
	VectorSpace::space_ptr_t    vec_space_comp_;  // Composite vector space for x = [ xD; xI ]
	Range1D                     var_dep_;         // Range for dependnet variables.
	Range1D                     var_indep_;       // Range for independent variables.

	bool         initialized_;            // flag for if initialized has been called.
	value_type   obj_scale_;              // default = 1.0;
	bool         has_bounds_;             // default = true
	bool         force_xinit_in_bounds_;  // default = true.

	size_type    n_;                      // Number of variables in the problem.
	VectorSpace::vec_mut_ptr_t  xinit_;   // Initial guess.
	VectorSpace::vec_mut_ptr_t  xl_;      // lower bounds.
	VectorSpace::vec_mut_ptr_t  xu_;      // upper bounds.

	// /////////////////////////////////////////
	// Private member functions

	///
	void assert_is_initialized() const;

};	// end class ExampleNLPObjGradient

// ///////////////////////////////////////////////
// Inline member functions

inline
void ExampleNLPObjGradient::assert_is_initialized() const
{
    using NLPInterfacePack::NLP;
	if( !is_initialized() )
		throw NLP::UnInitialized("ExampleNLPObjGradient::assert_is_initialized() : Error, "
			"ExampleNLPObjGradient::initialize() has not been called yet." );
}

}	// end namespace NLPInterfacePack

#endif	// EXAMPLE_NLP_OBJ_GRADIENT_H
