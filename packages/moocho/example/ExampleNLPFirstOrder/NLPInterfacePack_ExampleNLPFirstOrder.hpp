// ///////////////////////////////////////////////////////////////////
// ExampleNLPFirstOrderInfo.hpp
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

#ifndef EXAMPLE_NLP_FIRST_ORDER_INFO_H
#define EXAMPLE_NLP_FIRST_ORDER_INFO_H

#include "ExampleNLPFirstOrderDirect/ExampleNLPObjGradient.hpp"
#include "NLPInterfacePack/src/NLPFirstOrderInfo.hpp"

namespace NLPInterfacePack {

///
/** Simple example %NLP subclass to illustrate how to implement the
 * \c NLPFirstOrderInfo interface for a specialized \c NLP.
 *
 * The example %NLP we will use is a scalable problem where
 * the basis of the jacobian of the constraints is a diagonal
 * matrix.
 \verbatim

    min    f(x) = (1/2) * sum( x(i)^2, for i = 1..n )
    s.t.   c(x)(j) = x(j) * (x(m+j) -1) - 10 * x(m+j) = 0, for j = 1..m
          0.01 < x(i) < 20, for i = p...p+m

    where:
        m = n/2
        p = 1 if dep_bounded == true or m+1 if dep_bounded = false
 \endverbatim
 * This subclass inherits from the subclass
 * <tt>\ref NLPInterfacePack::ExampleNLPFirstOrderDirect "ExampleNLPFirstOrderDirect"</tt>
 * mostly out of lazyness but also to show how flexible these interfaces
 * can be using mutiple inheritance.
 *
 * ToDo: Finish documentation!
 */
class ExampleNLPFirstOrderInfo
	: virtual public NLPFirstOrderInfo
	, virtual public ExampleNLPObjGradient
{
public:

	///
	/** Constructor (see </tt>ExampleNLPFirstOrderDirect::ExampleNLPFirstOrderDirect()</tt>).
	 */
	ExampleNLPFirstOrderInfo(
		const VectorSpace::space_ptr_t&  vec_space
		,value_type                      xo
		,bool                            has_bounds
		,bool                            dep_bounded
		);

	/** @name Overridden public members from NLP */
	//@{

	///
	void initialize(bool test_setup);
	///
	bool is_initialized() const;

	//@}

	/** @name Overridden public members from NLPFirstOrderInfo */
	//@{

	/// Overridden to check the concrete type of Gc
	void set_Gc(MatrixWithOp* Gc);
	///
	const NLPFirstOrderInfo::mat_fcty_ptr_t factory_Gc() const;
	/// Return NULL
	const NLPFirstOrderInfo::mat_fcty_ptr_t factory_Gh() const;
	/// Returns an ExampleBasisSystem
	const basis_sys_ptr_t basis_sys() const;

	//@}

protected:

	/** @name Overridden protected members from NLPFirstOrderInfo */
	//@{

	///
	void imp_calc_Gc(const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const;
	///
	void imp_calc_Gh(const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const;

	//@}

private:

	// /////////////////////////////////////////
	// Private data members

	bool                                    initialized_;  // flag for if initialized has been called.
	NLPFirstOrderInfo::mat_fcty_ptr_t       factory_Gc_;   // Factory for Gc
	NLPFirstOrderInfo::basis_sys_ptr_t      basis_sys_;    // The basis system

	// /////////////////////////////////////////
	// Private member functions

	///
	void assert_is_initialized() const;

};	// end class ExampleNLPFirstOrderInfo

// ///////////////////////////////////////////////
// Inline member functions

inline
void ExampleNLPFirstOrderInfo::assert_is_initialized() const
{
    using NLPInterfacePack::NLP;
	if( !is_initialized() )
		throw NLP::UnInitialized("ExampleNLPFirstOrderInfo::assert_is_initialized() : Error, "
			"ExampleNLPFirstOrderInfo::initialize() has not been called yet." );
}

}	// end namespace NLPInterfacePack

#endif	// EXAMPLE_NLP_FIRST_ORDER_INFO_H
