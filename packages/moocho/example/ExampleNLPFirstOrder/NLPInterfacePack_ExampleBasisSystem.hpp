// //////////////////////////////////////////////////////////////////////////
// ExampleBasisSystem.h
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

#ifndef EXAMPLE_BASIS_SYSTEM_H
#define EXAMPLE_BASIS_SYSTEM_H

#include "NLPInterfacePack/include/NLPInterfacePackTypes.h"
#include "AbstractLinAlgPack/include/BasisSystem.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "Range1D.h"

namespace NLPInterfacePack {

///
/** Subclass of BasisSystem for example NLP.
 *
 * ToDo: Finish documentation!
 */
class ExampleBasisSystem : public BasisSystem {
public:

	/// Calls <tt>this->initialize()</tt>
	ExampleBasisSystem( const VectorSpace::space_ptr_t& space_x_DI = ReferenceCountingPack::null );
	
	///
	/** Initialize given the vector space for the dependent and independent variables.
	 *
	 * @param  space_x_DI
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->var_dep().size()   == space_x_DI->dim()</tt>
	 * <li><tt>this->var_indep().size() == space_x_DI->dim()</tt>
	 * </ul>
	 */
	void initialize( const VectorSpace::space_ptr_t& space_x_DI );

	/** @name Overridden from BasisSystem */
	//@{

	///
	const mat_nonsing_fcty_ptr_t factory_C() const;
	///
	const mat_fcty_ptr_t factory_D() const;
	///
	Range1D var_dep() const;
	///
	Range1D var_indep() const;
	///
	void update_basis(
		const MatrixWithOp*         Gc
		,const MatrixWithOp*        Gh
		,MatrixWithOpNonsingular*   C
		,MatrixWithOp*              D
		,MatrixWithOp*              GcUP
		,MatrixWithOp*              GhUP
		,EMatRelations              mat_rel
		) const;

	//@}

private:
	
	VectorSpace::space_ptr_t   space_x_DI_;
	Range1D                    var_dep_,
	                           var_indep_;
	mat_nonsing_fcty_ptr_t     factory_C_;
	mat_fcty_ptr_t             factory_D_;

}; // end class ExampleBasisSystem

} // end namespace NLPInterfacePack

#endif // EXAMPLE_BASIS_SYSTEM_H
