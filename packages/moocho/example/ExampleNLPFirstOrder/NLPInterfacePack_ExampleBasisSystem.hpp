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
	ExampleBasisSystem( const VectorSpace::space_ptr_t& space_x_DI = NULL );
	
	///
	/** Initialize given the vector space for the dependent and independent variables.
	 *
	 * @param  space_x_DI
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->var_dep().size()   == space_x_DI->dim()</tt>
	 * <li><tt>this->var_indep().size() == space_x_DI->dim()</tt>
	 * <li><tt>this->space_C()->space_rows().is_compatible(*space_x_DI) == true</tt>
	 * <li><tt>this->space_C()->space_cols().is_compatible(*space_x_DI) == true</tt>
	 * <li><tt>this->space_D()->space_rows().is_compatible(*space_x_DI) == true</tt>
	 * <li><tt>this->space_D()->space_cols().is_compatible(*space_x_DI) == true</tt>
	 * </ul>
	 */
	void initialize( const VectorSpace::space_ptr_t& space_x_DI );

	/** @name Overridden from BasisSystem */
	//@{

	///
	const space_C_ptr_t& space_C() const;
	///
	const space_D_ptr_t& space_D() const;
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
		);

	//@}

private:
	
	VectorSpace::space_ptr_t   space_x_DI_;
	Range1D                    var_dep_,
	                           var_indep_;
	space_C_ptr_t              space_C_;
	space_D_ptr_t              space_D_;

}; // end class ExampleBasisSystem

} // end namespace NLPInterfacePack

#endif // EXAMPLE_BASIS_SYSTEM_H
