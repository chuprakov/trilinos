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
#include "AbstractLinAlgPack/include/BasisSystemCompositeStd.h"

namespace NLPInterfacePack {

///
/** Subclass of BasisSystem for example NLP.
 *
 * ToDo: Finish documentation!
 */
class ExampleBasisSystem
	: public AbstractLinAlgPack::BasisSystemCompositeStd
{
public:

	/// Calls <tt>this->initialize()</tt>
	ExampleBasisSystem(
		const VectorSpace::space_ptr_t       &space_x
		,const Range1D                       &var_dep
		,const Range1D                       &var_indep
		);
	
	///
	/** Initialize given the vector space for the dependent and independent variables.
	 *
	 * @param  space_x   [in]
	 * @param  var_dep   [in]
	 * @param  var_indep [in]
	 */
	void initialize(
		const VectorSpace::space_ptr_t       &space_x
		,const Range1D                       &var_dep
		,const Range1D                       &var_indep
		);

	/** @name Overridden from BasisSystemCompositeStd */
	//@{

	///
	void update_D(
		const MatrixWithOpNonsingular&  C
		,const MatrixWithOp&            N
		,MatrixWithOp*                  D
		,EMatRelations                  mat_rel
		) const;

	//@}

}; // end class ExampleBasisSystem

} // end namespace NLPInterfacePack

#endif // EXAMPLE_BASIS_SYSTEM_H
