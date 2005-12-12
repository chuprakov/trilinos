// //////////////////////////////////////////////////////////////////////////
// NLPInterfacePack_ExampleBasisSystem.hpp
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

#include "NLPInterfacePack_Types.hpp"
#include "AbstractLinAlgPack_BasisSystemComposite.hpp"

namespace NLPInterfacePack {

///
/** Subclass of BasisSystem for example NLP.
 *
 * ToDo: Finish documentation!
 */
class ExampleBasisSystem
	: public AbstractLinAlgPack::BasisSystemComposite
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
	 *
	 * ToDo: Finish documentation!
	 */
	void initialize(
		const VectorSpace::space_ptr_t       &space_x
		,const Range1D                       &var_dep
		,const Range1D                       &var_indep
		);

	/** @name Overridden from BasisSystemComposite */
	//@{

	///
	void update_D(
		const MatrixOpNonsing       &C
		,const MatrixOp             &N
		,MatrixOp                   *D
		,EMatRelations               mat_rel
		) const;

	//@}

private:
	// Not defined and not to be called!
	ExampleBasisSystem();

}; // end class ExampleBasisSystem

} // end namespace NLPInterfacePack

#endif // EXAMPLE_BASIS_SYSTEM_H
