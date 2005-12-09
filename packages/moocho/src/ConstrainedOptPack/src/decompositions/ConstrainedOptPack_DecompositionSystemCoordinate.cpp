// /////////////////////////////////////////////////////////////////////////
// DecompositionSystemCoordinate.cpp
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

#include <assert.h>

#include "ConstrainedOptPack_DecompositionSystemCoordinate.hpp"
#include "ConstrainedOptPack_MatrixIdentConcatStd.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpSubView.hpp"
#include "AbstractLinAlgPack_MatrixZero.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TestForException.hpp"

namespace ConstrainedOptPack {

DecompositionSystemCoordinate::DecompositionSystemCoordinate(
	const VectorSpace::space_ptr_t     &space_x
	,const VectorSpace::space_ptr_t    &space_c
	,const basis_sys_ptr_t             &basis_sys
	,const basis_sys_tester_ptr_t      &basis_sys_tester
	,EExplicitImplicit                 D_imp
	,EExplicitImplicit                 Uz_imp
	)
	:DecompositionSystemVarReductImp(
		space_x, space_c, basis_sys, basis_sys_tester
		,D_imp, Uz_imp )
{}

// Overridden from DecompositionSystem

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemCoordinate::factory_Y() const
{
	namespace rcp = MemMngPack;
	return Teuchos::rcp(
		new Teuchos::AbstractFactoryStd<MatrixOp,MatrixIdentConcatStd>()
		);
}

const DecompositionSystem::mat_nonsing_fcty_ptr_t
DecompositionSystemCoordinate::factory_R() const
{
	if( basis_sys().get() )
		return basis_sys()->factory_C();
	return Teuchos::null;
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemCoordinate::factory_Uy() const
{
	return Teuchos::rcp(	new Teuchos::AbstractFactoryStd<MatrixOp,MatrixOpSubView>() );
}

// Overridden from DecompositionSystemVarReductImp

DecompositionSystem::mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t
DecompositionSystemCoordinate::uninitialize_matrices(
	std::ostream                                       *out
	,EOutputLevel                                      olevel
	,MatrixOp                                          *Y
	,MatrixOpNonsing                                   *R
	,MatrixOp                                          *Uy
	) const
{
	namespace rcp = MemMngPack;
	using Teuchos::dyn_cast;

	//
	// Get pointers to concreate matrices
	//
	
	MatrixOpSubView
		*Uy_sv = Uy ? &dyn_cast<MatrixOpSubView>(*Uy) : NULL;			

	//
	// Only uninitialize the sub_view matrices
	//

	if(Uy_sv)
		Uy_sv->initialize(Teuchos::null);

	//
	// Return the basis matrix object R == C as a smart pointer that is
	// not owned.  This matrix object (if R != NULL) will get updated
	// by the basis_sys object as C in the method
	// 
	//     DecompositionSystemVarReductImp::update_decomp(...)
	//

	return Teuchos::rcp(R,false);

}

void DecompositionSystemCoordinate::initialize_matrices(
	std::ostream                                           *out
	,EOutputLevel                                          olevel
	,const mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t &C
	,const mat_fcty_ptr_t::element_type::obj_ptr_t         &D
	,MatrixOp                                              *Y
	,MatrixOpNonsing                                       *R
	,MatrixOp                                              *Uy
	,EMatRelations                                         mat_rel
	) const
{
	namespace rcp = MemMngPack;
	using Teuchos::dyn_cast;

	const size_type
		n = this->n(),
		r = this->r();
	const Range1D
		var_dep(1,r),
		var_indep(r+1,n);

	//
	// Get pointers to concreate matrices
	//
	
	MatrixIdentConcatStd
		*Y_coor = Y ? &dyn_cast<MatrixIdentConcatStd>(*Y) : NULL;
	MatrixOpSubView
		*Uy_sv = Uy ? &dyn_cast<MatrixOpSubView>(*Uy) : NULL;			

	//
	// Initialize the matrices
	//

	if( Y_coor && Y_coor->D_ptr().get() == NULL ) {
		Y_coor->initialize(
			space_x()                                         // space_cols
			,space_x()->sub_space(var_dep)->clone()           // space_rows
			,MatrixIdentConcatStd::BOTTOM                     // top_or_bottom
			,0.0                                              // alpha
			,Teuchos::rcp(
				new MatrixZero(
					space_x()->sub_space(var_indep)->clone()
					,space_x()->sub_space(var_dep)->clone()
					) )                                       // D_ptr
			,BLAS_Cpp::no_trans                               // D_trans
			);
	}
	assert(Uy_sv == NULL); // ToDo: Implement for undecomposed equalities

	// The R = C matrix object should already be updateded

}

void DecompositionSystemCoordinate::print_update_matrices(
	std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Coordinate decompositon Y, R and Uy matrices (class DecompositionSystemCoordinate)\n"
		<< L << "Y = [ I; 0 ] (using class MatrixIdentConcatStd with MatrixZero)\n"
		<< L << "R = Gc(var_dep,equ_decomp)' = C\n"
		<< L << "Uy = Gc(var_dep,equ_undecomp)'\n"
		;
}
		
}	// end namespace ConstrainedOptPack
