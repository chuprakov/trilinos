// /////////////////////////////////////////////////////////////////////////
// DecompositionSystemOrthogonal.cpp
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

#include <typeinfo>

#include "ConstrainedOptimizationPack/include/DecompositionSystemOrthogonal.h"
#include "ConstrainedOptimizationPack/include/MatrixIdentConcatStd.h"
#include "ConstrainedOptimizationPack/include/MatrixDecompRangeOrthog.h"
#include "ConstrainedOptimizationPack/include/VarReductOrthog_Strategy.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "AbstractFactoryStd.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace ConstrainedOptimizationPack {

DecompositionSystemOrthogonal::DecompositionSystemOrthogonal(
	const VectorSpace::space_ptr_t           &space_x
	,const VectorSpace::space_ptr_t          &space_c
	,const VectorSpace::space_ptr_t          &space_h
	,const basis_sys_ptr_t                   &basis_sys
	,const basis_sys_tester_ptr_t            &basis_sys_tester
	,const var_reduct_orthog_strategy_ptr_t  &var_reduct_orthog_strategy
	,EExplicitImplicit                       D_imp
	,EExplicitImplicit                       Uz_imp
	,EExplicitImplicit                       Vz_imp
	)
	:DecompositionSystemVarReduct(
		space_x, space_c, space_h, basis_sys, basis_sys_tester
		,D_imp,Uz_imp,Vz_imp )
	 ,var_reduct_orthog_strategy_(var_reduct_orthog_strategy)
{}

// Overridden from DecompositionSystem

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemOrthogonal::factory_Y() const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp(
		new AbstractFactoryPack::AbstractFactoryStd<MatrixWithOp,MatrixIdentConcatStd>()
		);
}

const DecompositionSystem::mat_nonsing_fcty_ptr_t
DecompositionSystemOrthogonal::factory_R() const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp(
		new AbstractFactoryPack::AbstractFactoryStd<MatrixWithOpNonsingular,MatrixDecompRangeOrthog>()
		);
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemOrthogonal::factory_Uy() const
{
	assert(0); // ToDo: Return MatrixCompositeStd (which will use MatrixProductStd)
	return ReferenceCountingPack::null;
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemOrthogonal::factory_Vy() const
{
	assert(0); // ToDo: Return MatrixCompositeStd (which will use MatrixProductStd)
	return ReferenceCountingPack::null;
}

// Overridden from DecompositionSystemVarReduct

DecompositionSystem::mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t
DecompositionSystemOrthogonal::uninitialize_matrices(
	std::ostream                                           *out
	,EOutputLevel                                          olevel
	,MatrixWithOp                                          *Y
	,MatrixWithOpNonsingular                               *R
	,MatrixWithOp                                          *Uy
	,MatrixWithOp                                          *Vy
	) const
{
	namespace rcp = ReferenceCountingPack;
	using DynamicCastHelperPack::dyn_cast;
	typedef DecompositionSystem::mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t
		C_ptr_t;

	//
	// Get pointers to concreate matrices
	//
	
	MatrixIdentConcatStd
		*Y_orth = Y ? &dyn_cast<MatrixIdentConcatStd>(*Y)    : NULL;
	MatrixDecompRangeOrthog
		*R_orth = R ? &dyn_cast<MatrixDecompRangeOrthog>(*R) : NULL;
	MatrixCompositeStd
		*Uy_cpst = Uy ? &dyn_cast<MatrixCompositeStd>(*Uy)      : NULL;			
	MatrixCompositeStd
		*Vy_cpst = Vy ? &dyn_cast<MatrixCompositeStd>(*Vy)      : NULL;

	//
	// Get the smart pointer to the basis matrix object C and the
	// matrix S = I + D'*D
	//
	
	C_ptr_t C_ptr = rcp::null;
	if(R_orth) {
		C_ptr  = rcp::rcp_const_cast<MatrixWithOpNonsingular>(    R_orth->C_ptr() ); // This could be NULL!
		S_ptr_ = rcp::rcp_const_cast<MatrixSymWithOpNonsingular>( R_orth->S_ptr() ); // ""
	}
	
	//
	// Uninitialize all of the matrices to remove references to C, D etc.
	//

	if(Y_orth)
		Y_orth->set_uninitialized();
	if(R_orth)
		R_orth->set_uninitialized();
	if(Uy_cpst)
		Uy_cpst->reinitialize();
	if(Vy_cpst)
		Vy_cpst->reinitialize();

	//
	// Return the owned? basis matrix object C
	//

	return C_ptr;

}

void DecompositionSystemOrthogonal::initialize_matrices(
	std::ostream                                           *out
	,EOutputLevel                                          olevel
	,const mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t &C
	,const mat_fcty_ptr_t::element_type::obj_ptr_t         &D
	,MatrixWithOp                                          *Y
	,MatrixWithOpNonsingular                               *R
	,MatrixWithOp                                          *Uy
	,MatrixWithOp                                          *Vy
	,EMatRelations                                         mat_rel
	) const
{
	namespace rcp = ReferenceCountingPack;
	using DynamicCastHelperPack::dyn_cast;
	typedef DecompositionSystem::mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t
		C_ptr_t;
	typedef DecompositionSystem::mat_fcty_ptr_t::element_type::obj_ptr_t
		D_ptr_t;

	const size_type
		n = this->n(),
		m = this->m(),
		r = this->r();
	const Range1D
		var_dep(1,r),
		var_indep(r+1,n),
		con_decomp   = this->con_decomp(),
		con_undecomp = this->con_undecomp();

	//
	// Get pointers to concreate matrices
	//
	
	MatrixIdentConcatStd
		*Y_orth = Y ? &dyn_cast<MatrixIdentConcatStd>(*Y)    : NULL;
	MatrixDecompRangeOrthog
		*R_orth = R ? &dyn_cast<MatrixDecompRangeOrthog>(*R) : NULL;
	MatrixCompositeStd
		*Uy_cpst = Uy ? &dyn_cast<MatrixCompositeStd>(*Uy)      : NULL;			
	MatrixCompositeStd
		*Vy_cpst = Vy ? &dyn_cast<MatrixCompositeStd>(*Vy)      : NULL;

	//
	// Initialize the matrices
	//

	if(Y_orth) {
		D_ptr_t  D_ptr = D;
		if(mat_rel == MATRICES_INDEP_IMPS) {
			D_ptr = D->clone();
			THROW_EXCEPTION(
				D_ptr.get() == NULL, std::logic_error
				,"DecompositionSystemOrthogonal::update_decomp(...) : Error, "
				"The matrix class used for the direct sensitivity matrix D = inv(C)*N of type \'"
				<< typeid(*D).name() << "\' must return return.get() != NULL from the clone() method "
				"since mat_rel == MATRICES_INDEP_IMPS!" );
		}
		Y_orth->initialize(
			space_x()                                         // space_cols
			,space_x()->sub_space(var_dep)->clone()           // space_rows
			,MatrixIdentConcatStd::BOTTOM                     // top_or_bottom
			,-1.0                                             // alpha
			,D_ptr                                            // D_ptr
			,BLAS_Cpp::no_trans                               // D_trans
			);
	}
	if(R_orth) {
		C_ptr_t  C_ptr = C;
		if(mat_rel == MATRICES_INDEP_IMPS) {
			C_ptr = C->clone_mwons();
			THROW_EXCEPTION(
				C_ptr.get() == NULL, std::logic_error
				,"DecompositionSystemOrthogonal::update_decomp(...) : Error, "
				"The matrix class used for the basis matrix C of type \'"
				<< typeid(*C).name() << "\' must return return.get() != NULL from the clone_mwons() method "
				"since mat_rel == MATRICES_INDEP_IMPS!" );
		}
		D_ptr_t  D_ptr = D;
		if(mat_rel == MATRICES_INDEP_IMPS) {
			D_ptr = D->clone();
			THROW_EXCEPTION(
				D_ptr.get() == NULL, std::logic_error
				,"DecompositionSystemOrthogonal::update_decomp(...) : Error, "
				"The matrix class used for the direct sensitivity matrix D = inv(C)*N of type \'"
				<< typeid(*D).name() << "\' must return return.get() != NULL from the clone() method "
				"since mat_rel == MATRICES_INDEP_IMPS!" );
		}
		if(S_ptr_.get() == NULL) {
			S_ptr_ = var_reduct_orthog_strategy().factory_S()->create();
		}
		var_reduct_orthog_strategy().update_S(*D_ptr,S_ptr_.get());
		R_orth->initialize(C_ptr,D_ptr,S_ptr_);
	}
	assert(Uy_cpst == NULL); // ToDo: Implement for undecomposed equalities
	assert(Vy_cpst == NULL); // ToDo: Implement for general inequalities

}

void DecompositionSystemOrthogonal::print_update_matrices(
	std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Orthogonal decompositon Y, R, Uy and Vy matrices\n"
		<< L << "Y  = [ I; -D' ] (using class MatrixIdentConcatStd)\n"
		<< L << "R  = C*(I + D*D')\n"
		<< L << "Uy = E - F*D'\n"
		<< L << "Vy = Gh(var_dep,:)' - Gh(var_indep,:)'*D'\n"
		;
}
		
}	// end namespace ConstrainedOptimizationPack
