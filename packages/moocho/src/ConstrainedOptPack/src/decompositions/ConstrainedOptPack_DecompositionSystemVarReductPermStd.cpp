// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystemVarReductPermStd.cpp
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

#include "ConstrainedOptimizationPack/include/DecompositionSystemVarReductPermStd.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystemVarReductImp.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/BasisSystemPerm.h"
#include "ThrowException.h"

namespace ConstrainedOptimizationPack {

// Constructors / initializers

DecompositionSystemVarReductPermStd::DecompositionSystemVarReductPermStd(
	const decomp_sys_imp_ptr_t&        decomp_sys_imp
	,const basis_sys_ptr_t&            basis_sys
	,bool                              basis_selected
	)
{
	this->initialize(decomp_sys_imp,basis_sys,basis_selected);
}

void DecompositionSystemVarReductPermStd::initialize(
	const decomp_sys_imp_ptr_t&        decomp_sys_imp
	,const basis_sys_ptr_t&            basis_sys
	,bool                              basis_selected
	)
{
	decomp_sys_imp_ = decomp_sys_imp;
	basis_sys_      = basis_sys;
	basis_selected_ = basis_selected;
}

// Overridden from DecompositionSystem

size_type DecompositionSystemVarReductPermStd::n() const
{
	return decomp_sys_imp()->n();
}

size_type DecompositionSystemVarReductPermStd::m() const
{
	return decomp_sys_imp()->m();
}

size_type DecompositionSystemVarReductPermStd::r() const
{
	return decomp_sys_imp()->r();
}

Range1D DecompositionSystemVarReductPermStd::con_decomp() const
{
	return decomp_sys_imp()->con_decomp();
}

Range1D DecompositionSystemVarReductPermStd::con_undecomp() const
{
	return decomp_sys_imp()->con_undecomp();
}

const VectorSpace::space_ptr_t
DecompositionSystemVarReductPermStd::space_range() const
{
	return decomp_sys_imp()->space_range();
}

const VectorSpace::space_ptr_t
DecompositionSystemVarReductPermStd::space_null() const
{
	return decomp_sys_imp()->space_null();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Z() const
{
	return decomp_sys_imp()->factory_Z();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Y() const
{
	return decomp_sys_imp()->factory_Y();
}

const DecompositionSystem::mat_nonsing_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_R() const
{
	return decomp_sys_imp()->factory_R();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Uz() const
{
	return decomp_sys_imp()->factory_Uz();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Uy() const
{
	return decomp_sys_imp()->factory_Uy();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Vz() const
{
	return decomp_sys_imp()->factory_Vz();
}

const DecompositionSystem::mat_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_Vy() const
{
	return decomp_sys_imp()->factory_Vy();
}

void DecompositionSystemVarReductPermStd::update_decomp(
	std::ostream              *out
	,EOutputLevel             olevel
	,ERunTests                test_what
	,const MatrixWithOp       &Gc
	,const MatrixWithOp       *Gh
	,MatrixWithOp             *Z
	,MatrixWithOp             *Y
	,MatrixWithOpNonsingular  *R
	,MatrixWithOp             *Uz
	,MatrixWithOp             *Uy
	,MatrixWithOp             *Vz
	,MatrixWithOp             *Vy
	,EMatRelations            mat_rel
	) const
{
	assert_basis_selected();
	decomp_sys_imp()->update_decomp(
		out,olevel,test_what,Gc,Gh,Z,Y
		,R,Uz,Uy,Vz,Vy,mat_rel
		);
}

void DecompositionSystemVarReductPermStd::print_update_decomp(
	std::ostream& out, const std::string& L ) const
{
	// ToDo: Print basis permutation stuff also?
	decomp_sys_imp()->print_update_decomp(out,L);
}

// @name Overridden from DecompositionSystemVarReductPerm

const DecompositionSystemVarReductPerm::perm_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_P_var() const
{
	return basis_sys()->factory_P_var();
}

const DecompositionSystemVarReductPerm::perm_fcty_ptr_t
DecompositionSystemVarReductPermStd::factory_P_equ() const
{
	return basis_sys()->factory_P_equ();
}

bool DecompositionSystemVarReductPermStd::has_basis() const
{
	return basis_selected_;
}

void DecompositionSystemVarReductPermStd::set_decomp(
	std::ostream              *out
	,EOutputLevel             olevel
	,ERunTests                test_what
	,const Permutation        &P_var
	,const Range1D            &var_dep
	,const Permutation        *P_equ
	,const Range1D            *equ_decomp
	,const MatrixWithOp       &Gc
	,const MatrixWithOp       *Gh
	,MatrixWithOp             *Z
	,MatrixWithOp             *Y
	,MatrixWithOpNonsingular  *R
	,MatrixWithOp             *Uz
	,MatrixWithOp             *Uy
	,MatrixWithOp             *Vz
	,MatrixWithOp             *Vy
	,EMatRelations            mat_rel
	)
{
	// Forward these setting on to the implementation.
	decomp_sys_imp_->D_imp(  this->D_imp()  );
	decomp_sys_imp_->Uz_imp( this->Uz_imp() );
	decomp_sys_imp_->Vz_imp( this->Vz_imp() );
	// Get smart pointers to the basis matrix and the direct sensistivity matrices
	// and remove references to these matrix objects from the other decomposition
	// matrices by uninitializing them.
	MemMngPack::ref_count_ptr<MatrixWithOpNonsingular>  C_ptr;
	MemMngPack::ref_count_ptr<MatrixWithOp>             D_ptr;
	if( decomp_sys_imp_->basis_sys().get() ) {
		// It is assumed that the decomposition may have already
		// been updated so try to recycle storage.
		decomp_sys_imp_->get_basis_matrices(
			out, olevel, test_what
			,Z, Y, R, Uz, Uy, Vz, Vy
			,&C_ptr
			,this->D_imp() == MAT_IMP_EXPLICIT ? &D_ptr : NULL
			);
	}
	else {
		// It is assumed that the decomposition has not been
		// previously updated so we must allocate new storage.
		C_ptr = basis_sys_->factory_C()->create();
		if(this->D_imp() == MAT_IMP_EXPLICIT)
			D_ptr = basis_sys_->factory_D()->create();
	}
	// Tell the basis system object to set this basis
	try {
		basis_sys_->set_basis(
			P_var, var_dep
			,P_equ, equ_decomp
			,NULL, NULL
			,&Gc, Gh
			,C_ptr.get()
			,D_ptr.get()
			,this->Uz_imp() == MAT_IMP_EXPLICIT ? Uz : NULL
			,this->Vz_imp() == MAT_IMP_EXPLICIT ? Vz : NULL
			,(mat_rel == MATRICES_INDEP_IMPS
			  ? BasisSystem::MATRICES_INDEP_IMPS : BasisSystem::MATRICES_ALLOW_DEP_IMPS )
			,out
			);
	}
	catch( const BasisSystem::SingularBasis& except ) {
		if(out && olevel >= PRINT_BASIC_INFO)
			*out << "Passed in basis is singular, throwing SingularDecomposition: "
				 << except.what() << std::endl;
		THROW_EXCEPTION(
			true, SingularDecomposition
			,"DecompositionSystemVarReductPermStd::set_decomp(...): Passed in basis selection "
			"gave a singular basis matrix! : " << except.what() );
	}
	// If we get here the passed in basis selection is nonsingular and the basis matrices
	// are updated.  If the basis_sys object has not been set to the implementation then
	// do it.
	if( decomp_sys_imp_->basis_sys().get() == NULL ) {
		decomp_sys_imp_->initialize(
			decomp_sys_imp_->space_x()
			,decomp_sys_imp_->space_c()
			,decomp_sys_imp_->space_h()
			,basis_sys_
			);
	}
	// Now give them back to the decomp_sys_imp object and update the rest
	// of the decomposition matrices.
	decomp_sys_imp_->set_basis_matrices(
		out, olevel, test_what
		,C_ptr
		,this->D_imp() == MAT_IMP_EXPLICIT ? &D_ptr : NULL
		,this->Uz_imp() == MAT_IMP_EXPLICIT ? Uz : NULL
		,this->Vz_imp() == MAT_IMP_EXPLICIT ? Vz : NULL
		,decomp_sys_imp_->basis_sys().get() == NULL ? basis_sys_ : MemMngPack::null
		);
	decomp_sys_imp()->update_decomp(
		out,olevel,test_what,Gc,Gh,Z,Y
		,R,Uz,Uy,Vz,Vy,mat_rel
		);
}

void DecompositionSystemVarReductPermStd::select_decomp(
	std::ostream              *out
	,EOutputLevel             olevel
	,ERunTests                test_what
	,const VectorWithOp       *nu
	,MatrixWithOp             *Gc
	,MatrixWithOp             *Gh
	,Permutation              *P_var
	,Range1D                  *var_dep
	,Permutation              *P_equ
	,Range1D                  *equ_decomp
	,MatrixWithOp             *Z
	,MatrixWithOp             *Y
	,MatrixWithOpNonsingular  *R
	,MatrixWithOp             *Uz
	,MatrixWithOp             *Uy
	,MatrixWithOp             *Vz
	,MatrixWithOp             *Vy
	,EMatRelations            mat_rel
	)
{
/*
	basis_sys_->select_basis(
		nu, NULL, Gc, Gh
		,P_var, var_dep
		,P_equ, equ_decomp
		,NULL, NULL
		,???
		);
*/
	assert(0); // Todo: Implement!
	decomp_sys_imp_->initialize(
		decomp_sys_imp_->space_x()
		,decomp_sys_imp_->space_c()
		,decomp_sys_imp_->space_h()
		,basis_sys_
		);
}

// private

void DecompositionSystemVarReductPermStd::assert_basis_selected() const
{
	THROW_EXCEPTION(
		!basis_selected_, std::logic_error
		,"DecompositionSystemVarReductPermStd::assert_basis_selected(): Error, "
		"the methods set_decomp() or select_decomp() must be called first!" );
}

}	// end namespace ConstrainedOptimizationPack
