// ////////////////////////////////////////////////////////////
// BasisSystemPermDirectSparse.cpp
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

#include "SparseSolverPack/include/BasisSystemPermDirectSparse.h"
#include "SparseLinAlgPack/include/PermutationSerial.h"
#include "SparseLinAlgPack/include/MatrixConvertToSparseEncap.h"
#include "SparseLinAlgPack/include/MatrixExtractSparseElements.h"
#include "LinAlgPack/include/IVector.h"
#include "AbstractLinAlgPack/include/MatrixWithOp.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingularAggr.h"
#include "AbstractLinAlgPack/include/MatrixPermAggr.h"
#include "SparseLinAlgPack/include/MultiVectorMutableDense.h"
#include "AbstractFactoryStd.h"
#include "ThrowException.h"
#include "dynamic_cast_verbose.h"

namespace SparseSolverPack {

BasisSystemPermDirectSparse::BasisSystemPermDirectSparse(
	const direct_solver_ptr_t&   direct_solver
	)
{
	this->initialize(direct_solver);
}
	
void BasisSystemPermDirectSparse::initialize(
	const direct_solver_ptr_t&   direct_solver
	)
{
	direct_solver_ = direct_solver;
	// ToDo: Reinitialze things if we need to?
}

// Overridden from BasisSystem

const BasisSystem::mat_nonsing_fcty_ptr_t
BasisSystemPermDirectSparse::factory_C() const
{
	return MemMngPack::rcp(
		new MemMngPack::AbstractFactoryStd<MatrixWithOpNonsingular,MatrixWithOpNonsingularAggr>()
		);
}

const BasisSystem::mat_fcty_ptr_t
BasisSystemPermDirectSparse::factory_D() const
{
	return MemMngPack::rcp(
		new MemMngPack::AbstractFactoryStd<MatrixWithOp,MultiVectorMutableDense>()
		);
}

const BasisSystem::mat_fcty_ptr_t
BasisSystemPermDirectSparse::factory_GcUP() const
{
	return MemMngPack::rcp(
		new MemMngPack::AbstractFactoryStd<MatrixWithOp,MultiVectorMutableDense>()
		);
}

const BasisSystem::mat_fcty_ptr_t
BasisSystemPermDirectSparse::factory_GhUP() const
{
	return MemMngPack::rcp(
		new MemMngPack::AbstractFactoryStd<MatrixWithOp,MultiVectorMutableDense>()
		);
}

Range1D BasisSystemPermDirectSparse::var_dep() const
{
	return Range1D(1,r_);
}

Range1D BasisSystemPermDirectSparse::var_indep() const
{
	return r_ < n_ ? Range1D(r_+1,n_) : Range1D::Invalid;
}

void BasisSystemPermDirectSparse::update_basis(
	const MatrixWithOp*         Gc
	,const MatrixWithOp*        Gh
	,MatrixWithOpNonsingular*   C
	,MatrixWithOp*              D
	,MatrixWithOp*              GcUP
	,MatrixWithOp*              GhUP
	,EMatRelations              mat_rel
	,std::ostream               *out
	) const
{
	assert(0); // ToDo: Implement!
}

// Overridded from BasisSystemPerm

const AbstractLinAlgPack::BasisSystemPerm::perm_fcty_ptr_t
BasisSystemPermDirectSparse::factory_P_var() const
{
	assert(0); // ToDo: Implement using PermutationSerial
	return MemMngPack::null;
}

const AbstractLinAlgPack::BasisSystemPerm::perm_fcty_ptr_t
BasisSystemPermDirectSparse::factory_P_equ() const
{
	assert(0); // ToDo: Implement using PermutationSerial
	return MemMngPack::null;
}

const AbstractLinAlgPack::BasisSystemPerm::perm_fcty_ptr_t
BasisSystemPermDirectSparse::factory_P_inequ() const
{
	assert(0); // ToDo: Implement using PermutationSerial
	return MemMngPack::null;
}

void BasisSystemPermDirectSparse::set_basis(
	const Permutation          &P_var
	,const Range1D             &var_dep
	,const Permutation         *P_equ
	,const Range1D             *equ_decomp
	,const Permutation         *P_inequ
	,const Range1D             *inequ_decomp
	,const MatrixWithOp        *Gc
	,const MatrixWithOp        *Gh
	,MatrixWithOpNonsingular   *C
	,MatrixWithOp              *D
	,MatrixWithOp              *GcUP
	,MatrixWithOp              *GhUP
	,EMatRelations              mat_rel
	,std::ostream               *out
	) const
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	if(out)
		*out << "\nUsing a direct sparse solver to set the basis ...\n";
#ifdef _DEBUG
	// Validate input
	THROW_EXCEPTION(
		Gc == NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error, "
		"Must have equality constriants in this current implementation! " );
	THROW_EXCEPTION(
		P_equ == NULL || equ_decomp == NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error!" );
	if(Gh) THROW_EXCEPTION(
		P_inequ != NULL && inequ_decomp != NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error, "
		"Can not handle decomposed inequalities yet!" );
	THROW_EXCEPTION(
		C == NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error!" );
	if(Gh) THROW_EXCEPTION(
		var_dep.size() != equ_decomp->size(), std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error!" );
#endif
	// Get the aggreate matrix object for Gc
	const MatrixPermAggr	
		&Gc_pa = dyn_cast<const MatrixPermAggr>(*Gc);
	// Get the basis matrix object from the aggregate or allocate one
	MatrixWithOpNonsingularAggr
		&C_aggr = dyn_cast<MatrixWithOpNonsingularAggr>(*C);
	mmp::ref_count_ptr<DirectSparseSolver::BasisMatrix>
		C_bm;
	if( C_aggr.mns().get() ) {
		C_bm = mmp::rcp_dynamic_cast<DirectSparseSolver::BasisMatrix>(
			mmp::rcp_const_cast<MatrixNonsingular>(C_aggr.mns() ) );
		if(C_bm.get() == NULL)
			dyn_cast<const DirectSparseSolver::BasisMatrix>(*C_aggr.mns()); // Throws exception!
	}
	else {
		C_bm = direct_solver_->basis_matrix_factory()->create();
	}
	// Get the concreate type of the direct sensitivity matrix (if one was passed in)
	MultiVectorMutableDense
		*D_mvd = NULL;
	if( D ) {
		D_mvd = &dyn_cast<MultiVectorMutableDense>(*D);
	}
	// Get at the concreate permutation vectors
	const PermutationSerial
		&P_var_s = dyn_cast<const PermutationSerial>(P_var),
		&P_equ_s = dyn_cast<const PermutationSerial>(*P_equ);
	// Setup the encapulated convert-to-sparse matrix object
	MatrixConvertToSparseEncap
		C_mctse(
			mmp::rcp_dynamic_cast<const MatrixExtractSparseElements>(Gc_pa.mat_orig())
			,P_var_s.inv_perm()
			,var_dep
			,P_equ_s.inv_perm()
			,*equ_decomp
			,BLAS_Cpp::trans
			);
	// Analyze and factor this basis (it had better be full rank)!
	IVector row_perm_ds, col_perm_ds; // Must store these even though we don't want them!
	size_type rank = 0;
	direct_solver_->analyze_and_factor(
		C_mctse
		,&row_perm_ds
		,&col_perm_ds
		,&rank
		,C_bm.get()
		,out
		);
	if( rank < var_dep.size() ) {
		assert(0); // ToDo: Throw an exception with a good error message!
	}
	// Compute the auxiliary projected matrices
	if( D_mvd ) {
		// D = -inv(C) * N
		Range1D
			var_indep = ( var_dep.lbound() == 1
						  ? Range1D(var_dep.size()+1,P_var.dim())
						  : Range1D(1,P_var.dim()-var_dep.size()) ); 
		D_mvd->initialize(var_dep.size(),var_indep.size());
		AbstractLinAlgPack::M_StInvMtM(
			D_mvd, -1.0, *C_bm, BLAS_Cpp::no_trans
			,*Gc->sub_view(var_indep,*equ_decomp),BLAS_Cpp::trans // N = Gc(var_indep,equ_decomp)'
			);
	}
	if( GcUP ) {
		assert(0); // ToDo: Implement!
	}
	if( GhUP ) {
		assert(0); // ToDo: Implement!
	}
	// Initialize the aggregate basis matrix object.
	C_aggr.initialize(
		Gc->sub_view(var_dep,*equ_decomp) // handled by the MatrixPermAggr subclass!
		,BLAS_Cpp::trans
		,C_bm
		,BLAS_Cpp::no_trans
		);
}

void BasisSystemPermDirectSparse::select_basis(
	const VectorWithOp          *nu
	,const VectorWithOp         *lambdaI
	,MatrixWithOp               *Gc
	,MatrixWithOp               *Gh
	,Permutation                *P_var
	,Range1D                    *var_dep
	,Permutation                *P_equ
	,Range1D                    *equ_decomp
	,Permutation                *P_inequ
	,Range1D                    *inequ_decomp
	,MatrixWithOpNonsingular    *C
	,MatrixWithOp               *D
	,MatrixWithOp               *GcUP
	,MatrixWithOp               *GhUP
	,EMatRelations              mat_rel
	,std::ostream               *out
	) const
{
	assert(0); // ToDo: Implement!
}
	
} // end namespace SparseSolverPack
