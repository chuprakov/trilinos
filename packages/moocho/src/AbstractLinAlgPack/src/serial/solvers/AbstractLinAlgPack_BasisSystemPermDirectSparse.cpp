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
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingularAggr.h"
#include "SparseLinAlgPack/include/MultiVectorMutableDense.h"
#include "AbstractFactoryStd.h"

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
	) const
{
	assert(0); // ToDo: Implement!
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
	) const
{
	assert(0); // ToDo: Implement!
}
	
} // end namespace SparseSolverPack
