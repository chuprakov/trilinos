// /////////////////////////////////////////////////////////////
// DirectSparseSolverMA28.cpp
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

#ifdef SPARSE_SOLVER_PACK_USE_MA28

#include <assert.h>

#include <ostream>

#include "SparseSolverPack/include/DirectSparseSolverMA28.h"
#include "AbstractFactoryStd.h"
#include "ThrowException.h"

namespace SparseSolverPack {

// //////////////////////////////////////////////////
// DirectSparseSolverMA28::BasisMatrixMA28

// Overridden from BasisMatrixImp

MemMngPack::ref_count_ptr<DirectSparseSolverImp::BasisMatrixImp>
DirectSparseSolverMA28::BasisMatrixMA28::create_matrix() const
{
	return MemMngPack::rcp(new BasisMatrixMA28);
}

void DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(
	VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
	,const VectorWithOp& v_rhs2) const 
{
	assert(0); // ToDo: Implement!
}

// //////////////////////////////////////////////////
// DirectSparseSolverMA28

// Constructors/initializers

DirectSparseSolverMA28::DirectSparseSolverMA28()
	: estimated_fillin_ratio_(10.0)
{}

// Overridden from DirectSparseSolver

const DirectSparseSolver::basis_matrix_factory_ptr_t
DirectSparseSolverMA28::basis_matrix_factory() const
{
	namespace mmp = MemMngPack;
	return mmp::rcp(new mmp::AbstractFactoryStd<BasisMatrix,BasisMatrixMA28>());
}

void DirectSparseSolverMA28::estimated_fillin_ratio(
	value_type estimated_fillin_ratio
	)
{
	estimated_fillin_ratio_ = estimated_fillin_ratio;
}

// Overridden from DirectSparseSolverImp

const MemMngPack::ref_count_ptr<DirectSparseSolver::FactorizationStructure>
DirectSparseSolverMA28::create_fact_struc() const
{
	return MemMngPack::rcp(new FactorizationStructureMA28);
}

const MemMngPack::ref_count_ptr<DirectSparseSolverImp::FactorizationNonzeros>
DirectSparseSolverMA28::create_fact_nonzeros() const
{
	return MemMngPack::rcp(new FactorizationNonzerosMA28);
}

void DirectSparseSolverMA28::imp_analyze_and_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
	,LinAlgPack::IVector                            *row_perm
	,LinAlgPack::IVector                            *col_perm
	,size_type                                      *rank
	,std::ostream                                   *out
	)
{
	if(out)
		*out << "\nUsing MA28 to analyze and factor a new matrix ...\n";
	THROW_EXCEPTION(
		true, std::logic_error
		,"DirectSparseSolverMA28::imp_analyze_and_factor(...) : Error, not implemented yet!" );
}

void DirectSparseSolverMA28::imp_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
	,std::ostream                                   *out
	)
{
	assert(0); // ToDo: Implement!
}

}	// end namespace SparseSolverPack 

#endif // SPARSE_SOLVER_PACK_USE_MA28
