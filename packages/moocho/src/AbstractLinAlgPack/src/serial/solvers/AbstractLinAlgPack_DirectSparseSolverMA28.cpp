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

#include <assert.h>

#include "SparseSolverPack/include/DirectSparseSolverMA28.h"

namespace SparseSolverPack {

// Overridden from DirectSparseSolver

const DirectSparseSolver::basis_matrix_factory_ptr_t
DirectSparseSolverMA28::basis_matrix_factory() const
{
	assert(0); // ToDo: Implement!
	return MemMngPack::null;
}

void DirectSparseSolverMA28::estimated_fillin_ratio(
	value_type estimated_fillin_ratio
	)
{
	assert(0); // ToDo: Implement!
}

// Overridden from DirectSparseSolverImp

const MemMngPack::ref_count_ptr<DirectSparseSolver::FactorizationStructure>
DirectSparseSolverMA28::create_fact_struc() const
{
	assert(0); // ToDo: Implement!
	return MemMngPack::null;
}

const MemMngPack::ref_count_ptr<DirectSparseSolverImp::FactorizationNonzeros>
DirectSparseSolverMA28::create_fact_nonzeros() const
{
	assert(0); // ToDo: Implement!
	return MemMngPack::null;
}

void DirectSparseSolverMA28::imp_analyze_and_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
	,LinAlgPack::IVector                            *row_perm
	,LinAlgPack::IVector                            *col_perm
	,size_type                                      *rank
	,BasisMatrixImp                                 *basis_matrix
	,std::ostream                                   *out            = NULL
	)
{
	assert(0); // ToDo: Implement!
}

void DirectSparseSolverMA28::imp_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
	,BasisMatrixImp                                 *basis_matrix
	,std::ostream                                   *out            = NULL
	)
{
	assert(0); // ToDo: Implement!
}

}	// end namespace SparseSolverPack 
