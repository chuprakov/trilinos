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

void DirectSparseSolverMA28::analyze_and_factor(
	const MatrixConvertToSparse     &A
	,IVector                        *row_perm
	,IVector                        *col_perm
	,size_type                      *rank
	,BasisMatrix                    *basis_matrix
	,std::ostream                   *out
	)
{
	assert(0); // ToDo: Implement!
}

void DirectSparseSolverMA28::factor(
	const MatrixConvertToSparse              &A
	,const BasisMatrix::fact_struc_ptr_t     &fact_struc
	,BasisMatrix                             *basis_matrix
	,std::ostream                            *out
	)
{
	assert(0); // ToDo: Implement!
}

const DirectSparseSolverMA28::basis_matrix_ptr_t&
DirectSparseSolverMA28::get_basis_matrix() const
{
	assert(0); // ToDo: Implement!
	return MemMngPack::null;
}

void DirectSparseSolverMA28::release_memory()
{
	assert(0); // ToDo: Implement!
}

}	// end namespace SparseSolverPack 
