// /////////////////////////////////////////////////////////////
// DirectSparseSolverImp.cpp
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

#include "SparseSolverPack/include/DirectSparseSolverImp.h"
#include "AbstractFactoryStd.h"
#include "ThrowException.h"
#include "dynamic_cast_verbose.h"

namespace SparseSolverPack {

// /////////////////////////////////////////
// DirectSparseSolverImp::BasisMatrixImp

DirectSparseSolverImp::BasisMatrixImp::BasisMatrixImp(
	size_type                      dim
	,const fact_struc_ptr_t        &fact_struc
	,const fact_nonzeros_ptr_t     &fact_nonzeros
	)
{
	this->initialize(dim,fact_struc,fact_nonzeros);
}

void DirectSparseSolverImp::BasisMatrixImp::initialize(
	size_type                      dim
	,const fact_struc_ptr_t        &fact_struc
	,const fact_nonzeros_ptr_t     &fact_nonzeros
	)
{
#ifdef _DEBUG
	const char msg_err[] = "DirectSparseSolverImp::BasisMatrixImp::initialize(...): Error!";
	THROW_EXCEPTION( dim < 0, std::logic_error, msg_err );
	THROW_EXCEPTION( fact_struc.get() == NULL, std::logic_error, msg_err );
	THROW_EXCEPTION( fact_nonzeros.get() == NULL, std::logic_error, msg_err );
#endif
	dim_            = dim;
	fact_struc_     = fact_struc;
	fact_nonzeros_  = fact_nonzeros;
}

void DirectSparseSolverImp::BasisMatrixImp::set_uninitialized()
{
	dim_            = 0;
	fact_struc_     = MemMngPack::null;
	fact_nonzeros_  = MemMngPack::null;
}

const DirectSparseSolverImp::BasisMatrixImp::fact_nonzeros_ptr_t&
DirectSparseSolverImp::BasisMatrixImp::get_fact_nonzeros() const
{
	return fact_nonzeros_;
}

// Overridden from MatrixBase

size_type DirectSparseSolverImp::BasisMatrixImp::rows() const
{
	return dim_;
}

size_type DirectSparseSolverImp::BasisMatrixImp::cols() const
{
	return dim_;
}

// Overridden from BasisMatrix

const DirectSparseSolver::BasisMatrix::fact_struc_ptr_t&
DirectSparseSolverImp::BasisMatrixImp::get_fact_struc() const
{
	return fact_struc_;
}

// //////////////////////////
// DirectSparseSolverImp

// Overridden from DirectSparseSolver

void DirectSparseSolverImp::analyze_and_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,LinAlgPack::IVector                            *row_perm
	,LinAlgPack::IVector                            *col_perm
	,size_type                                      *rank
	,BasisMatrix                                    *basis_matrix
	,std::ostream                                   *out
	)
{
	using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	const char msg_err[] = "DirectSparseSolverImp::analyze_and_factor(...): Error!";
	THROW_EXCEPTION( row_perm == NULL, std::logic_error, msg_err );
	THROW_EXCEPTION( col_perm == NULL, std::logic_error, msg_err );
	THROW_EXCEPTION( rank == NULL, std::logic_error, msg_err );
	THROW_EXCEPTION( basis_matrix == NULL, std::logic_error, msg_err );
#endif
	BasisMatrixImp
		&basis_matrix_imp = dyn_cast<BasisMatrixImp>(*basis_matrix);
	// Get references to factorization structure and factorizaton nonzeros
	// objects that are not being shared by other client.
	BasisMatrix::fact_struc_ptr_t         fact_struc;
	BasisMatrixImp::fact_nonzeros_ptr_t   fact_nonzeros;
	assert(0); // ToDo: Implement!
	// Now ask the subclass to do the work
	this->imp_analyze_and_factor(
		A,fact_struc,fact_nonzeros,row_perm,col_perm,rank,&basis_matrix_imp,out
		);
}

void DirectSparseSolverImp::factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,BasisMatrix                                    *basis_matrix
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,std::ostream                                   *out
	)
{
	assert(0); // ToDo: Implement!
}

const DirectSparseSolver::BasisMatrix::fact_struc_ptr_t&
DirectSparseSolverImp::get_fact_struc() const
{
	return fact_struc_;
}

void DirectSparseSolverImp::set_uninitialized()
{
	fact_struc_ = MemMngPack::null;
}

}	// end namespace SparseSolverPack 
