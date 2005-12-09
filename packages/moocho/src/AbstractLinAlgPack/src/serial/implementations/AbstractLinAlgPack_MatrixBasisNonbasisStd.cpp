// //////////////////////////////////////////////////////////////////////////////////
// MatrixBasisNonbasisStd.cpp
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

#include <stdexcept>
#include <string>

#include "AbstractLinAlgPack_MatrixBasisNonbasisStd.hpp"

namespace AbstractLinAlgPack {

MatrixBasisNonbasisStd::MatrixBasisNonbasisStd()
	: C_(NULL), C_nonsingular_(NULL), N_(NULL)
{}

void MatrixBasisNonbasisStd::initialize(
	  const C_ptr_t&                  C
	, const C_nonsingular_ptr_t&      C_nonsingular
	, const N_ptr_t&                  N
	)
{
	if( C.get() == NULL && N.get() != NULL )
		throw std::logic_error(
			"MatrixBasisNonbasisStd::initialize(...) : Error , "
			"If N!=NULL then C must be non-NULL" );
	if( C.get() == NULL && C_nonsingular.get() != NULL )
		throw std::logic_error(
			"MatrixBasisNonbasisStd::initialize(...) : Error , "
			"If C_nonsingular!=NULL then C must be non-NULL" );
	if( C.get() != NULL ) {
		if( C->rows() != C->cols() )
			throw std::length_error(
				"MatrixBasisNonbasisStd::initialize(...) : Error , "
				"C->rows() must equal C->cols()" );
	 	if( N.get()!=NULL && C->rows() != N->rows() ) {
			throw std::length_error(
				"MatrixBasisNonbasisStd::initialize(...) : Error , "
				"N->rows() must equal C->cols()" );
		}
	 	if( C_nonsingular.get()!=NULL
			&& (C_nonsingular->rows() != C->rows() || C_nonsingular->cols() != C->cols() ) )
		{
			throw std::length_error(
				"MatrixBasisNonbasisStd::initialize(...) : Error , "
				"C_nonsingular must be the same dimension as C" );
		}
	}
	C_             = C;             // Will delete the Old C if C_->count() == 1
	C_nonsingular_ = C_nonsingular; // Same here
	N_             = N;             // Same here
} 

// Overridden from MatrixOp

size_type MatrixBasisNonbasisStd::rows() const
{
	return C_.get() ? C_->cols() + (N_.get() ? N_->cols() : 0) : 0;
}

size_type MatrixBasisNonbasisStd::cols() const
{
	return C_.get() ? C_->rows() : 0;
}
	
// Overridden from MatrixBasisNonbasis

const MatrixOp& MatrixBasisNonbasisStd::C() const
{
	if(C_.get()==NULL)
		throw std::logic_error(
			"MatrixBasisNonbasisStd::C() : Error, "
			"Matrix not initialized!" );
	return *C_;
}

const MatrixFactorized& MatrixBasisNonbasisStd::C_nonsingular() const
{
	if(C_.get()==NULL)
		throw std::logic_error(
			"MatrixBasisNonbasisStd::C_nonsingular() : Error, "
			"Matrix not set!" );
	return *C_nonsingular_;
}

const MatrixOp& MatrixBasisNonbasisStd::N() const
{
	if(N_.get()==NULL)
		throw std::logic_error(
			"MatrixBasisNonbasisStd::N() : Error, "
			"Matrix has not been initialized or there is no N matrix!" );
	return *N_;
}

}	// end namespace AbstractLinAlgPack 
