// ////////////////////////////////////////////////////////////////////////
// MatrixIdentConcatStd.cpp
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

#include "ConstrainedOptPack_MatrixIdentConcatStd.hpp"
#include "Teuchos_TestForException.hpp"

namespace ConstrainedOptPack {

// Setup and representation access

MatrixIdentConcatStd::MatrixIdentConcatStd()
{
	this->set_uninitialized();
}

void MatrixIdentConcatStd::initialize(
		const VectorSpace::space_ptr_t&    space_cols
		,const VectorSpace::space_ptr_t&   space_rows
		,ETopBottom                        top_or_bottom
		,value_type                        alpha
		,const D_ptr_t                     &D_ptr
		,BLAS_Cpp::Transp                  D_trans
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		space_cols.get() == NULL, std::invalid_argument
		,"MatrixIdentConcatStd::initialize(...): Error, "
		"space_cols.get() can not be NULL!" );
	TEST_FOR_EXCEPTION(
		space_rows.get() == NULL, std::invalid_argument
		,"MatrixIdentConcatStd::initialize(...): Error, "
		"space_rows.get() can not be NULL!" );
	TEST_FOR_EXCEPTION(
		D_ptr.get() == NULL, std::invalid_argument
		,"MatrixIdentConcatStd::initialize(...): Error, "
		"D_ptr.get() can not be NULL!" );
#endif
	const size_type
		D_rows   = D_ptr->rows(),
		D_cols   = D_ptr->cols(),
		opD_rows = BLAS_Cpp::rows( D_rows, D_cols, D_trans ),
		opD_cols = BLAS_Cpp::cols( D_rows, D_cols, D_trans ),
		rows     = opD_rows + opD_cols;
	space_cols_ = space_cols;
	space_rows_ = space_rows;
	alpha_      = alpha;
	D_ptr_      = D_ptr;
	D_trans_    = D_trans;
	D_rng_      = top_or_bottom == TOP ? Range1D(1,opD_rows)      : Range1D(opD_cols+1,rows);
	I_rng_      = top_or_bottom == TOP ? Range1D(opD_rows+1,rows) : Range1D(1,opD_cols);
}

void MatrixIdentConcatStd::set_uninitialized()
{
	namespace rcp = MemMngPack;
	space_cols_ = Teuchos::null;
	space_rows_ = Teuchos::null;
	alpha_      = 0.0;
	D_ptr_      = Teuchos::null;
	D_trans_    = BLAS_Cpp::no_trans;
	D_rng_      = Range1D::Invalid;
	I_rng_      = Range1D::Invalid;
}

const MatrixIdentConcatStd::D_ptr_t& MatrixIdentConcatStd::D_ptr() const
{
	return D_ptr_;
}

// Overridden form MatrixIdentConcat

Range1D MatrixIdentConcatStd::D_rng() const
{
	return D_rng_;
}

Range1D MatrixIdentConcatStd::I_rng() const
{
	return I_rng_;
}

value_type MatrixIdentConcatStd::alpha() const
{
	return alpha_;
}

const MatrixOp& MatrixIdentConcatStd::D() const
{
	return *D_ptr_;
}

BLAS_Cpp::Transp MatrixIdentConcatStd::D_trans() const
{
	return D_trans_;
}

// Overridden from MatrixOp

const VectorSpace& MatrixIdentConcatStd::space_cols() const
{
	return *space_cols_;
}

const VectorSpace& MatrixIdentConcatStd::space_rows() const
{
	return *space_rows_;
}

MatrixOp& MatrixIdentConcatStd::operator=(const MatrixOp& m)
{
	assert(0); // Finish!
	return *this;
}

// private

void MatrixIdentConcatStd::assert_initialized() const {
	TEST_FOR_EXCEPTION(
		space_cols_.get() == NULL, std::logic_error
		,"Error, the MatrixIdentConcatStd object has not been initialized!" );
}

} // end namespace ConstrainedOptPack
