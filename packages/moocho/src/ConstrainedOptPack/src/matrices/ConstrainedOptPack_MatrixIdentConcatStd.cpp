// ////////////////////////////////////////////////////////////////////////
// MatrixIdentConcatStd.cpp

#include "ConstrainedOptimizationPack/include/MatrixIdentConcatStd.h"

namespace ConstrainedOptimizationPack {

// Setup and representation access

void MatrixIdentConcatStd::initialize(
	ETopBottom top_or_bottom, value_type alpha, const D_ptr_t& D_ptr, BLAS_Cpp::Transp D_trans)
{
	// validate input
	if( !D_ptr.get() )
		throw std::invalid_argument(
			"MatrixIdentConcatStd::initialize(...): Error, "
			"D_ptr.get() can not be NULL!" );
	const size_type
		D_rows   = D_ptr->rows(),
		D_cols   = D_ptr->cols(),
		opD_rows = BLAS_Cpp::rows( D_rows, D_cols, D_trans ),
		opD_cols = BLAS_Cpp::cols( D_rows, D_cols, D_trans ),
		rows     = opD_rows + opD_cols;
	alpha_    = alpha;
	D_ptr_    = D_ptr;
	D_trans_  = D_trans;
	D_rng_    = top_or_bottom == TOP ? Range1D(1,opD_rows)      : Range1D(opD_cols+1,rows);
	I_rng_    = top_or_bottom == TOP ? Range1D(opD_rows+1,rows) : Range1D(1,opD_cols);
}

void MatrixIdentConcatStd::set_uninitialized()
{
	alpha_    = 0.0;
	D_ptr_    = NULL;
	D_trans_  = BLAS_Cpp::no_trans;
	D_rng_    = Range1D();
	I_rng_    = Range1D();
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

const MatrixWithOp& MatrixIdentConcatStd::D() const
{
	return *D_ptr_;
}

BLAS_Cpp::Transp MatrixIdentConcatStd::D_trans() const
{
	return D_trans_;
}

// Overridden from MatrixWithOp

MatrixWithOp& MatrixIdentConcatStd::operator=(const MatrixWithOp& m)
{
	assert(0); // Finish!
	return *this;
}

} // end namespace ConstrainedOptimizationPack
