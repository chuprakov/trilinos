// ////////////////////////////////////////////////////////
// MatrixHessianSuperBasic.cpp

#include "ConstrainedOptimizationPack/include/MatrixHessianSuperBasic.h"
#include "SparseLinAlgPack/include/MatrixWithOpOut.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"

namespace ConstrainedOptimizationPack {

void MatrixHessianSuperBasic::initialize(
	size_type            n
	,size_type           n_R
	,const size_type     i_x_free[]
	,const size_type     i_x_fixed[]
	,const EBounds       bnd_fixed[]
	,const B_RR_ptr_t&   B_RR_ptr
	,const B_RX_ptr_t&   B_RX_ptr
	,BLAS_Cpp::Transp    B_RX_trans
	,const B_XX_ptr_t&   B_XX_ptr
	)
{
	using LinAlgPack::Mp_M_assert_sizes;
	using BLAS_Cpp::no_trans;

	const size_type
		n_X = n - n_R;

    // Validate input arguments

	// B_RR
	if(n_R > 0 ) {
		if( !B_RR_ptr.get() )
			throw std::invalid_argument(
				"MatrixHessianSuperBasic::initialize(...) : Error, "
				"B_RR_ptr.get() can not be NULL when n_R > 0" );
		Mp_M_assert_sizes( n_R, n_R, no_trans, B_RR_ptr->rows(), B_RR_ptr->cols(), no_trans );
	}
	else {
		if( B_RR_ptr.get() )
			throw std::invalid_argument(
				"MatrixHessianSuperBasic::initialize(...) : Error, "
				"B_RR_ptr.get() must be NULL when n_R == 0" );
	}
	// op(B_RX)
	if( n_R < n ) {
		if( B_RX_ptr.get() ) {
			Mp_M_assert_sizes( n_R, n_X, no_trans, B_RX_ptr->rows(), B_RX_ptr->cols(), B_RX_trans );
		}
	}
	else {
		if( B_RX_ptr.get() )
			throw std::invalid_argument(
				"MatrixHessianSuperBasic::initialize(...) : Error, "
				"B_RX_ptr.get() must be NULL when n_R == n" );
	}
	// B_XX
	if( n_R < n ) {
		if( !B_XX_ptr.get() )
			throw std::invalid_argument(
				"MatrixHessianSuperBasic::initialize(...) : Error, "
				"B_XX_ptr.get() can not be NULL if n_R < n" );
		Mp_M_assert_sizes( n_X, n_X, no_trans, B_XX_ptr->rows(), B_XX_ptr->cols(), no_trans );
	}
	else {
		if( B_XX_ptr.get() )
			throw std::invalid_argument(
				"MatrixHessianSuperBasic::initialize(...) : Error, "
				"B_XX_ptr.get() must be NULL when n_R == n" );
	}

	// Setup Q_R and Q_X and validate i_x_free[] and i_x_fixed[]
	assert(0); // ToDo: Implement!

	// Setup bnd_fixed
	bnd_fixed_.resize(n_X);
	{for(size_type i = 0; i < n_X; ++i) bnd_fixed_[i] = bnd_fixed[i]; }

	// Set the rest of the arguments
	B_RR_ptr_    = B_RR_ptr;
	B_RX_ptr_    = B_RX_ptr;
	B_RX_trans_  = B_RX_trans;
	B_XX_ptr_    = B_XX_ptr;

}

// Overridden from Matrix

size_type MatrixHessianSuperBasic::rows() const
{
	return (B_RR_ptr_.get()?B_RR_ptr_->rows():0)+(B_XX_ptr_.get()?B_XX_ptr_->rows():0);
}

// Overridden from MatrixWithOp

MatrixWithOp& MatrixHessianSuperBasic::operator=(const MatrixWithOp& m)
{
	assert_initialized();
	assert(0);
	return *this;
}

std::ostream& MatrixHessianSuperBasic::output(std::ostream& out) const
{
	assert_initialized();
	assert(0);
}

void MatrixHessianSuperBasic::Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2, value_type beta) const
{
	assert_initialized();
	assert(0);
}

void MatrixHessianSuperBasic::Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	assert_initialized();
	assert(0);
}

void MatrixHessianSuperBasic::Vp_StPtMtV(
	VectorSlice* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const VectorSlice& vs_rhs3, value_type beta) const
{
	assert_initialized();
	assert(0);
}

void MatrixHessianSuperBasic::Vp_StPtMtV(
	VectorSlice* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta) const
{
	assert_initialized();
	assert(0);
}

value_type MatrixHessianSuperBasic::transVtMtV(
	const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
	, const SpVectorSlice& sv_rhs3) const
{
	assert_initialized();
	assert(0);
	return 0.0;
}

// Private

void MatrixHessianSuperBasic::assert_initialized() const
{
	if( !B_RR_ptr_.get() && !B_XX_ptr_.get() )
		throw std::logic_error(
			"MatrixHessianSuperBasic::assert_initialized() : Error, "
			"The matrix is not initialized yet" );
}

} // end namespace ConstrainedOptimizationPack
