// ///////////////////////////////////////////////////////
// MatrixGenBanded.cpp

#include <assert.h>
#include <sstream>

#include "ConstrainedOptimizationPack/include/MatrixGenBanded.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/BLAS_Cpp.h"
#include "Misc/include/WorkspacePack.h"

namespace ConstrainedOptimizationPack {

MatrixGenBanded::MatrixGenBanded(
	size_type                         m
	,size_type                        n
	,size_type                        kl
	,size_type                        ku
	,GenMatrixSlice                   *MB
	,const release_resource_ptr_t&    MB_release_resource_ptr
	)
{
	initialize(m,n,kl,ku,MB,MB_release_resource_ptr);
}

void MatrixGenBanded::initialize(
	size_type                         m
	,size_type                        n
	,size_type                        kl
	,size_type                        ku
	,GenMatrixSlice                   *MB
	,const release_resource_ptr_t&    MB_release_resource_ptr
	)
{
	// Validate input

	if( m == 0 ) {
		if( n != 0 )
			throw std::invalid_argument(
				"MatrixGenBanded::initialize(...): Error, "
				"n must be 0 if m == 0" );
		if( kl != 0 )
			throw std::invalid_argument(
				"MatrixGenBanded::initialize(...): Error, "
				"kl must be 0 if m == 0" );
		if( ku != 0 )
			throw std::invalid_argument(
				"MatrixGenBanded::initialize(...): Error, "
				"ku must be 0 if m == 0" );
		if( MB != NULL )
			throw std::invalid_argument(
				"MatrixGenBanded::initialize(...): Error, "
				"MB must be NULL if m == 0" );
		if( MB_release_resource_ptr.get() != NULL )
			throw std::invalid_argument(
				"MatrixGenBanded::initialize(...): Error, "
				"MB_release_resource_ptr.get() must be NULL if m == 0" );
	}
	else {
		if( kl + 1 > m )
			throw std::invalid_argument(
				"MatrixGenBanded::initialize(...): Error, "
				"kl + 1 can not be larger than m" );
		if( ku + 1 > n )
			throw std::invalid_argument(
				"MatrixGenBanded::initialize(...): Error, "
				"ku + 1 can not be larger than n" );
		if( MB == NULL )
			throw std::invalid_argument(
				"MatrixGenBanded::initialize(...): Error, "
				"MB must not be NULL if n > 0" );
	}

	// Set the members

	if( m == 0 ) {
		m_                        = 0;
		n_                        = 0;
		kl_                       = 0;
		ku_                       = 0;
		MB_.bind(GenMatrixSlice());
		MB_release_resource_ptr_  = NULL;
	}
	else {
		// Set the members
		m_                        = m;
		n_                        = n;
		kl_                       = kl;
		ku_                       = ku;
		MB_.bind(*MB);
	}
}

// Overridden from MatrixWithOp

size_type MatrixGenBanded::rows() const
{
	return m_;
}

size_type MatrixGenBanded::cols() const
{
	return n_;
}

size_type MatrixGenBanded::nz() const
{
	return (ku_ + kl_ + 1) * n_ - ( (ku_+1) * (ku_+1) - (ku_+1) )/2  - ( (kl_+1) * (kl_+1) - (kl_+1) )/2; // Is correct?
}

std::ostream& MatrixGenBanded::output(std::ostream& out) const
{
	return MatrixWithOp::output(out); // ToDo: Implement specialized version later!
}

void MatrixGenBanded::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
	, const VectorSlice& x, value_type b) const
{
	assert_initialized();
	LinAlgPack::Vp_MtV_assert_sizes( y->size(), n_, n_, BLAS_Cpp::no_trans, x.size() );
	BLAS_Cpp::gbmv(M_trans,m_,n_,kl_,ku_,a,MB_.col_ptr(1),MB_.max_rows(),x.raw_ptr(),x.stride()
				   ,b,y->raw_ptr(),y->stride());
}

void MatrixGenBanded::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b) const
{
	assert_initialized();
	MatrixWithOp::Vp_StMtV(y,a,M_trans,x,b); // ToDo: Implement spacialized operation when needed!
}

void MatrixGenBanded::Vp_StPtMtV(
	VectorSlice* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp M_trans
	, const VectorSlice& x, value_type b) const
{
	assert_initialized();
	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,b); // ToDo: Implement spacialized operation when needed!
}

void MatrixGenBanded::Vp_StPtMtV(
	VectorSlice* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b) const
{
	assert_initialized();
	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,b); // ToDo: Implement spacialized operation when needed!
}

// Private member functions

void MatrixGenBanded::assert_initialized() const
{
	if( m_ == 0 )
		throw std::logic_error("MatrixGenBanded::assert_initialized(): Error, "
							   "not initialized!" );
}

} // end namespace ConstrainedOptimizationPack
