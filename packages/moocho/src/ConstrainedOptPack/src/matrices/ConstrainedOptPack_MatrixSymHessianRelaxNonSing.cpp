// ////////////////////////////////////////////////////////////////
// MatrixSymHessianRelaxNonSing.cpp

#include <assert.h>

#include "ConstrainedOptimizationPack/include/MatrixSymHessianRelaxNonSing.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "Misc/include/profile_hack.h"

namespace ConstrainedOptimizationPack {

MatrixSymHessianRelaxNonSing::MatrixSymHessianRelaxNonSing(
	const G_ptr_t       &G_ptr
	,const VectorSlice  &M_diag
	)
{
	initialize(G_ptr,M_diag);
}

void MatrixSymHessianRelaxNonSing::initialize(
	const G_ptr_t       &G_ptr
	,const VectorSlice  &M_diag
	)
{
	if( G_ptr.get() != NULL ) {
		if( G_ptr->rows() == 0 )
			throw std::invalid_argument(
				"MatrixSymHessianRelaxNonSing::initialize(...): Error, if G_ptr.get() != NULL "
				"then G_ptr->rows() > 0 must be true" );
		if( M_diag.size() == 0 )
			throw std::invalid_argument(
				"MatrixSymHessianRelaxNonSing::initialize(...): Error, if G_ptr.get() != NULL "
				"then M_diag.size() > 0 must be true" );
		G_ptr_ = G_ptr;
		M_.init_diagonal(M_diag);
	}
	else {
		if( M_diag.size() > 0 )
			throw std::invalid_argument(
				"MatrixSymHessianRelaxNonSing::initialize(...): Error, if G_ptr.get() == NULL "
				"then M_diag.size() == 0 must be true" );
		G_ptr_ = NULL;
		M_.init_identity(0,0.0);
	}
}
	
const MatrixSymHessianRelaxNonSing::G_ptr_t& MatrixSymHessianRelaxNonSing::G_ptr() const
{
	return G_ptr_;
}

const MatrixSymWithOpFactorized& MatrixSymHessianRelaxNonSing::G() const
{
	return *G_ptr_;
}

const SparseLinAlgPack::MatrixSymDiagonalStd& MatrixSymHessianRelaxNonSing::M() const
{
	return M_;
}
	
// Overridden from Matrix

size_type MatrixSymHessianRelaxNonSing::rows() const
{
	return G_ptr_.get() ? G_ptr_->rows() + M_.rows() : 0;
}

// Overridden from MatrixWithOp

void MatrixSymHessianRelaxNonSing::Mp_StM(
	GenMatrixSlice* C, value_type a, BLAS_Cpp::Transp H_trans
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Mp_StM(...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	SparseLinAlgPack::Mp_StM( &(*C)(1,nG,1,nG), a, *G_ptr_, H_trans);
	SparseLinAlgPack::Mp_StM( &(*C)(nG+1,nG+nM,nG+1,nG+nM), a, M_, H_trans);
}

void MatrixSymHessianRelaxNonSing::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp H_trans
	, const VectorSlice& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StMtV(...VectorSlice...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	LinAlgPack::Vt_S(y,b); // Takes care of b == 0.0 and y uninitialized
	SparseLinAlgPack::Vp_StMtV( &(*y)(1,nG), a, *G_ptr_, H_trans, x(1,nG) );
	SparseLinAlgPack::Vp_StMtV( &(*y)(nG+1,nG+nM), a, M_, H_trans, x(nG+1,nG+nM) );
}

void MatrixSymHessianRelaxNonSing::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp H_trans
	, const SpVectorSlice& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StMtV(...SpVectorSlice...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	LinAlgPack::Vt_S(y,b); // Takes care of b == 0.0 and y uninitialized
	SparseLinAlgPack::Vp_StMtV( &(*y)(1,nG), a, *G_ptr_, H_trans, x(1,nG) );
	SparseLinAlgPack::Vp_StMtV( &(*y)(nG+1,nG+nM), a, M_, H_trans, x(nG+1,nG+nM) );
}

void MatrixSymHessianRelaxNonSing::Vp_StPtMtV(
	VectorSlice* y, value_type a, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp H_trans, const VectorSlice& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StPtMtV(...VectorSlice...)" );
#endif
	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // called for profiling
	// ToDo: Implement when needed!  After profiling
}

void MatrixSymHessianRelaxNonSing::Vp_StPtMtV(
	VectorSlice* y, value_type a, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp H_trans, const SpVectorSlice& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StPtMtV(...SpVectorSlice...)" );
#endif
	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // called for profiling
	// ToDo: Implement when needed!  After profiling
}

// Overridden from MatrixWithOpFactorized

void MatrixSymHessianRelaxNonSing::V_InvMtV(
	VectorSlice* y, BLAS_Cpp::Transp H_trans, const VectorSlice& x
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::V_InvMtV(...VectorSlice...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	SparseLinAlgPack::V_InvMtV( &(*y)(1,nG), *G_ptr_, H_trans, x(1,nG) );
	SparseLinAlgPack::V_InvMtV( &(*y)(nG+1,nG+nM), M_, H_trans, x(nG+1,nG+nM) );
}

void MatrixSymHessianRelaxNonSing::V_InvMtV(
	VectorSlice* y, BLAS_Cpp::Transp H_trans, const SpVectorSlice& x
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::V_InvMtV(...SpVectorSlice...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	SparseLinAlgPack::V_InvMtV( &(*y)(1,nG), *G_ptr_, H_trans, x(1,nG) );
	SparseLinAlgPack::V_InvMtV( &(*y)(nG+1,nG+nM), M_, H_trans, x(nG+1,nG+nM) );
}

// private

void MatrixSymHessianRelaxNonSing::assert_initialized() const
{
	if( G_ptr_.get() == NULL )
		throw std::logic_error(
			"MatrixSymHessianRelaxNonSing::assert_initialized(): Error, Not initalized yet!" );
}

} // end namespace ConstrainedOptimizationPack
