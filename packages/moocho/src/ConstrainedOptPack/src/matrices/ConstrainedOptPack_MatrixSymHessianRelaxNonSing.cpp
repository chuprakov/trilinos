// ////////////////////////////////////////////////////////////////
// MatrixSymHessianRelaxNonSing.cpp

#include <assert.h>

#include "ConstrainedOptimizationPack/include/MatrixSymHessianRelaxNonSing.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/GenPermMatrixSlice.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"
#include "Misc/include/profile_hack.h"

namespace {

//
template<class V>
void Vp_StPtMtV_imp( 
	LinAlgPack::VectorSlice* y, LinAlgPack::value_type a
	, const SparseLinAlgPack::GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, const ConstrainedOptimizationPack::MatrixSymHessianRelaxNonSing& H, BLAS_Cpp::Transp H_trans
	, const V& x, LinAlgPack::value_type b
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::Vp_StPtMtV;
	using SparseLinAlgPack::GenPermMatrixSlice;
	namespace GPMSIP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;

	const LinAlgPack::size_type
		no = H.G().rows(),  // number of original variables
		nr = H.M().rows(),      // number of relaxation variables
		nd = no + nr;        // total number of variables

	LinAlgPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans
		, BLAS_Cpp::rows( nd, nd, H_trans) );
	LinAlgPack::Vp_MtV_assert_sizes( BLAS_Cpp::cols( P.rows(), P.cols(), P_trans)
		, nd, nd, H_trans, x.size() );

	//
	// y = b*y + a * op(P) * H * x
	//
	// y = b*y + a * [op(P1)  op(P2) ] * [ G  0 ] * [ x1 ]
	//                                   [ 0  M ]   [ x2 ]
	//
	// =>
	//
	// y = b*y + a*op(P1)*G*x1 + a*op(P2)*H*x2
	//
	// For this to work op(P) must be sorted by column.
	//
	if( 	( P.ordered_by() == GPMSIP::BY_ROW && P_trans == BLAS_Cpp::no_trans )
	    || 	( P.ordered_by() == GPMSIP::BY_COL && P_trans == BLAS_Cpp::trans )
		||  ( P.ordered_by() == GPMSIP::UNORDERED ) )
	{
		// Call the default implementation
		H.MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b);
		return;
	}
	const LinAlgPack::Range1D
		o_rng(1,no),
		r_rng(no+1,no+nr);
	const SparseLinAlgPack::GenPermMatrixSlice
		P1 = ( P.is_identity() 
			   ? GenPermMatrixSlice(
				   P_trans == no_trans ? nd : no 
				   ,P_trans == no_trans ? no : nd
				   ,GenPermMatrixSlice::IDENTITY_MATRIX )
			   : P.create_submatrix(o_rng,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
			),
		P2 = ( P.is_identity()
			   ? GenPermMatrixSlice(
				   P_trans == no_trans ? nd : nr
				   ,P_trans == no_trans ? nr : nd
				   ,GenPermMatrixSlice::ZERO_MATRIX )
			   : P.create_submatrix(r_rng,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
			);
	const V
		x1 = x(o_rng),
		x2 = x(r_rng);
	// y = b*y
	LinAlgPack::Vt_S(y,b); // Handles b == 0.0 properly!

	// y += a*op(P1)*G*x1
	if( P1.nz() )
		SparseLinAlgPack::Vp_StPtMtV( y, a, P1, P_trans, H.G(), H_trans, x1, b );
	// y2 += a*op(P2)*M*x2
	if( P2.nz() )
		SparseLinAlgPack::Vp_StPtMtV( y, a, P2, P_trans, H.M(), H_trans, x2, 1.0 );
}

} // end namespace

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
	assert_initialized();
	return G_ptr_;
}

const MatrixSymWithOpFactorized& MatrixSymHessianRelaxNonSing::G() const
{
	assert_initialized();
	return *G_ptr_;
}

const SparseLinAlgPack::MatrixSymDiagonalStd& MatrixSymHessianRelaxNonSing::M() const
{
	assert_initialized();
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
	assert_initialized();
//	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // Uncomment for this default implementation
	Vp_StPtMtV_imp(y,a,P,P_trans,*this,H_trans,x,b);
}

void MatrixSymHessianRelaxNonSing::Vp_StPtMtV(
	VectorSlice* y, value_type a, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp H_trans, const SpVectorSlice& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StPtMtV(...SpVectorSlice...)" );
#endif
	assert_initialized();
//	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // Uncomment for this default implementation
	Vp_StPtMtV_imp(y,a,P,P_trans,*this,H_trans,x,b);
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
