// ////////////////////////////////////////////////////////////////
// MatrixSymHessianRelaxNonSing.cpp
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

#include "ConstrainedOptimizationPack/src/MatrixSymHessianRelaxNonSing.h"
#include "AbstractLinAlgPack/src/SpVectorClass.h"
#include "AbstractLinAlgPack/src/GenPermMatrixSlice.h"
#include "AbstractLinAlgPack/src/VectorSpaceCompositeStd.h"
#include "AbstractLinAlgPack/src/LinAlgOpPack.h"
#include "profile_hack.h"
#include "ThrowException.h"

namespace {

/* ToDo: Finish updating below code!

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
	using SparseLinAlgPack::MatrixWithOp;
	namespace GPMSIP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;

	const LinAlgPack::size_type
		no = H.G().rows(),  // number of original variables
		nr = H.M().rows(),  // number of relaxation variables
		nd = no + nr;       // total number of variables

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
	LinAlgPack::Vt_S(y,b);
	// y += a*op(P1)*G*x1
	if( P1.nz() )
		SparseLinAlgPack::Vp_StPtMtV( y, a, P1, P_trans, H.G(), H_trans, x1, b );
	// y2 += a*op(P2)*M*x2
	if( P2.nz() )
		SparseLinAlgPack::Vp_StPtMtV( y, a, P2, P_trans, H.M(), H_trans, x2, 1.0 );
}

*/

} // end namespace

namespace ConstrainedOptimizationPack {

MatrixSymHessianRelaxNonSing::MatrixSymHessianRelaxNonSing()
	: vec_space_(NULL,0)
{}

MatrixSymHessianRelaxNonSing::MatrixSymHessianRelaxNonSing(
	const G_ptr_t         &G_ptr
	,const vec_mut_ptr_t  &M_diag_ptr
	,const space_ptr_t    &space
	)
	: vec_space_(NULL,0)
{
	initialize(G_ptr,M_diag_ptr,space);
}

void MatrixSymHessianRelaxNonSing::initialize(
	const G_ptr_t         &G_ptr
	,const vec_mut_ptr_t  &M_diag_ptr
	,const space_ptr_t    &space
	)
{
	namespace mmp = MemMngPack;
#ifdef _DEBUG
	const char err_msg_head[] = "MatrixSymHessianRelaxNonSing::initialize(...) : Error!";
	THROW_EXCEPTION(G_ptr.get()==NULL, std::invalid_argument, err_msg_head);
	THROW_EXCEPTION(M_diag_ptr.get()==NULL, std::invalid_argument, err_msg_head);
	const size_type G_rows = G_ptr->rows(), M_diag_dim = M_diag_ptr->dim();
	THROW_EXCEPTION(G_rows==0, std::invalid_argument, err_msg_head);
	THROW_EXCEPTION(M_diag_dim==0, std::invalid_argument, err_msg_head);
#endif
	if( space.get() ) {
#ifdef _DEBUG
		const size_type space_dim = space->dim();
		THROW_EXCEPTION(space_dim != G_rows + M_diag_dim, std::invalid_argument, err_msg_head);
#endif
		vec_space_ = space;
	}
	else {
		VectorSpace::space_ptr_t spaces[]
			= { mmp::rcp(&G_ptr->space_cols(),false), mmp::rcp(&M_diag_ptr->space(),false) };
		vec_space_ = mmp::rcp(new VectorSpaceCompositeStd( spaces, 2 ) );
	}
	G_ptr_ = G_ptr;
	M_.initialize(M_diag_ptr);
}
	
// Overridden from MatrixWithOp

const VectorSpace& MatrixSymHessianRelaxNonSing::space_cols() const
{
	assert_initialized();
	return *vec_space_;
}

bool MatrixSymHessianRelaxNonSing::Mp_StM(
	MatrixWithOp* C, value_type a, BLAS_Cpp::Transp H_trans
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Mp_StM(...)" );
#endif
	assert_initialized();
	return MatrixWithOp::Mp_StM(C,a,H_trans); // ToDo: Update below code!
/*
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	SparseLinAlgPack::Mp_StM( &(*C)(1,nG,1,nG), a, *G_ptr_, H_trans);
	SparseLinAlgPack::Mp_StM( &(*C)(nG+1,nG+nM,nG+1,nG+nM), a, M_, H_trans);
*/
}

void MatrixSymHessianRelaxNonSing::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp H_trans
	,const VectorWithOp& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StMtV(...VectorWithOp...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	AbstractLinAlgPack::Vt_S(y,b);
	AbstractLinAlgPack::Vp_StMtV( y->sub_view(1,nG).get(), a, *G_ptr_, H_trans, *x.sub_view(1,nG) );
	AbstractLinAlgPack::Vp_StMtV( y->sub_view(nG+1,nG+nM).get(), a, M_, H_trans, *x.sub_view(nG+1,nG+nM) );
}

void MatrixSymHessianRelaxNonSing::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp H_trans
	,const SpVectorSlice& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StMtV(...SpVectorSlice...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	AbstractLinAlgPack::Vt_S(y,b); // Takes care of b == 0.0 and y uninitialized
	AbstractLinAlgPack::Vp_StMtV( y->sub_view(1,nG).get(), a, *G_ptr_, H_trans, x(1,nG) );
	AbstractLinAlgPack::Vp_StMtV( y->sub_view(nG+1,nG+nM).get(), a, M_, H_trans, x(nG+1,nG+nM) );
}

void MatrixSymHessianRelaxNonSing::Vp_StPtMtV(
	VectorWithOpMutable* y, value_type a, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,BLAS_Cpp::Transp H_trans, const VectorWithOp& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StPtMtV(...VectorWithOp...)" );
#endif
	assert_initialized();
	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // Uncomment for this default implementation
/* ToDo: Update below code!
	Vp_StPtMtV_imp(y,a,P,P_trans,*this,H_trans,x,b);
*/
}

void MatrixSymHessianRelaxNonSing::Vp_StPtMtV(
	VectorWithOpMutable* y, value_type a, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,BLAS_Cpp::Transp H_trans, const SpVectorSlice& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Vp_StPtMtV(...SpVectorSlice...)" );
#endif
	assert_initialized();
	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,H_trans,x,b); // Uncomment for this default implementation
/* ToDo: Update below code!
	Vp_StPtMtV_imp(y,a,P,P_trans,*this,H_trans,x,b);
*/
}

// Overridden form MatrixSymWithOp

void MatrixSymHessianRelaxNonSing::Mp_StPtMtP(
	MatrixSymWithOp* S, value_type a
	,EMatRhsPlaceHolder dummy_place_holder
	,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,value_type b
	) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::Mp_StPtMtP(...)" );
#endif
	assert_initialized();

	MatrixSymWithOp::Mp_StPtMtP(S,a,dummy_place_holder,P,P_trans,b); // ToDo: Override when needed!
	return;
/* ToDo: Update below code!
	const LinAlgPack::size_type
		no = G().rows(),     // number of original variables
		nr = M().rows(),     // number of relaxation variables
		nd = no + nr;        // total number of variables

	LinAlgPack::Mp_MtM_assert_sizes( S->rows(), S->cols(), no_trans
									 , P.rows(), P.cols(), trans_not(P_trans)
									 , P.rows(), P.cols(), P_trans );
	LinAlgPack::Vp_V_assert_sizes( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans), nd );

	//
	// S = b*S + a * op(P)' * H * op(P)
	//
	// S = b*S + a * [op(P1)'  op(P2)' ] * [ G  0 ] * [ op(P1) ]
	//                                     [ 0  M ]   [ op(P2) ]
	//
	// =>
	//
	// S = b*S
	// S1 += op(P1)' * G * op(P1)
	// S2 += op(P2)' * M * op(P2)
	//
	// For this to work op(P) must be sorted by row.
	//
	if( 	( P.ordered_by() == GPMSIP::BY_ROW && P_trans == BLAS_Cpp::trans )
	    || 	( P.ordered_by() == GPMSIP::BY_COL && P_trans == BLAS_Cpp::no_trans )
		||  ( P.ordered_by() == GPMSIP::UNORDERED ) )
	{
		// Call the default implementation
		MatrixSymWithOp::Mp_StPtMtP(S,a,dummy_place_holder,P,P_trans,b);
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
			   : P.create_submatrix(o_rng,P_trans==no_trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
			),
		P2 = ( P.is_identity()
			   ? GenPermMatrixSlice(
				   P_trans == no_trans ? nd : nr
				   ,P_trans == no_trans ? nr : nd
				   ,GenPermMatrixSlice::ZERO_MATRIX )
			   : P.create_submatrix(r_rng,P_trans==no_trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
			);
	// S = b*S
	LinAlgPack::Mt_S( &tri_ele_gms(S->gms(),S->uplo()),b); // Handles b == 0.0 properly!

	// S1 += a*op(P1)'*G*op(P1)
	if( P1.nz() )
		SparseLinAlgPack::Mp_StPtMtP(
			&sym_gms( S->gms()(1,no,1,no), S->uplo() )
			, a, dummy_place_holder, G(), P1, P_trans );
	// S2 += a*op(P2)'*M*op(P2)
	if( P2.nz() )
		SparseLinAlgPack::Mp_StPtMtP(
			&sym_gms( S->gms()(no+1,nd,no+1,nd), S->uplo() )
			, a, dummy_place_holder, M(), P2, P_trans );
*/
}

// Overridden from MatrixWithOpNonsingular

void MatrixSymHessianRelaxNonSing::V_InvMtV(
	VectorWithOpMutable* y, BLAS_Cpp::Transp H_trans, const VectorWithOp& x
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::V_InvMtV(...VectorWithOp...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	AbstractLinAlgPack::V_InvMtV( y->sub_view(1,nG).get(), *G_ptr_, H_trans, *x.sub_view(1,nG) );
	AbstractLinAlgPack::V_InvMtV( y->sub_view(nG+1,nG+nM).get(), M_, H_trans, *x.sub_view(nG+1,nG+nM) );
}

void MatrixSymHessianRelaxNonSing::V_InvMtV(
	VectorWithOpMutable* y, BLAS_Cpp::Transp H_trans, const SpVectorSlice& x
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixSymHessianRelaxNonSing::V_InvMtV(...SpVectorSlice...)" );
#endif
	assert_initialized();
	const size_type
		nG = G_ptr_->rows(),
		nM = M_.rows();
	AbstractLinAlgPack::V_InvMtV( y->sub_view(1,nG).get(), *G_ptr_, H_trans, x(1,nG) );
	AbstractLinAlgPack::V_InvMtV( y->sub_view(nG+1,nG+nM).get(), M_, H_trans, x(nG+1,nG+nM) );
}

// private

void MatrixSymHessianRelaxNonSing::assert_initialized() const
{
	THROW_EXCEPTION(
		G_ptr_.get() == NULL, std::logic_error
		,"MatrixSymHessianRelaxNonSing::assert_initialized(): Error, Not initalized yet!" );
}

} // end namespace ConstrainedOptimizationPack
