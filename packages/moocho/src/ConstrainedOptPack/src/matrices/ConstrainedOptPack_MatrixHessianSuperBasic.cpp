// ////////////////////////////////////////////////////////
// MatrixHessianSuperBasic.cpp

#include "ConstrainedOptimizationPack/include/MatrixHessianSuperBasic.h"
#include "ConstrainedOptimizationPack/include/initialize_Q_R_Q_X.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/MatrixWithOpOut.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"

namespace ConstrainedOptimizationPack {

MatrixHessianSuperBasic::MatrixHessianSuperBasic()
	: n_(0)
{}

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

	// i_x_free
	if( 0 < n_R && n_R < n && i_x_free == NULL ) {
		throw std::invalid_argument(
			"MatrixHessianSuperBasic::initialize(...) : Error, "
			"i_x_free can not be NULL when 0 < n_R < n" );
	}
	// i_x_fixed
	if( 0 < n_X && n_X < n && i_x_fixed == NULL ) {
		throw std::invalid_argument(
			"MatrixHessianSuperBasic::initialize(...) : Error, "
			"i_x_fixed can not be NULL when 0 < n-n_R < n" );
	}
	// bnd_fixed
	if( 0 < n_X && bnd_fixed == NULL ) {
		throw std::invalid_argument(
			"MatrixHessianSuperBasic::initialize(...) : Error, "
			"bnd_fixed can not be NULL when 0 < n-n_R" );
	}
	// B_RR
	if(n_R > 0 ) {
		if( !B_RR_ptr.get() )
			throw std::invalid_argument(
				"MatrixHessianSuperBasic::initialize(...) : Error, "
				"B_RR_ptr.get() can not be NULL when n_R > 0" );
		Mp_M_assert_sizes( n_R, n_R, no_trans, B_RR_ptr->rows(), B_RR_ptr->cols(), no_trans );
	}
	// op(B_RX)
	if( n_R < n ) {
		if( B_RX_ptr.get() ) {
			Mp_M_assert_sizes( n_R, n_X, no_trans, B_RX_ptr->rows(), B_RX_ptr->cols(), B_RX_trans );
		}
	}
	// B_XX
	if( n_R < n ) {
		if( !B_XX_ptr.get() )
			throw std::invalid_argument(
				"MatrixHessianSuperBasic::initialize(...) : Error, "
				"B_XX_ptr.get() can not be NULL if n_R < n" );
		Mp_M_assert_sizes( n_X, n_X, no_trans, B_XX_ptr->rows(), B_XX_ptr->cols(), no_trans );
	}

	// Setup Q_R and Q_X and validate i_x_free[] and i_x_fixed[]
	Q_R_row_i_.resize(n_R);
	Q_R_col_j_.resize(n_R);
	Q_X_row_i_.resize(n_X);
	Q_X_col_j_.resize(n_X);
	bool test_setup = true;  // ToDo: Make this an input parameter!
	initialize_Q_R_Q_X(
		n_R,n_X,i_x_free,i_x_fixed,test_setup
		,n_R ? &Q_R_row_i_[0] : NULL
		,n_R ? &Q_R_col_j_[0] : NULL
		,&Q_R_
		,n_X ? &Q_X_row_i_[0] : NULL
		,n_X ? &Q_X_col_j_[0] : NULL
		,&Q_X_
	);

	// Setup bnd_fixed
	bnd_fixed_.resize(n_X);
	{for(size_type i = 0; i < n_X; ++i) bnd_fixed_[i] = bnd_fixed[i]; }

	// Set the rest of the arguments
	n_           = n;
	n_R_         = n_R;
	B_RR_ptr_    = B_RR_ptr;
	B_RX_ptr_    = B_RX_ptr;
	B_RX_trans_  = B_RX_trans;
	B_XX_ptr_    = B_XX_ptr;

}

// Overridden from Matrix

size_type MatrixHessianSuperBasic::rows() const
{
	return n_;
}

// Overridden from MatrixWithOp

void MatrixHessianSuperBasic::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp B_trans
	, const VectorSlice& x, value_type b
	) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using SparseLinAlgPack::V_MtV;
	using LinAlgOpPack::V_MtV;
	assert_initialized();
	LinAlgPack::Vp_MtV_assert_sizes( y->size(), n_, n_, B_trans, x.size() );
	if( n_ == n_R_ ) {
		//
		// B = Q_R*B_RR*Q_R'
		//
		// y = b*y + a*Q_R*B_RR*Q_R'*x
		//
		if( Q_R().is_identity() ) {
			SparseLinAlgPack::Vp_StMtV(y,a,*this->B_RR_ptr(),no_trans,x,b);
		}
		else {
			Vector Q_R_x;
			V_MtV( &Q_R_x, Q_R(), trans, x );
			SparseLinAlgPack::Vp_StPtMtV(y,a,Q_R(),no_trans,*this->B_RR_ptr(),no_trans,Q_R_x(),b);
		}
	}
	else if( n_R_ == 0 ) {
		//
		// B = Q_X *B_XX * Q_X'
		//
		assert(0); // ToDo: Implement this!
	}
	else {
		//
		// B = [ Q_R  Q_X  ] * [   B_RR      op(B_RX) ] * [ Q_R' ]
		//                     [ op(B_RX')      B_XX  ]   [ Q_X' ]
		//
		// y = b*y + a*op(B)*x
		//
		// y = b*y + a * [ Q_R  Q_X ] * [   B_RR      op(B_RX) ] * [ Q_R' ] * x
		//                              [ op(B_RX')      B_XX  ]   [ Q_X' ]
		//
		// y = b*y + a*Q_R*B_RR*x_R      + a*Q_R*op(B_RX)*x_X
		//         + a*Q_X*op(B_RX')*x_R + a*Q_X*B_XX*x_X
		// where:
		//     x_R = Q_R'*x
		//     x_X = Q_X'*x
		//
		SpVector
			x_R,
			x_X;
		// x_R = Q_R'*x
		V_MtV( &x_R, Q_R(), trans, x );
		// x_X = Q_X'*x
		V_MtV( &x_X, Q_X(), trans, x );
		// y = b*y + a*Q_R*B_RR*x_R
		SparseLinAlgPack::Vp_StPtMtV(
			y, a, Q_R(), no_trans, *B_RR_ptr(), no_trans, x_R(), b );
		// y += a*Q_R*op(B_RX)*x_X + a*Q_X*op(B_RX')*x_R
		if( B_RX_ptr().get() ) {
			SparseLinAlgPack::Vp_StPtMtV(
				y, a, Q_R(), no_trans, *B_RX_ptr(), B_RX_trans(), x_X() );
			SparseLinAlgPack::Vp_StPtMtV(
				y, a, Q_X(), no_trans, *B_RX_ptr(), trans_not(B_RX_trans()), x_R() );
		}
		// y += a*Q_X*B_XX*x_X
		SparseLinAlgPack::Vp_StPtMtV(
			y, a, Q_X(), no_trans, *B_XX_ptr(), no_trans, x_X() );
	}
}

void MatrixHessianSuperBasic::Vp_StMtV(
	VectorSlice* y, value_type a, BLAS_Cpp::Transp B_trans
	, const SpVectorSlice& x, value_type b
	) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using SparseLinAlgPack::V_MtV;
	using LinAlgOpPack::V_MtV;
	assert_initialized();
	LinAlgPack::Vp_MtV_assert_sizes( y->size(), n_, n_, B_trans, x.size() );
	if( n_ == n_R_ ) {
		//
		// B = Q_R*B_RR*Q_R'
		//
		// y = b*y + a*Q_R*B_RR*Q_R'*x
		//
		if( Q_R().is_identity() ) {
			SparseLinAlgPack::Vp_StMtV(y,a,*this->B_RR_ptr(),no_trans,x,b);
		}
		else {
			SpVector Q_R_x;
			SparseLinAlgPack::V_MtV( &Q_R_x, Q_R(), trans, x );
			SparseLinAlgPack::Vp_StPtMtV(y,a,Q_R(),no_trans,*this->B_RR_ptr(),no_trans,Q_R_x(),b);
		}
	}
	else if( n_R_ == 0 ) {
		//
		// B = Q_X *B_XX * Q_X'
		//
		assert(0); // ToDo: Implement this!
	}
	else {
		//
		// B = [ Q_R  Q_X  ] * [   B_RR      op(B_RX) ] * [ Q_R' ]
		//                     [ op(B_RX')      B_XX  ]   [ Q_X' ]
		//
		// y = b*y + a*op(B)*x
		//
		// y = b*y + a * [ Q_R  Q_X ] * [   B_RR      op(B_RX) ] * [ Q_R' ] * x
		//                              [ op(B_RX')      B_XX  ]   [ Q_X' ]
		//
		// y = b*y + a*Q_R*B_RR*x_R      + a*Q_R*op(B_RX)*x_X
		//         + a*Q_X*op(B_RX')*x_R + a*Q_X*B_XX*x_X
		// where:
		//     x_R = Q_R'*x
		//     x_X = Q_X'*x
		//
		SpVector
			x_R,
			x_X;
		// x_R = Q_R'*x
		V_MtV( &x_R, Q_R(), trans, x );
		// x_X = Q_X'*x
		V_MtV( &x_X, Q_X(), trans, x );
		// y = b*y + a*Q_R*B_RR*x_R
		SparseLinAlgPack::Vp_StPtMtV(
			y, a, Q_R(), no_trans, *B_RR_ptr(), no_trans, x_R(), b );
		// y += a*Q_R*op(B_RX)*x_X + a*Q_X*op(B_RX')*x_R
		if( B_RX_ptr().get() ) {
			SparseLinAlgPack::Vp_StPtMtV(
				y, a, Q_R(), no_trans, *B_RX_ptr(), B_RX_trans(), x_X() );
			SparseLinAlgPack::Vp_StPtMtV(
				y, a, Q_X(), no_trans, *B_RX_ptr(), trans_not(B_RX_trans()), x_R() );
		}
		// y += a*Q_X*B_XX*x_X
		SparseLinAlgPack::Vp_StPtMtV(
			y, a, Q_X(), no_trans, *B_XX_ptr(), no_trans, x_X() );
	}
}

value_type MatrixHessianSuperBasic::transVtMtV(
	const SpVectorSlice& x1, BLAS_Cpp::Transp B_trans
	, const SpVectorSlice& x2 ) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	assert_initialized();
	LinAlgPack::Vp_MtV_assert_sizes( x1.size(), rows(), cols(), B_trans, x1.size() );
	if( n_ == n_R_ ) {
		//
		// B = Q_R*B_RR*Q_R'
		//
		// a = x1'*Q_R*B_RR*Q_R'*x2
		//
		if( Q_R().is_identity() ) {
			return SparseLinAlgPack::transVtMtV( x1, *B_RR_ptr(), no_trans, x2 );
		}
		else {
			if( x1.overlap(x2) == LinAlgPack::SAME_MEM ) {
				SpVector Q_RT_x2;
				SparseLinAlgPack::V_MtV( &Q_RT_x2, Q_R(), trans, x2 );
				SpVectorSlice Q_RT_x2_slc = Q_RT_x2();
				return SparseLinAlgPack::transVtMtV(
					Q_RT_x2_slc, *B_RR_ptr(), no_trans, Q_RT_x2_slc );
 			}
			else {
				SpVector Q_RT_x2;
				SparseLinAlgPack::V_MtV( &Q_RT_x2, Q_R(), trans, x2 );
				SpVector Q_RT_x1;
				SparseLinAlgPack::V_MtV( &Q_RT_x1, Q_R(), trans, x1 );
				return SparseLinAlgPack::transVtMtV(
					Q_RT_x1(), *B_RR_ptr(), no_trans, Q_RT_x2() );
			}
		}
	}
	else if( n_R_ == 0 ) {
		//
		// B = Q_X *B_XX * Q_X'
		//
		assert(0); // ToDo: Implement this!
	}
	else {
		//
		// B = [ Q_R  Q_X  ] * [   B_RR      op(B_RX) ] * [ Q_R' ]
		//                     [ op(B_RX')      B_XX  ]   [ Q_X' ]
		//
		//
		// a = x1'*B*x2
		// =>
		// a = x1' * [ Q_R  Q_X  ] * [   B_RR      op(B_RX) ] * [ Q_R' ] * x2
		//                           [ op(B_RX')      B_XX  ]   [ Q_X' ]
		//
		// a = x1'*Q_R*B_RR*Q_R'*x2 + 2*x1'*Q_R*op(B_RX)*Q_X'*x2 + x1'*Q_X*B_XX*Q_X'*x2
		//
		if( x1.overlap(x2) == LinAlgPack::SAME_MEM ) {
			// a = x1'*Q_R*B_RR*Q_R'*x1 + 2*x1'*Q_R*op(B_RX)*Q_X'*x1 + x1'*Q_X*B_XX*Q_X'*x1
			SpVector Q_RT_x1;
			if( Q_R().nz() )
				SparseLinAlgPack::V_MtV( &Q_RT_x1, Q_R(), trans, x1 );
			SpVector Q_XT_x1;
			if( Q_X().nz() )
				SparseLinAlgPack::V_MtV( &Q_XT_x1, Q_X(), trans, x1 );
			SpVectorSlice Q_RT_x1_slc = Q_RT_x1();
			SpVectorSlice Q_XT_x1_slc = Q_XT_x1();
			return
				( Q_R().nz()
				  ? SparseLinAlgPack::transVtMtV(
					  Q_RT_x1_slc, *B_RR_ptr(), no_trans, Q_RT_x1_slc )
				  : 0.0
					)
				+ 2*(  B_RX_ptr().get() && Q_R().nz() && Q_X().nz()
					   ? SparseLinAlgPack::transVtMtV(
						   Q_RT_x1_slc, *B_RX_ptr(), B_RX_trans(), Q_XT_x1_slc )
					   : 0.0
					)
				+ ( Q_X().nz()
					? SparseLinAlgPack::transVtMtV(
						Q_XT_x1_slc, *B_XX_ptr(), no_trans, Q_XT_x1_slc )
					: 0.0
					);
		}
		else {
			assert(0); // ToDo: Implement this!
		}
	}
}

// Private

void MatrixHessianSuperBasic::assert_initialized() const
{
	if( !n_ )
		throw std::logic_error(
			"MatrixHessianSuperBasic::assert_initialized() : Error, "
			"The matrix is not initialized yet" );
}

} // end namespace ConstrainedOptimizationPack
