// ////////////////////////////////////////////////////////////////////////////////////////////////////
// MatrixIdentConcat.cpp

#include "ConstrainedOptimizationPack/include/MatrixIdentConcat.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/include/MatrixWithOpOut.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
}

namespace {

template<class V>
void mat_vec(
	LinAlgPack::VectorSlice                *y
	,LinAlgPack::value_type                a
	,LinAlgPack::value_type                alpha
	,const SparseLinAlgPack::MatrixWithOp  &D
	,BLAS_Cpp::Transp                      D_trans
	,const LinAlgPack::Range1D             &D_rng
	,const LinAlgPack::Range1D             &I_rng
	,BLAS_Cpp::Transp                      M_trans
	,const V                               &x
	,LinAlgPack::value_type                b
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	//
	// M = [ alpha*op(D) ] or [      I       ]
	//     [     I       ]    [ alpha*op(D)  ]
	//
	if( M_trans == no_trans ) {
		// 
		// y = b*y + a*M*x
		// =>
		// y(D_rng) = b*y(D_rng) + a*alpha*op(D)*x
		// y(I_rng) = b*y(I_rng) + a*x
		//
		LinAlgPack::VectorSlice
			y_D = (*y)(D_rng),
			y_I = (*y)(I_rng);
		// y(D_rng) = b*y(D_rng) + a*alpha*op(D)*x
		SparseLinAlgPack::Vp_StMtV( &y_D, a*alpha, D, D_trans, x, b );
		// y(I_rng) = b*y(I_rng) + a*x
		if( b == 0.0 )     y_I = 0.0;
		else if( b!=1.0 )  LinAlgOpPack::Vt_S(&y_I,b);
		LinAlgOpPack::Vp_StV( &y_I, a, x );
	}
	else {
		//
		// y = b*y + a*M'*x
		// =>
		// y = b*y + a*alpha*op(D')*x(D_rng) + a*x(I_rng)
		//
		const V
			x_D = x(D_rng),
			x_I = x(I_rng);
		// y = b*y + a*alpha*op(D')*x(D_rng) + a*x(I_rng)
		SparseLinAlgPack::Vp_StMtV( y, a*alpha, D, trans_not(D_trans), x_D, b );
		LinAlgOpPack::Vp_StV( y, a, x_I );
	}
}

} // end namespace

namespace ConstrainedOptimizationPack {

// Overridden from Matrix

size_type MatrixIdentConcat::rows() const
{
	const MatrixWithOp& D = this->D();
	return D.rows() + D.cols();
}

size_type MatrixIdentConcat::cols() const
{
	const MatrixWithOp& D = this->D();
	return BLAS_Cpp::cols(D.rows(),D.cols(),D_trans());
}

size_type MatrixIdentConcat::nz() const
{
	const MatrixWithOp& D = this->D();
	return D.nz() + BLAS_Cpp::cols(D.rows(),D.cols(),D_trans()); // D and I
}

// Overridden from MatrixWithOp

std::ostream& MatrixIdentConcat::output(std::ostream& out) const
{
	const Range1D           D_rng   = this->D_rng();
	const BLAS_Cpp::Transp  D_trans = this->D_trans();
	if( D_rng.lbound() == 1 ) {
		if( D_trans == BLAS_Cpp::no_trans )
			out << "[ alpha*D; I ]";
		else
			out << "[ alpha*D'; I ]";
	}
	else {
		if( D_trans == BLAS_Cpp::no_trans )
			out << "[ I; alpha*D ]";
		else
			out << "[ I; alpha*D' ]";
	}
	out << "\nalpha = " << this->alpha();
	out << "\nD =\n" << D();
	return out;
}

void MatrixIdentConcat::Vp_StMtV(VectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
	, const VectorSlice& x, value_type b) const
{
	LinAlgPack::Vp_MtV_assert_sizes( y->size(), rows(), cols(), M_trans, x.size() );
	mat_vec( y, a, alpha(), D(), D_trans(), D_rng(), I_rng(), M_trans, x, b );
}

void MatrixIdentConcat::Vp_StMtV(VectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b) const
{
	LinAlgPack::Vp_MtV_assert_sizes( y->size(), rows(), cols(), M_trans, x.size() );
	mat_vec( y, a, alpha(), D(), D_trans(), D_rng(), I_rng(), M_trans, x, b );
}

} // end namespace ConstrainedOptimizationPack
