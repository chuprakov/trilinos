// /////////////////////////////////////
// MatrixHessianRelaxed.cpp

#include <assert.h>

#include "ConstrainedOptimizationPack/include/MatrixHessianRelaxed.h"
#include "SparseLinAlgPack/include/MatrixSymWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptimizationPack {

MatrixHessianRelaxed::MatrixHessianRelaxed()
	:
		n_(0)
		,H_(NULL)
		,bigM_(0.0)
{}

void MatrixHessianRelaxed::initialize(
	  const MatrixSymWithOp	&H
	, value_type			bigM
	)
{
	n_	= H.rows();
	H_	= &H;
	bigM_	= bigM;
}

// Overridden from Matrix

size_type MatrixHessianRelaxed::rows() const
{
	return n_ + 1;
}

// Overridden from MatrixWithOp

MatrixWithOp& MatrixHessianRelaxed::operator=(const MatrixWithOp& m)
{
	// ToDo: Finish this!
	assert(0);
	return *this;
}

void MatrixHessianRelaxed::Vp_StMtV(
	  VectorSlice* y, value_type a, BLAS_Cpp::Transp G_trans
	, const VectorSlice& x, value_type b) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using SparseLinAlgPack::Vp_StMtV;
	
	LinAlgOpPack::Vp_MtV_assert_sizes(y->size(),rows(),cols(),G_trans,x.size());

	//
	// y = b*y + a * M * x
	// 
	//   = b*y + a * [ H  0    ] * [ x1 ]
	//               [ 0  bigM ]   [ x2 ]
	//               
	// =>              
	//               
	// y1 = b*y1 + a*H*x1
	// 
	// y2 = b*y2 + bigM * x2
	//

	VectorSlice
		y1 = (*y)(1,n_);
	value_type
		&y2 = (*y)(n_+1);
	const VectorSlice
		x1 = x(1,n_);
	const value_type
		x2 = x(n_+1);

	// y1 = b*y1 + a*H*x1
	Vp_StMtV( &y1, a, *H_, no_trans, x1, b );

	// y2 = b*y2 + bigM * x2
	if( b == 0.0 )
		y2 = bigM_ * x2;
	else
		y2 = b*y2 + bigM_ * x2;
	
}

}	// end namespace ConstrainedOptimizationPack
