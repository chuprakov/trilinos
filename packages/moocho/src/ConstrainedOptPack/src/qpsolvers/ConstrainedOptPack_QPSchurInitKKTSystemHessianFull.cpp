// /////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianFull.cpp

#include "ConstrainedOptimizationPack/include/QPSchurInitKKTSystemHessianFull.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpFactorized.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace ConstrainedOptimizationPack {

void QPSchurInitKKTSystemHessianFull::initialize_kkt_system(
	const VectorSlice&    g
	,const MatrixWithOp&  G
	,value_type           etaL
	,const SpVectorSlice& dL
	,const SpVectorSlice& dU
	,const MatrixWithOp*  F
	,BLAS_Cpp::Transp     trans_F
	,const VectorSlice*   f
	,const VectorSlice&   d
	,const SpVectorSlice& nu
	,size_type*           n_R
	,i_x_free_t*          i_x_free
	,i_x_fixed_t*         i_x_fixed
	,bnd_fixed_t*         bnd_fixed
	,j_f_decomp_t*        j_f_decomp
	,Vector*              b_X
	,Ko_ptr_t*            Ko
	,Vector*              fo
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgOpPack::V_mV;

	// Validate type of and convert G
#ifdef _WINDOWS
	const MatrixSymWithOpFactorized&
		G_sym = dynamic_cast<const MatrixSymWithOpFactorized&>(G);
#else
	const MatrixSymWithOpFactorized&
		G_sym = dyn_cast<const MatrixSymWithOpFactorized>(G);
#endif

	const size_type nd = g.size();
	
	// n_R
	*n_R = nd;
	// i_x_free[i-1] = i, i = 1...nd
	i_x_free->resize(0);
	// i_x_fixed[0] = nd+1
	i_x_fixed->resize(1);
	(*i_x_fixed)[0] = nd+1;
	// bnd_fixed[0] = LOWER
	bnd_fixed->resize(1);
	(*bnd_fixed)[0] = LOWER;
	// j_f_decomp[] = empty
	j_f_decomp->resize(0);
	// b_X = etaL
	b_X->resize(1);
	(*b_X)[0] = etaL;
	// Ko = G
	*Ko = &G_sym;
	Ko->release(); // Not dynamically allocated so don't delete!
	// fo = -g
	V_mV(fo,g);
}

} // end namesapce ConstrainedOptimizationPack
