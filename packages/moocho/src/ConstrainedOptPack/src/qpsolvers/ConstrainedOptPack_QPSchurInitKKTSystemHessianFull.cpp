// /////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianFull.cpp
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

#include "ConstrainedOptimizationPack/src/QPSchurInitKKTSystemHessianFull.hpp"
#include "AbstractLinAlgPack/src/MatrixSymWithOpNonsingular.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"
#include "SparseLinAlgPack/src/VectorDenseEncap.hpp"
#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "dynamic_cast_verbose.hpp"

namespace ConstrainedOptimizationPack {

void QPSchurInitKKTSystemHessianFull::initialize_kkt_system(
	const VectorWithOp    &g
	,const MatrixWithOp   &G
	,value_type           etaL
	,const VectorWithOp   *dL
	,const VectorWithOp   *dU
	,const MatrixWithOp   *F
	,BLAS_Cpp::Transp     trans_F
	,const VectorWithOp   *f
	,const VectorWithOp   *d
	,const VectorWithOp   *nu
	,size_type            *n_R
	,i_x_free_t           *i_x_free
	,i_x_fixed_t          *i_x_fixed
	,bnd_fixed_t          *bnd_fixed
	,j_f_decomp_t         *j_f_decomp
	,DVector               *b_X
	,Ko_ptr_t             *Ko
	,DVector               *fo
	) const
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgOpPack::V_mV;

	// Validate type of and convert G
	const MatrixSymWithOpNonsingular&
		G_sym = dyn_cast<const MatrixSymWithOpNonsingular>(G);

	const size_type nd = g.dim();
	
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
	*Ko = mmp::rcp(&G_sym,false); // Not dynamically allocated so don't delete!
	// fo = -g
	V_mV(fo,VectorDenseEncap(g)());
}

} // end namesapce ConstrainedOptimizationPack
