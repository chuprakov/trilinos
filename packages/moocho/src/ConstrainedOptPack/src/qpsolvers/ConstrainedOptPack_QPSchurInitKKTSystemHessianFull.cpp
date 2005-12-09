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

#include "ConstrainedOptPack_QPSchurInitKKTSystemHessianFull.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace ConstrainedOptPack {

void QPSchurInitKKTSystemHessianFull::initialize_kkt_system(
	const Vector    &g
	,const MatrixOp   &G
	,value_type           etaL
	,const Vector   *dL
	,const Vector   *dU
	,const MatrixOp   *F
	,BLAS_Cpp::Transp     trans_F
	,const Vector   *f
	,const Vector   *d
	,const Vector   *nu
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
	using Teuchos::dyn_cast;
	using LinAlgOpPack::V_mV;

	// Validate type of and convert G
	const MatrixSymOpNonsing&
		G_sym = dyn_cast<const MatrixSymOpNonsing>(G);

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
	*Ko = Teuchos::rcp(&G_sym,false); // Not dynamically allocated so don't delete!
	// fo = -g
	V_mV(fo,VectorDenseEncap(g)());
}

} // end namesapce ConstrainedOptPack
