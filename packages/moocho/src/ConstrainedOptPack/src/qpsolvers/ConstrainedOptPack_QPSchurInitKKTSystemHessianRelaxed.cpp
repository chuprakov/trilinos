// /////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianRelaxed.cpp
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

#include "ConstrainedOptPack_QPSchurInitKKTSystemHessianRelaxed.hpp"
#include "ConstrainedOptPack_MatrixSymHessianRelaxNonSing.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Midynamic_cast_verbose.h"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StV;
	using AbstractLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptPack {

void QPSchurInitKKTSystemHessianRelaxed::initialize_kkt_system(
	const DVectorSlice&    g
	,const MatrixOp&  G
	,value_type           etaL
	,const SpVectorSlice& dL
	,const SpVectorSlice& dU
	,const MatrixOp*  F
	,BLAS_Cpp::Transp     trans_F
	,const DVectorSlice*   f
	,const DVectorSlice&   d
	,const SpVectorSlice& nu
	,size_type*           n_R
	,i_x_free_t*          i_x_free
	,i_x_fixed_t*         i_x_fixed
	,bnd_fixed_t*         bnd_fixed
	,j_f_decomp_t*        j_f_decomp
	,DVector*              b_X
	,Ko_ptr_t*            Ko
	,DVector*              fo
	) const
{
	using BLAS_Cpp::trans;

	// Validate type of and convert G
	const MatrixSymHessianRelaxNonSing
		*G_relax_ptr = dynamic_cast<const MatrixSymHessianRelaxNonSing*>(&G);

	if( G_relax_ptr == NULL ) {
		init_kkt_full_.initialize_kkt_system(
			g,G,etaL,dL,dU,F,trans_F,f,d,nu,n_R,i_x_free,i_x_fixed,bnd_fixed
			,j_f_decomp,b_X,Ko,fo);
		return;
	}

	const MatrixSymHessianRelaxNonSing
		&G_relax = *G_relax_ptr;

	// get some stuff
	const MatrixSymWithOpFactorized
		&G_orig = G_relax.G(),
		&M      = G_relax.M();
	const size_type
		nd = g.size(),
		no = G_orig.rows(),
		nr = M.rows();
	assert( no + nr == nd );

	// Setup output arguments

	// n_R = nd_R
	*n_R = no;
	// i_x_free.size() == 0 and i_x_free is implicitly identity
	i_x_free->resize(no);
	{for(size_type l = 1; l <= no; ++l ) {
		(*i_x_free)[l-1] = l;
	}}
	// i_x_fixed[]
	i_x_fixed->resize(nr+1);
	if(nr) {
		// i_x_fixed[l-1] = no + l, l = 1...nr
		for( size_type l = 1; l <= nr; ++l )
			(*i_x_fixed)[l-1] = no+l;
	}
	(*i_x_fixed)[nr] = nd+1; // extra relaxation is always initially active
	// bnd_fixed[]
	bnd_fixed->resize(nr+1);
	if(nr) {
		// bnd_fixed[l-1] = LOWER, l = 1...nr
		std::fill_n( bnd_fixed->begin(), nr, LOWER );
	}
	(*bnd_fixed)[nr] = LOWER; // relaxation is always initially active
	// j_f_decomp[]
	j_f_decomp->resize(0);
	// b_X
	b_X->resize(nr+1);
	if(nr) {
		// b_X[l-1] = dL(no+l), l = 1...nr
		LinAlgOpPack::assign( &(*b_X)(1,nr), dL(no+1,no+nr) );
	}
	(*b_X)[nr] = etaL; // relaxation is always initially active
	// Ko = G.G
	*Ko = G_relax.G_ptr(); // now B_RR is a shared object
	// fo = - *g(1:no)
	LinAlgOpPack::V_StV( fo, -1.0, g(1,no) );

}

} // end namesapce ConstrainedOptPack
