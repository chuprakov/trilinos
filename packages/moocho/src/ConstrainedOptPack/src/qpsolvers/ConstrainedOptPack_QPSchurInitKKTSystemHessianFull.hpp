// ///////////////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianFull.h
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

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H

#include "QPSolverRelaxedQPSchur.h"

namespace ConstrainedOptimizationPack {

///
/** Implementation of initial KKT system for all variables initially free
 * and <tt>Ko = G</tt>.
 *
 * In this implementation, \c G must support the \c MatrixSymWithOpNonsingular
 * interface.  Using this initial KKT system essentially make QPSchur use
 * a range space approach (Nocedal & Wright, 1999) for factorizing the KKT
 * system for the current active set.
 */
class QPSchurInitKKTSystemHessianFull
	: public QPSolverRelaxedQPSchur::InitKKTSystem 
{
public:

	/** @name Overridden from InitKKTSystem */
	//@{

	///
	/** Initialize the KKT system where all variables (except the relaxation variable)
	 * are initially free and no constraints are in Ko.
	 *
	 * For this implementation:
	 *
	 * <tt>n_R = nd</tt>\\
	 * <tt>i_x_free = emply (it is identity)</tt>\\
	 * <tt>i_x_fixed[0] = nd+1</tt>\\
	 * <tt>bnd_fixed[0] = LOWER</tt>\\
	 * <tt>j_f_decomp[] = empty</tt>\\
	 * <tt>b_X = etaL</tt>\\
	 * <tt>Ko = G</tt>\\
	 * <tt>fo = -g</tt>\\
	 */
	void initialize_kkt_system(
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
		,Vector               *b_X
		,Ko_ptr_t             *Ko
		,Vector               *fo
		) const;

	//@}

}; // end class QPSchurInitKKTSystemHessianFull

} // end namesapce ConstrainedOptimizationPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H
