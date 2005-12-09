// ///////////////////////////////////////////////////////////////////
// ConstrainedOptPack_QPSchurInitKKTSystemHessianFixedFree.hpp
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

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FIXED_FREE_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FIXED_FREE_H

#include "ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp"

namespace ConstrainedOptPack {

///
/** Implementation of initial KKT system using the Hessian for the free
 * variables only.
 *
 * In this implementation, #G# must support the #MatrixSymOp#
 * interface.
 */
class QPSchurInitKKTSystemHessianFixedFree
	: public QPSolverRelaxedQPSchur::InitKKTSystem 
{
public:

	// ////////////////////////////////
	// Overridden from InitKKTSystem

	///
	/** Initialize the KKT system where initially fixed variables are removed and
	 * no equality constraints are included in Ko.
	 *
	 * For this implementation:
	 *
	 * ToDo: Finish documentation!
	 */
	void initialize_kkt_system(
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
		) const;

}; // end class QPSchurInitKKTSystemHessianFixedFree

} // end namesapce ConstrainedOptPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FIXED_FREE_H
