// ///////////////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianFull.h

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H

#include "QPSolverRelaxedQPSchur.h"

namespace ConstrainedOptimizationPack {

///
/** Implementation of initial KKT system for all variables initially free
 * and Ko = G.
 *
 * In this implementation, #G# must support the #MatrixSymWithOpFactorized#
 * interface.  Using this initial KKT system essentially make QPSchur use
 * a Range Space approach (Nocedal & Wright, 1999) for factorizing the KKT
 * system for the current active set.
 */
class QPSchurInitKKTSystemHessianFull
	: public QPSolverRelaxedQPSchur::InitKKTSystem 
{

	// ////////////////////////////////
	// Overridden from InitKKTSystem

	///
	/** Initialize the KKT system where all variables (except the relaxation variable)
	 * are initially free and no constraints are in Ko.
	 *
	 * For this implementation:
	 *
	 * #n_R = nd#\\
	 * #i_x_free = emply (it is identity)#\\
	 * #i_x_fixed[0] = nd+1#\\
	 * #bnd_fixed[0] = LOWER#\\
	 * #j_f_decomp[] = empty#\\
	 * #b_X = etaL#\\
	 * #Ko = G#\\
	 * #fo = -g#\\
	 */
	void initialize_kkt_system(
		const VectorSlice&    g
		,const MatrixWithOp&  G
		,value_type           etaL
		,const SpVectorSlice& dL
		,const SpVectorSlice& dU
		,const MatrixWithOp*  F
		,BLAS_Cpp::Transp     trans_F
		,const VectorSlice*   f
		,size_type*           n_R
		,i_x_free_t*          i_x_free
		,i_x_fixed_t*         i_x_fixed
		,bnd_fixed_t*         bnd_fixed
		,j_f_decomp_t*        j_f_decomp
		,Vector*              b_X
		,Ko_ptr_t*            Ko
		,Vector*              fo
		) const;

}; // end class QPSchurInitKKTSystemHessianFull

} // end namesapce ConstrainedOptimizationPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H
