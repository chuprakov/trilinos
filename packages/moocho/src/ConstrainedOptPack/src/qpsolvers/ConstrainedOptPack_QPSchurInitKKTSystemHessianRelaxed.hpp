// ///////////////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianRelaxed.h

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_RELAXED_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_RELAXED_H

#include "QPSolverRelaxedQPSchur.h"
#include "QPSchurInitKKTSystemHessianFull.h"

namespace ConstrainedOptimizationPack {

///
/** Implementation of initial KKT system where all original variables
 * are free and all the relaxation variables are fixed.
 *
 * In this implementation, #G# should support the \Ref{MatrixSymHessianRelaxNonSing}
 * interface.   Otherwise, it will try the \Ref{MatrixSymWithOpFactorized} interface
 * using the base implementation of \Ref{QPSchurInitKKTSystemHessianFull}.
 */
class QPSchurInitKKTSystemHessianRelaxed
	: public QPSchurInitKKTSystemHessianFull
{
public:

	// ////////////////////////////////
	// Overridden from InitKKTSystem

	///
	/** Initialize the KKT system where the original variables are initiallly 
	 * free and all the relaxation variables are fixed and their are no
	 * constraints in Ko.
	 *
	 * The Hessian for the QP without the relaxation #G# is represented as
	 * a \Ref{MatrixSymHessianRelaxNonSing} object and is:
     \begin{verbatim}
	 G = [ G_orig     ]
	     [          M ]
	 \end{verbatim}
	 * If #G# does not support the interface #MatrixSymHessianRelaxNonSing# then
	 * the function #QPSchurInitKKTSystemHessianFull::initialize_kkt_system(...)#
	 * will be called.
	 *
	 * Given the above parts of #G#, define: #[no,no] = size(G.G)# and
	 * #[nr,nr] = size(G.M)#.  Then initial KKT system is defined as:
	 *
	 * #n_R = no#\\
	 * #i_x_free.size() == 0# and #i_x_free is implicitly identity#\\
	 * #i_x_fixed[l-1] = no + l, l = 1...nr#\\
	 * #i_x_fixed[nr] = no+nr+1#\\
	 * #bnd_fixed[l-1] = LOWER, l = 1...nr#\\
	 * #bnd_fixed[nr] = LOWER#\\
	 * #j_f_decomp[] = empty#\\
	 * #b_X[l-1] = dL(no+l), l = 1...nr#\\
	 * #b_X[nr] = etaL#\\
	 * #Ko = G.G#\\
	 * #fo = - g(1:no)#\\\
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
		) const;

private:
	QPSchurInitKKTSystemHessianFull  init_kkt_full_;

}; // end class QPSchurInitKKTSystemHessianRelaxed

} // end namesapce ConstrainedOptimizationPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_RELAXED_H
