// ///////////////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianSuperBasic.h

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_SUPER_BASIC_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_SUPER_BASIC_H

#include "QPSolverRelaxedQPSchur.h"

namespace ConstrainedOptimizationPack {

///
/** Implementation of initial KKT system for all variables initially fixed
 * and free where #Ko = B_RR#.
 *
 * In this implementation, #G# must support the \Ref{MatrixHessianSuperBasic}
 * interface.
 */
class QPSchurInitKKTSystemHessianSuperBasic
	: public QPSolverRelaxedQPSchur::InitKKTSystem 
{

	// ////////////////////////////////
	// Overridden from InitKKTSystem

	///
	/** Initialize the KKT system where the variables are initiallly 
	 * fixed and free and no constraints are in Ko.
	 *
	 * The Hessian for the QP without the relaxation #G# is represented as
	 * a \Ref{MatrixHessianSuperBasic} object and is:
     *
	 * #G = Q_R*B_RR*Q_R' + Q_R*op(B_RX)*Q_X' + Q_X*op(B_RX')*Q_R + Q_X*B_XX*Q_X'#
	 *
	 * If #G# does not support the interface #MatrixHessianSuperBasic# then
	 * an exception will be thrown.
	 *
	 * Given the above parts of #G#, define: #[nd,nd_R] = size(G.Q_R)# and
	 * #[nd,nd_X] = size(G.Q_X)#.  Then initial KKT system is defined as:
	 *
	 * #i_x_free[l-1] = (G.Q_R.begin()+l-1)->row_i(), l = 1...nd_R#\\
	 * #i_x_fixed[l-1] = (G.Q_X.begin()+l-1)->row_i(), l = 1...nd_X#\\
	 * #i_x_fixed[nd_X] = nd+1#\\
	 * #bnd_fixed[l-1] = G.bnd_fixed[l-1], l = 1...nd_X#\\
	 * #bnd_fixed[nd_X] = LOWER#\\
	 * #j_f_decomp[] = empty#\\
	 * #b_X[l-1] = { dL(i) if bnd_fixed[l-1] == LOWER or EQUALITY, dU(i) if bnd_fixed[l-1] == UPPER }#
	 * #, l = 1...nd_X (where i = i_x_fixed[l-1])#\\
	 * #b_X[nd_X] = etaL#\\
	 * #Ko = G.B_RR#\\
	 * #fo = - G.Q_R'*g - op(G.B_RX)*b_X(1:nd_X)#\\\
	 *
	 * Above, it is assumed that if #G.bnd_fixed[l-1] == EQUALITY#, that
	 * #dL(G.i_x_fixed[l-1]) == dU(G.i_x_fixed[l-1]# but this may not be
	 * inforced by this class.
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
		,i_x_free_t*          i_x_free
		,i_x_fixed_t*         i_x_fixed
		,bnd_fixed_t*         bnd_fixed
		,j_f_decomp_t*        j_f_decomp
		,Vector*              b_X
		,Ko_ptr_t*            Ko
		,Vector*              fo
		) const;

}; // end class QPSchurInitKKTSystemHessianSuperBasic

} // end namesapce ConstrainedOptimizationPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_SUPER_BASIC_H
