// //////////////////////////////////////////////////////////////////////////////////
// MatrixVarReductImplicit.h

#ifndef MATRIX_VAR_REDUCT_IMPLICIT_H
#define MATRIX_VAR_REDUCT_IMPLICIT_H

#include <vector>

#include "DecompositionSystemVarReduct.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "Misc/include/ref_count_ptr.h"
#include "Misc/include/ReleaseResource.h"

namespace ConstrainedOptimizationPack {

///
/** Implements #D = - inv(C) * N# for a variable reduction projection.
  *
  * This class is used to implement the adjoint factorization for the
  * null space matrix #Z = [ D; I ]#.  It implements #y = op(D)*x# through 
  * a referenced \Ref{DecompositionSystemVarReduct} object.
  *
  * The operations #y = op(D)*x# are implemented as:
  \begin{verbatim}
	y = D * x
	  = -inv(C) * (N * x)

	y = D' * x
	  = - N' * (inv(C') * x)
  \end{verbatim}
  *
  * The client can also setup a precomputed dense matrix #D = -inv(C)*N# which
  * may be used in many operations if convenient.
  *
  * This implementation is designed to deal efficiently with the case where matrix
  * vector multiplications will only be performed with subsets of rows of 
  * inv(C)*N or columns of N'*inv(C').  This primarily affects two types of operations:
  *
  * y = b*y + a*[-N'*inv(C')]*x       (x is a sparse vector)
  *
  * y = b*y + a*op(P)*[-inv(C)*N]*x   ( P has few nonzeros, x is sparse or dense )
  *
  * When D_dense has been set, then no sparse linear algebra is needed.  However, when
  * D_dense is not set then the needed rows of inv(C)*N are generated on the fly
  * and stored away for later use.  When \Ref{initialize}#(...)# is called then all
  * of these computed rows are discarded and they must be generated again.
  */
class MatrixVarReductImplicit : public MatrixWithOp
{
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		ResourceManagementPack::ReleaseResource>  release_resource_ptr_t;
	///
	MatrixVarReductImplicit();
	///
	/** Initialize the matrix object.
	 *
	 * @param  decomp_sys  [in] Pointer to the decomposition system object that
	 *                     represents #C# and #N#.  This must not be #NULL# and
	 *                     the object must not be altered while #this# matrix
	 *                     object is in use.
	 * @param  D_dense     [in] Pointer to a matrix #D = -inv(C)*N# already
	 *                     computed.  The matrix object #D_dense# will not be modifed
	 *                     by #this# and must not be altered while #this# matrix object
	 *                     is in use.  #D_dense == NULL# is allowed and #this#
	 *                     matrix object will just have to do without (okay).
	 * @param  release_resource_ptr
	 *                     [in] Points to a resource to be released when this
	 *                     matrix object is deleted or when \Ref{set_uninitialized}#()#
	 *                     is called.
	 */
	virtual void initialize(
		const DecompositionSystemVarReduct     *decomp_sys
		,const GenMatrixSlice                  *D_dense
		,const release_resource_ptr_t          &release_resource_ptr
		);
	///
	/** Set the matrix to uninitialized.
	 */
	virtual void set_uninitialized();
	///
	const DecompositionSystemVarReduct& decomp_sys() const;
	///
	const release_resource_ptr_t& release_resource_ptr() const;

	/** @name Overridden from Matrix. */
	//@{
	///
	size_type rows() const;
	///
	size_type cols() const;
	//@}

	/** @name Overridden from MatrixWithOp. */
	//@{
	///
	MatrixWithOp& operator=(const MatrixWithOp& M);
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorSlice& vs_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;
	//@}

public:

	// Not the users concern!
	typedef std::vector<value_type*>    InvCtN_rows_t;

private:

	// //////////////////////////
	// Private data members

	const DecompositionSystemVarReduct* decomp_sys_;
	GenMatrixSlice                      D_dense_;
	bool                                use_dense_mat_vec_;  // If true then D_dense_ will be used
	                                                         // for matrix vector multiplication if available
	release_resource_ptr_t              release_resource_ptr_;
	mutable InvCtN_rows_t               InvCtN_rows_;
	// InvCtN_rows_ keeps track of a set pointers of computed rows of inv(C)*N.
	// If D_dense_ is setup then InvCtN_rows_ is not necessary.  However, if
	// not then #InvCtN_rows_[j-1]# will be #!=NULL# if it points to the precomputed
	// jth row of #inv(C)*N# and will be #==NULL# if this row has not been computed
	// yet.  Each time initialize(...) is called, these rows are deallocated and
	// #InvCtN_rows_[j-1], j=1...this->rows()# is set to #NULL#.

	// //////////////////////////////////
	// Private member functions

	///
	void assert_initialized() const;

};	// end class MatrixVarReductImplicit

}	// end namespace ConstrainedOptimizationPack 

#endif	// MATRIX_VAR_REDUCT_IMPLICIT_H
