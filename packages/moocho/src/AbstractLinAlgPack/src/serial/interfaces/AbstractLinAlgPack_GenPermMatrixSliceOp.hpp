// /////////////////////////////////////////////////////////
// GenPermMatrixSliceOp.h

#ifndef GEN_PERM_MATRIX_SLICE_OP_H
#define GEN_PERM_MATRIX_SLICE_OP_H

#include "GenPermMatrixSlice.h"

namespace SparseLinAlgPack {

/** @name Operations for GenPermMatrixSlice.
  *
  * ToDo: Finish documentation!
  */
//@{

///
/** sv_lhs = alpha * op(P_rhs1) * vs_rhs2.
  * 
  * This function will resize the sparse vector lhs and add only the
  * nonzero elements in the rhs.
  * 
  * If op(P_rhs1) is sorted by row (i.e. op(P_rhs1) = P_rhs1 sorted by row
  * or op(P_rhs1) = P_rhs1' sorted by column) then sv_lhs->assume_sorted(true)
  * is called.
  * 
  * This function will execute in O(P_rhs1.nz()) time.
  */ 
void V_StMtV( SpVector* sv_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
	, BLAS_Cpp::Transp P_rhs1_trans, const VectorSlice& vs_rhs2 );

inline
/// sv_lhs = op(P_rhs1) * vs_rhs2
void V_MtV( SpVector* sv_lhs, const GenPermMatrixSlice& P_rhs1
	, BLAS_Cpp::Transp P_rhs1_trans, const VectorSlice& vs_rhs2 )
{
	V_StMtV(sv_lhs,1.0,P_rhs1,P_rhs1_trans,vs_rhs2);
}

///
/** sv_lhs = alpha * op(P_rhs1) * sv_rhs2.
  * 
  * This function will resize the sparse vector lhs and add only the
  * nonzero elements in the rhs.
  * 
  * If op(P_rhs1) is sorted by row (i.e. op(P_rhs1) = P_rhs1 sorted by row
  * or op(P_rhs1) = P_rhs1' sorted by column) then sv_lhs->assume_sorted(true)
  * is called.
  * 
  * Let's assume that op(P_rhs1) is sorted by row and sv_rhs2 is also sorted.
  * In this case a linear search will have to be performed to match up
  * elements.  P_rhs1 will be iterated through sequentially and
  * the corresponding nonzero element in sv_rhs2 searched for (binary search).
  * Therefore, the runtime in this case will be:
  *
  *   O( P_rhs1.nz() * log(sv_rhs2.nz()) )
  *
  * If P_rhs1 and sv_rhs2 are unsorted, then the runtime will be:
  * 
  *   O( P_rhs1.nz() * sv_rhs2.nz() )
  *
  * If op(P_rhs1) is sorted by column and sv_rhs2 is sorted then the runtime
  * will be:
  * 
  *   O( max( P_rhs1.nz(), sv_rhs2.nz() ) )
  * 
  * Of course if op(P_rhs1) is not sorted by row then the output vector will
  * not be assumed sorted.
  */ 
void V_StMtV( SpVector* sv_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
	, BLAS_Cpp::Transp P_rhs1_trans, const SpVectorSlice& sv_rhs2 );

inline
/// sv_lhs = op(P_rhs1) * sv_rhs2
void V_MtV( SpVector* sv_lhs, const GenPermMatrixSlice& P_rhs1
	, BLAS_Cpp::Transp P_rhs1_trans, const SpVectorSlice& sv_rhs2 )
{
	V_StMtV(sv_lhs,1.0,P_rhs1,P_rhs1_trans,sv_rhs2);
}


///
/** sv_lhs += alpha * op(P_rhs1) * vs_rhs2.
  * 
  * This function will not resize the sparse vector rhs and will add
  * new elements for the nonzero elements in the rhs.  Therefore it is
  * up to the client to ensure that there is sufficient storage for
  * these elements.  This function
  * will not check to see if elements with duplicate indices are added.
  * It is up to the client to determine that.
  * 
  * This function will execute in O(P_rhs1.nz()) time.
  */ 
void Vp_StMtV( SpVector* sv_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
	, BLAS_Cpp::Transp P_rhs1_trans, const VectorSlice& vs_rhs2 );

inline
///
/** sv_lhs += op(P_rhs1) * vs_rhs2.
  */ 
void Vp_MtV( SpVector* sv_lhs, const GenPermMatrixSlice& P_rhs1
	, BLAS_Cpp::Transp P_rhs1_trans, const VectorSlice& vs_rhs2 )
{
	Vp_StMtV(sv_lhs,1.0,P_rhs1,P_rhs1_trans,vs_rhs2);
}

/// vs_lhs = alpha * op(P_rhs1) * vs_rhs2 + beta * vs_lhs
void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
	, BLAS_Cpp::Transp P_rhs1_trans, const VectorSlice& vs_rhs2, value_type beta = 1.0);

/// vs_lhs = alpha * op(P_rhs1) * sv_rhs2 + beta * vs_lhs
void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, const GenPermMatrixSlice& P_rhs1
	, BLAS_Cpp::Transp P_rhs1_trans, const SpVectorSlice& sv_rhs2, value_type beta = 1.0);

///
/** Find the intersection between two GenPermMatrixSlice objects.
 *
 * This subroutine has two modes.  In the first mode (#Q_max_nz == 0#) it just
 * computes the number of nonzero entries in the matrix:
 *
 * #Q = op(P1)*op(P1)#
 *
 * In the second mode (#Q_max_nz > 0 && Q_row_i != NULL && Q_col_j != NULL#) it also
 * computes the row and column arrays for the resultant matrix #Q#.  In addition
 * if #Q != NULL# then a GenPermMatrixSlice object will be setup with as:
 *
 * #Q->initialize( op(P1).rows(), op(P2).cols(), Q_nz, 0, 0, Q_ordered_by, Q_row_i, Q_col_j, false )#
 *
 * Above #Q_ordered_by# will be determined of the fly.
 *
 * This operation might require  O(min(op(P1).cols(),op(P2).rows()) temporary storage
 * but will always be executed in:
 *
 * O(P1.nz()) + O(P2.nz()) + O(min(op(P1).cols(),op(P2).rows())
 *
 * @param  P1         [in] Right hand side permutation matrix.
 * @param  P1_trans   [in] If no_trans then op(P1) = P1 otherwise op(P1) = P1'.
 * @param  P1         [in] Left hand side permutation matrix.
 * @param  P1_trans   [in] If no_trans then op(P2) = P2 otherwise op(P2) = P2'.
 * @param  Q_nz       [out] On return will contain the number of nonzeros in the
 *                    resultant matrix #Q#
 * @param  Q_max_nz   [in] If #Q_max_nz > 0# then the resultant row #Q_row_i# and column #Q_col_j#
 *                    indices will be set.    If it turns out that #Q_nz# will be larger than
 *                    #Q_max_nz# then the exception #std::length_error# will be thrown and #Q_row_i#
 *                    and #Q_col_j# may be left in an inconsistent state. If #Q_max_nz == 0# then the
 *                    rest of the return arguments are ignored and the resultant matrix will not be returned.
 * @param  Q_row_i    [out] Array (length #Q_max_nz#): If #Q_max_nz > 0# then out return #Q_row_i[k], k=0,,,Q_nz-1#
 *                    will contain the row indices for the resultant matrix #Q#.  If #Q_max_nz == 0# then
 *                    #Q_row_i# can be #NULL#.
 * @param  Q_row_i    [out] Array (length #Q_max_nz#): If #Q_max_nz > 0# then out return #Q_col_j[k], k=0,,,Q_nz-1#
 *                    will contain the column indices for the resultant matrix #Q#.  If #Q_max_nz == 0# then
 *                    #Qd_col_j# can be #NULL#.
 * @param  Q          [out] If #Q_max_nz > 0 && Q != NULL# then #Q# will be initialized as described above.
 *                    It is allowd for #Q == NULL#.
 */
void intersection(
	const GenPermMatrixSlice     &P1
	,BLAS_Cpp::Transp            P1_trans
	,const GenPermMatrixSlice    &P2
	,BLAS_Cpp::Transp            P2_trans
	,size_type                   *Q_nz
	,const size_type             Q_max_nz     = 0
	,size_type                   Q_row_i[]    = NULL
	,size_type                   Q_col_j[]    = NULL
	,GenPermMatrixSlice          *Q           = NULL
	);

//@}

}	// end namespace SparseLinAlgPack

#endif   // GEN_PERM_MATRIX_SLICE_OP_H
