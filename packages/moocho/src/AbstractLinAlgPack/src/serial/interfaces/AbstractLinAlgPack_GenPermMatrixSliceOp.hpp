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

//@}

}	// end namespace SparseLinAlgPack

#endif   // GEN_PERM_MATRIX_SLICE_OP_H
