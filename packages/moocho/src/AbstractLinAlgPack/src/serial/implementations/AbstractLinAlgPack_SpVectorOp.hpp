// //////////////////////////////////////////////////////////////////////
// SpVectorOp.h

#ifndef SP_VECTOR_OP_H
#define SP_VECTOR_OP_H

#include "SparseElement.h"
#include "SparseVectorSliceOp.h"

namespace SparseLinAlgPack {

///
/** Add elements from a dense vector to a sparse vector.
 *
 * Here sv_lhs is not resized and only elements are added.
 * The purpose of this function is to add elements from
 * a dense vector to a sparse vector.
 *
 * Postconditions:\begin{itemize}
 * \item sv_lhs->nz() == vs_rhs->size() + sv_lhs_before->size()
 * \item [sv_lhs_before->is_sorted() || sv_lhs_before->nz() == 0] sv_lhs->is_sorted() == true
 * \end{itemize}
 */
void add_elements( SpVector* sv_lhs, const VectorSlice& vs_rhs, size_type offset = 0 );

///
/** Add elements from a sparse vector to another sparse vector.
 *
 * Here sv_lhs is not resized and only elements are added.
 * The purpose of this function is to add elements from
 * a sparse vector to a sparse vector.
 *
 * Postconditions:\begin{itemize}
 * \item sv_lhs->nz() == sv_rhs->nz() + sv_lhs_before->size()
 * \item [(sv_lhs_before->is_sorted() || sv_lhs_before->nz() == 0)
 *        && (sv_rhs.is_sorted() || sv_rhs.nz() == 0)] sv_lhs->is_sorted() == true
 * \end{itemize}
 */
void add_elements( SpVector* sv_lhs, const SpVectorSlice& sv_rhs, size_type offset = 0 );

inline
///
/** Create a dense representation of a sparse vector.
 *
 * The primary use if the function is to create a VectorSlice
 * object that only represents the nonzero values of the
 * sparse vector.  This could have several different uses
 * but one of the most significant examples is when you want
 * to discard the indices when sv_rhs->size() == sv_rhs->nz() and
 * sv_rhs->is_sorted() == true.
 */
VectorSlice dense_view( SpVectorSlice& sv_rhs );

inline
///
const VectorSlice dense_view( const SpVectorSlice& sv_rhs );

} // end namespace SparseLinAlgPack

// /////////////////////////////////////////
// Inline function definitions

inline
LinAlgPack::VectorSlice
SparseLinAlgPack::dense_view( SpVectorSlice& sv_rhs )
{
	return sv_rhs.nz()
		? VectorSlice( &sv_rhs.begin()->value(), sv_rhs.nz(), 2 )
		: VectorSlice( NULL, 0, 0 );
}

inline
const LinAlgPack::VectorSlice
SparseLinAlgPack::dense_view( const SpVectorSlice& sv_rhs )
{
	return sv_rhs.nz()
		? VectorSlice( &const_cast<SpVectorSlice&>(sv_rhs).begin()->value(), sv_rhs.nz(), 2 )
		: VectorSlice( NULL, 0, 0 );
}

#endif // SP_VECTOR_CLASS_H
