// //////////////////////////////////////////////////////////////////////////////////
// VectorAssign.h

#ifndef VECTOR_ASSIGN_H
#define VECTOR_ASSIGN_H

#include "LinAlgPackTypes.h"

namespace LinAlgPack {

/** @name {\bf Vector/VectorSlice assignment functions}.
  */
//@{

/// v_lhs = alpha (elementwise)
void assign(Vector* v_lhs, value_type alpha);
/// v_lhs = vs_rhs.
void assign(Vector* v_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = alpha (elementwise)
void assign(VectorSlice* vs_lhs, value_type alpha);
/// vs_lhs = vs_rhs
void assign(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);

//@}


} // end namespace LinAlgPack

#endif	// VECTOR_ASSIGN_H
