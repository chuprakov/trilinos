// //////////////////////////////////////////////////////////////////////////////////
// VectorAssign.h
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
