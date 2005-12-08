// //////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_DVectorAssign.hpp
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

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

/* * @name {\bf DVector/DVectorSlice assignment functions}.
  */
// @{

/// v_lhs = alpha (elementwise)
void assign(DVector* v_lhs, value_type alpha);
/// v_lhs = vs_rhs.
void assign(DVector* v_lhs, const DVectorSlice& vs_rhs);
/// vs_lhs = alpha (elementwise)
void assign(DVectorSlice* vs_lhs, value_type alpha);
/// vs_lhs = vs_rhs
void assign(DVectorSlice* vs_lhs, const DVectorSlice& vs_rhs);

// @}


} // end namespace DenseLinAlgPack

#endif	// VECTOR_ASSIGN_H
