// /////////////////////////////////////////////
// RTOp_TOp_assign_vectors.h
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

#ifndef RTOP_TOP_ASSIGN_VECTORS_H
#define RTOP_TOP_ASSIGN_VECTORS_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_assign_vectors.h Assign the elements of one vector to another.
  *
  * <tt>targ_vec(i) <- vec[0](i), for i = 1...n</tt>
  *
  * This operation is only defined to allow vector argument
  * (<tt>num_vecs == 1</tt>) <tt>v[0]</tt>, but can handle a dense or
  * sparse vector.
  */
//@{

/// Name of this transformation operator class
extern const char RTOp_TOp_assign_vectors_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_assign_vectors_vtbl;

/// Constructor
int RTOp_TOp_assign_vectors_construct( struct RTOp_RTOp* op );

/// Destructor
int RTOp_TOp_assign_vectors_destroy( struct RTOp_RTOp* op );

//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOP_TOP_ASSIGN_VECTORS_H
