// /////////////////////////////////////////////
// RTOp_TOp_axpy.h
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

#ifndef RTOP_TOP_AXPY_H
#define RTOP_TOP_AXPY_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_axpy.h Perform an "axpy" operation.
  *
  * <tt>targ_vec(i) <- alpha * vec[0](i) + targ_vec(i), for i = 1...n</tt>.
  *
  * This operation is defined for dense and sparse sub-vectors
  * <tt>v[0]</tt> (<tt>num_vecs == 1</tt>).
  */
//@{

/// Name of this transformation operator class
extern const char RTOp_TOp_axpy_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_axpy_vtbl;

/// Constructor
int RTOp_TOp_axpy_construct( RTOp_value_type alpha, struct RTOp_RTOp* op );

/// Destructor
int RTOp_TOp_axpy_destroy( struct RTOp_RTOp* op );

/// Reset alpha
int RTOp_TOp_axpy_set_alpha( RTOp_value_type alpha, struct RTOp_RTOp* op );

//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOP_TOP_AXPY_H
