// /////////////////////////////////////////////
// RTOp_ROp_find_nan_inf.h
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

#ifndef RTOP_ROP_FIND_NAN_INF_H
#define RTOP_ROP_FIND_NAN_INF_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_find_nan_inf.h Reduction operator that looks for the first element that is
 * NaN or Inf and returns this index!
 *
 * <tt>targ_obj <- { (v0_i,i) | RTOp_in_nan_inf(v[0](i)) }</tt>
 *
 * This operator is defined to allow exactly one vecto arguments
 * (<tt>num_vecs == 1</tt>) <tt>v[0]</tt> but can handle sparse or dense vectors.
 * The element with the lowest index is selected so that the
 * reduction object returned will be unique for a given vector.
 */
//@{

///
struct RTOp_ROp_find_nan_inf_reduct_obj_t {
  RTOp_value_type v0_i;
  RTOp_index_type i;
};

/// Name of this reduction operator class
extern const char RTOp_ROp_find_nan_inf_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_find_nan_inf_vtbl;

/// Constructor
int RTOp_ROp_find_nan_inf_construct( struct RTOp_RTOp* op );

/// Destructor
int RTOp_ROp_find_nan_inf_destroy( struct RTOp_RTOp* op );

///
/** Extract the number offending element.
 *
 * If <tt>return.i == 0</tt> then no element was found to be NaN or Inf.
 */
struct RTOp_ROp_find_nan_inf_reduct_obj_t
RTOp_ROp_find_nan_inf_val(RTOp_ReductTarget targ_obj);

//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOP_ROP_FIND_NAN_INF_H
