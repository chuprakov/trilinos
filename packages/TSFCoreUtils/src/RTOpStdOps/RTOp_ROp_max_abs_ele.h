/* /////////////////////////////////////////////
// RTOp_ROp_max_abs_ele.h
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
*/

#ifndef RTOP_ROP_MAX_ABS_ELE_H
#define RTOP_ROP_MAX_ABS_ELE_H

#include "RTOp.h"
#include "RTOp_obj_value_index_vtbl.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_max_abs_ele.h   Returns the element with the maximum
 * absolute value and its index.
 *
 * This reduction operator finds <tt>value</tt> and <tt>index</tt> such that
 \verbatim

 value = v0(index) s.t. |v0(index)| >= v0(i), for i = 1...n
 \endverbatim
 */
/*@{*/

/** Name of this reduction operator class */
extern const char RTOp_ROp_max_abs_ele_name[];

/** Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_abs_ele_vtbl;

/** Constructor */
int RTOp_ROp_max_abs_ele_construct( struct RTOp_RTOp* op );

/** Destructor */
int RTOp_ROp_max_abs_ele_destroy( struct RTOp_RTOp* op );

/** Extract the concrete reduction target object from its pointer (handle). */
struct RTOp_value_index_type
RTOp_ROp_max_abs_ele_val(RTOp_ReductTarget targ_obj);

/*@}*/

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_ROP_MAX_ABS_ELE_H */
