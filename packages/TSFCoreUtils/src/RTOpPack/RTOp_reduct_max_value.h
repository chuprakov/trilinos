// ///////////////////////////////////////////////////////////////
// RTOp_reduct_max_value.h
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

#ifndef RTOP_REDUCT_MAX_VALUE_H
#define RTOP_REDUCT_MAX_VALUE_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @name Definitions of reduction functions for max(.,.) scalar values.
 *
 * These functions perform a simple <tt>max()</tt> of a scalar reduction objects
 * as defined by the virtual function table \Ref{RTOp_obj_value_vtbl}.
 */
//@{

///
/** Use this function for <tt>reduce_reduct_objs</tt> in the RTOp_RTOp_vtbl_t virtual
 * function table.
 */
int RTOp_reduct_max_value(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, RTOp_ReductTarget in_targ_obj, RTOp_ReductTarget inout_targ_obj );

///
/** Use this function for <tt>get_reduct_op</tt> in the RTOp_RTOp_vtbl_t virtual
 * function table.
 */
int RTOp_get_reduct_max_value_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr );

//@}

#ifdef __cplusplus
}
#endif

#endif // RTOP_REDUCT_MAX_VALUE_H
