/* /////////////////////////////////////////////
// RTOp_ROp_norms.h
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

#ifndef RTOP_ROP_NORMS_H
#define RTOP_ROP_NORMS_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_norms.h Reduction operator classes for common norms.
  */
/*@{*/

/** @name One norm reduction operator class.
 *
 * <tt>||v[0]||_1 -> targ_obj</tt>
 */
/*@{*/

/** Name of this reduction operator class */
extern const char RTOp_ROp_norm_1_name[];

/** Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_1_vtbl;

/** Constructor */
int RTOp_ROp_norm_1_construct( struct RTOp_RTOp* op );

/** Extract the value of the norm */
RTOp_value_type RTOp_ROp_norm_1_val(RTOp_ReductTarget targ_obj);

/*@}*/

/** @name Two (Euclidean) norm reduction operator class.
 *
 * <tt>||v[0]||_2 -> targ_obj</tt>
 */
/*@{*/

/** Name of this reduction operator class */
extern const char RTOp_ROp_norm_2_name[];

/** Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_2_vtbl;

/** Constructor */
int RTOp_ROp_norm_2_construct( struct RTOp_RTOp* op );

/** Extract the value of the norm */
RTOp_value_type RTOp_ROp_norm_2_val(RTOp_ReductTarget targ_obj);

/*@}*/

/** @name Infinity norm reduction operator class.
 *
 * <tt>||v[0]||_inf -> targ_obj</tt>
 */
/*@{*/

/** Name of this reduction operator class */
extern const char RTOp_ROp_norm_inf_name[];

/** Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_inf_vtbl;

/** Constructor */
int RTOp_ROp_norm_inf_construct( struct RTOp_RTOp* op );

/**  Extract the value of the norm */
RTOp_value_type RTOp_ROp_norm_inf_val(RTOp_ReductTarget targ_obj);

/*@}*/

/** Destructor (for all three norms) */
int RTOp_ROp_norm_destroy( struct RTOp_RTOp* op );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_ROP_NORMS_H */
