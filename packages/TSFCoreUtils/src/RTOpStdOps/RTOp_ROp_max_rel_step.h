/* /////////////////////////////////////////////
// RTOp_ROp_max_rel_step.h
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

//
// Note: This file was created automatically by 'new_rtop.pl'
//       on 6/24/2002 at 21:2
//
*/

#ifndef RTOp_ROp_max_rel_step_H
#define RTOp_ROp_max_rel_step_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_max_rel_step.h
 *
 * This reduction operator finds the maximum relative step change in \c v0 for:
 \verbatim

 v0 = v0 + v1;
 \endverbatim
 *
 * This is found by the reduciton:
 \verbatim

 element-wise reduction      : gamma = max{ |v1(i)| / ( 1.0 + |v0| ), i = 1...n }
 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 */
/*@{*/

/** Name of this transformation operator class */
extern const char RTOp_ROp_max_rel_step_name[];

/** Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_rel_step_vtbl;

/** Constructor */
int RTOp_ROp_max_rel_step_construct(  struct RTOp_RTOp* op );

/** Destructor */
int RTOp_ROp_max_rel_step_destroy( struct RTOp_RTOp* op );

/** Extract the value of the reduction object \c gamma */
RTOp_value_type RTOp_ROp_max_rel_step_val(RTOp_ReductTarget reduct_obj);

/*@}*/

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_ROp_max_rel_step_H */
