/* /////////////////////////////////////////////
// RTOp_TOp_inv_of_difference.h
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
//       on 6/27/2002 at 20:41
//
*/

#ifndef RTOp_TOp_inv_of_difference_H
#define RTOp_TOp_inv_of_difference_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_inv_of_difference.h
 *
 \verbatim


 element-wise transformation : z0 = alpha/(v0-v1);
 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * Calculates a vector that is the inverse of the
 *  element-wise difference of two other vectors
 *
 */
/*@{*/

/** Name of this transformation operator class */
extern const char RTOp_TOp_inv_of_difference_name[];

/** Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_inv_of_difference_vtbl;

/** Constructor */
int RTOp_TOp_inv_of_difference_construct( RTOp_value_type alpha,  struct RTOp_RTOp* op );

/** Destructor */
int RTOp_TOp_inv_of_difference_destroy( struct RTOp_RTOp* op );

/** Initialize the state of the operator object */
int RTOp_TOp_inv_of_difference_init( RTOp_value_type alpha, struct RTOp_RTOp* op );



/*@}*/

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_TOp_inv_of_difference_H */
