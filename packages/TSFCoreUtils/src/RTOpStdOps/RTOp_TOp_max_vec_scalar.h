/* /////////////////////////////////////////////
// RTOp_TOp_max_vec_scalar.h

//
// Note: This file was created automatically by 'new_rtop.pl'
//       on 7/13/2002 at 13:44
//
*/

#ifndef RTOp_TOp_max_vec_scalar_H
#define RTOp_TOp_max_vec_scalar_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_max_vec_scalar.h
 *
 \verbatim

element-wise transformation:
    z0 = max(z0,min_ele);
 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * ToDo: Write the documentation for this class!
 */
/*@{*/

/** Name of this transformation operator class */
extern const char RTOp_TOp_max_vec_scalar_name[];

/** Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_max_vec_scalar_vtbl;

/** Constructor */
int RTOp_TOp_max_vec_scalar_construct( RTOp_value_type min_ele,  struct RTOp_RTOp* op );

/** Destructor */
int RTOp_TOp_max_vec_scalar_destroy( struct RTOp_RTOp* op );

/** Initialize the state of the operator object */
int RTOp_TOp_max_vec_scalar_init( RTOp_value_type min_ele, struct RTOp_RTOp* op );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_TOp_max_vec_scalar_H */
