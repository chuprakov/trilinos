/* ///////////////////////////////////////////// */
/* RTOp_TOp_multiplier_step.h */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/10/2002 at 1:19 */
/* */

#ifndef RTOp_TOp_multiplier_step_H
#define RTOp_TOp_multiplier_step_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_multiplier_step.h
 *
 \verbatim


element-wise transformation:
    z0 = -v1 + mu*v0 + alpha*v0*v1*v2;

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * This operator calculates the multiplier steps for
 *  an interior point algorithm
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_multiplier_step_vtbl;

/* Constructor */
int RTOp_TOp_multiplier_step_construct( RTOp_value_type mu, RTOp_value_type alpha,  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_multiplier_step_destroy( struct RTOp_RTOp* op );

/* Initialize the state of the operator object */
int RTOp_TOp_multiplier_step_init( RTOp_value_type mu, RTOp_value_type alpha, struct RTOp_RTOp* op );



/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_TOp_multiplier_step_H */
