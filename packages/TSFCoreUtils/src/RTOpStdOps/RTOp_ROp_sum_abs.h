/* ///////////////////////////////////////////// */
/* RTOp_ROp_sum_abs.h */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 8/12/2003 at 19:31 */
/* */

#ifndef RTOp_ROp_sum_abs_H
#define RTOp_ROp_sum_abs_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_sum_abs.h
 *
 \verbatim

element-wise reduction:
    abs_sum_ith = abs(v0);


 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * ToDo: Write the documentation for this class!
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_sum_abs_vtbl;

/* Constructor */
int RTOp_ROp_sum_abs_construct(  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_sum_abs_destroy( struct RTOp_RTOp* op );


/* Extract the value of the reduction object */
RTOp_value_type RTOp_ROp_sum_abs_val(RTOp_ReductTarget reduct_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_ROp_sum_abs_H */
