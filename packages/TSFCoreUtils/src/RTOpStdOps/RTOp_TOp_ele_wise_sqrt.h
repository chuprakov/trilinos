/* ///////////////////////////////////////////// */
/* RTOp_TOp_ele_wise_sqrt.h */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/2/2002 at 19:1 */
/* */

#ifndef RTOp_TOp_ele_wise_sqrt_H
#define RTOp_TOp_ele_wise_sqrt_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_ele_wise_sqrt.h
 *
 \verbatim


element-wise transformation:
    z0 = sqrt(z0);

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * ToDo: Write the documentation for this class!
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_ele_wise_sqrt_vtbl;

/* Constructor */
int RTOp_TOp_ele_wise_sqrt_construct(  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_ele_wise_sqrt_destroy( struct RTOp_RTOp* op );




/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_TOp_ele_wise_sqrt_H */
