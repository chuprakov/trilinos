/* /////////////////////////////////////////////
// RTOp_ROp_comp_err_with_mu.h

//
// Note: This file was created automatically by 'new_rtop.pl'
//       on 7/24/2002 at 23:46
//
*/

#ifndef RTOp_ROp_comp_err_with_mu_H
#define RTOp_ROp_comp_err_with_mu_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_comp_err_with_mu.h
 *
 \verbatim

element-wise reduction:
    comp_err_ith = max(comp_err_ith, fabs(v3*(v0-v1)-mu));
    comp_err_ith = max(comp_err_ith, fabs(v4*(v2-v1)-mu));


 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * ToDo: Write the documentation for this class!
 */
/*@{*/

/** Name of this transformation operator class */
extern const char RTOp_ROp_comp_err_with_mu_name[];

/** Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_comp_err_with_mu_vtbl;

/** Constructor */
int RTOp_ROp_comp_err_with_mu_construct( RTOp_value_type mu, RTOp_value_type inf_bound,  struct RTOp_RTOp* op );

/** Destructor */
int RTOp_ROp_comp_err_with_mu_destroy( struct RTOp_RTOp* op );

/** Initialize the state of the operator object */
int RTOp_ROp_comp_err_with_mu_init( RTOp_value_type mu, RTOp_value_type inf_bound, struct RTOp_RTOp* op );

/** Extract the value of the reduction object */
RTOp_value_type RTOp_ROp_comp_err_with_mu_val(RTOp_ReductTarget reduct_obj);

/*@}*/

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_ROp_comp_err_with_mu_H */
