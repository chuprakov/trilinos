// /////////////////////////////////////////////
// RTOp_ROp_fraction_to_boundary.h

//
// Note: This file was created automatically by 'new_rtop.pl'
//       on 7/8/2002 at 18:22
//

#ifndef RTOp_ROp_fraction_to_boundary_H
#define RTOp_ROp_fraction_to_boundary_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  // vtbl for reduction object data

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_fraction_to_boundary.h
 *
 \verbatim

element-wise reduction:
    if (v1 < 0)
       { alpha_max_ith = tau*(v0-v2)/(-v1); }
    else
       { alpha_max_ith = tau*(v3-v0)/v1; }


 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * This operator returns alpha_max as calculated by the
 *  fraction to boundary rule
 */
//@{

/// Name of this transformation operator class
extern const char RTOp_ROp_fraction_to_boundary_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_fraction_to_boundary_vtbl;

/// Constructor
int RTOp_ROp_fraction_to_boundary_construct( RTOp_value_type tau,  struct RTOp_RTOp* op );

/// Destructor
int RTOp_ROp_fraction_to_boundary_destroy( struct RTOp_RTOp* op );

/// Initialize the state of the operator object
int RTOp_ROp_fraction_to_boundary_init( RTOp_value_type tau, struct RTOp_RTOp* op );

/// Extract the value of the reduction object
RTOp_value_type RTOp_ROp_fraction_to_boundary_val(RTOp_ReductTarget reduct_obj);

//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOp_ROp_fraction_to_boundary_H
