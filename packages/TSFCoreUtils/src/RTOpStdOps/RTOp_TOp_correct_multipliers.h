// /////////////////////////////////////////////
// RTOp_TOp_Correct_Multipliers.h

//
// Note: This file was created automatically by 'new_rtop.pl'
//       on 7/1/2002 at 18:11
//

#ifndef RTOp_TOp_correct_multipliers_H
#define RTOp_TOp_correct_multipliers_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  // vtbl for reduction object data

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_Correct_Multipliers.h
 *
 \verbatim


element-wise transformation:
    if (lower_or_upper == 0) // lower bound
        { z0 = (v0 <= inf_bound_limit) ? 0.0 : z0; }
    else // upper bound
        { z0 = (v0 >= inf_bound_limit) ? 0.0 : z0; }

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * This class sets the corresponding multiplier value
 *  to zero if the bound is equal or outside
 *  inf_bound_limit
 *  lower_or_upper : pass 0 for lower bound check
 *                   pass 1 for upper bound check
 */
//@{

/// Name of this transformation operator class
extern const char RTOp_TOp_Correct_Multipliers_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_Correct_Multipliers_vtbl;

/// Constructor
int RTOp_TOp_Correct_Multipliers_construct( RTOp_value_type inf_bound_limit, RTOp_index_type lower_or_upper,  struct RTOp_RTOp* op );

/// Destructor
int RTOp_TOp_Correct_Multipliers_destroy( struct RTOp_RTOp* op );

/// Initialize the state of the operator object
int RTOp_TOp_Correct_Multipliers_init( RTOp_value_type inf_bound_limit, RTOp_index_type lower_or_upper, struct RTOp_RTOp* op );



//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOp_TOp_Correct_Multipliers_H
