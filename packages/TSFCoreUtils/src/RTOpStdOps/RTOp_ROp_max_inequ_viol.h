// /////////////////////////////////////////////
// RTOp_ROp_max_inequ_viol.h
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

#ifndef RTOP_ROP_MAX_INEQU_VIOL_H
#define RTOP_ROP_MAX_INEQU_VIOL_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_max_inequ_viol.h Computes the maximum violation of a set of
 * inequality constraints.
 *
 * Determines the maximum violation from a set of double-sided inequality
 * constriants of the form:
 \verbatim

 vL(i) <= v(i) <= vU(i), for i = 1...n
 \endverbatim
 * Constraint violations are scaled by <tt>1.0/(1.0+|v(i)|)</tt> in the
 * determination of the maximum violation.  Any ties are broken by
 * returning the lowest <tt>i</tt>.
 *
 * The input vectors are passed in the order:
 \verbatim

 v0 = v, v1 = vL, v2 = vU
 \endverbatim
 */
//@{

///
/** Reduction target object for this max_inequ_viol operation.
 */
struct RTOp_ROp_max_inequ_viol_reduct_obj_t {
  RTOp_value_type   max_viol;   ///< <tt>max_viol = |v_i-vLU|/(1.0+|v_i|) > 0</tt>
  RTOp_value_type   v_i;        ///< <tt>v_i = v(max_viol_i)</tt>
  RTOp_value_type   vLU_i;      ///< <tt>vLU_i = vL(i)</tt> or <tt>vU(i)</tt>
  RTOp_index_type   max_viol_i; ///< <tt>max_viol_i > 0</tt> if a inequality is violated
  RTOp_index_type   bnd_type;   ///< -1 : LOWER, 0 : EQUALITY, +1 : UPPER
};

/// Name of this reduction operator class
extern const char RTOp_ROp_max_inequ_viol_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_inequ_viol_vtbl;

/// Constructor
int RTOp_ROp_max_inequ_viol_construct( struct RTOp_RTOp* op );

/// Destructor
int RTOp_ROp_max_inequ_viol_destroy( struct RTOp_RTOp* op );

/// Extract the concrete reduction target object from its pointer (handle).
struct RTOp_ROp_max_inequ_viol_reduct_obj_t
RTOp_ROp_max_inequ_viol_val(RTOp_ReductTarget targ_obj);

//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOP_ROP_MAX_INEQU_VIOL_H
