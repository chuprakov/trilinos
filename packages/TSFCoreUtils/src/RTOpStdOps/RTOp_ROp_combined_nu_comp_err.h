// /////////////////////////////////////////////
// RTOp_ROp_combined_nu_comp_err.h
//
// Copyright (C) 2001
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
//       on 6/27/2002 at 15:7
//

#ifndef RTOp_ROp_combined_nu_comp_err_H
#define RTOp_ROp_combined_nu_comp_err_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  // vtbl for reduction object data

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_combined_nu_comp_err.h
 *
 \verbatim

 element-wise reduction      : comp_err = max(comp_err, v0(i)*(v3(i)-v1(i), -v0(i)*(v1(i)-v2(i))));

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * This operator calculates an estimate of the complementarity error using a combined
 *  nu for the upper and lower bound inequality constraints. At the soln, a positive
 *  nu indicates the upper bound is active, while a negative nu indicates the lower
 *  bound is active. If nu is zero, then neither bound is active.
 *
 *
 \verbatim
      for every i..
       comp_err = max(comp_err, v(i)*(xu(i)-x(i), -v(i)*(x(i)-xl(i))));
\endverbatim
 */
//@{

/// Name of this transformation operator class
extern const char RTOp_ROp_combined_nu_comp_err_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_combined_nu_comp_err_vtbl;

/// Constructor
int RTOp_ROp_combined_nu_comp_err_construct(  struct RTOp_RTOp* op );

/// Destructor
int RTOp_ROp_combined_nu_comp_err_destroy( struct RTOp_RTOp* op );


/// Extract the value of the reduction object
RTOp_value_type RTOp_ROp_combined_nu_comp_err_val(RTOp_ReductTarget reduct_obj);

//@}


/**
 *
 \verbatim

 element-wise reduction      : comp_err = max(comp_err, v0(i)*(v1(i)-v2(i)));

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * This operator calculates the comp_err for a single bound vector
 *  (upper or lower) - it is to be used when either one of the bounds
 *  vectors are all infinite (-/+)
 */
//@{

/// Name of this transformation operator class
extern const char RTOp_ROp_combined_nu_comp_err_one_only_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_combined_nu_comp_err_one_only_vtbl;

/// Constructor
int RTOp_ROp_combined_nu_comp_err_one_only_construct(  struct RTOp_RTOp* op );

/// Destructor
int RTOp_ROp_combined_nu_comp_err_one_only_destroy( struct RTOp_RTOp* op );


/// Extract the value of the reduction object
RTOp_value_type RTOp_ROp_combined_nu_comp_err_one_only_val(RTOp_ReductTarget reduct_obj);

//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOp_ROp_combined_nu_comp_err_H
