/* ///////////////////////////////////////////////////////////////
// RTOp_reduct_sum_values.c
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
*/

#include <assert.h>

#include "RTOp_reduct_sum_values.h"

int RTOp_reduct_sum_values(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget in_targ_obj, RTOp_ReductTarget inout_targ_obj )
{
  int num_values, k;
#ifdef RTOp_DEBUG
  assert(obj_data);
#endif
  num_values = *(RTOp_index_type*)obj_data;
  /* inout_dot_prod += in_dot_prod */
  for( k = 0; k < num_values; ++k )
    ((RTOp_value_type*)inout_targ_obj)[k] += ((RTOp_value_type*)in_targ_obj)[k];
  return 0;
}

static void external_reduct_op( void* in_targ_array, void* inout_targ_array
  , int* len, RTOp_Datatype* datatype )
{
  /* inout_dot_prod += in_dot_prod */
  RTOp_index_type
    num_values   = *(RTOp_value_type*)in_targ_array; /* num_values is first size member */
  RTOp_value_type /* index past the size members */
    *in_targs    = (RTOp_value_type*)in_targ_array    + 3,
    *inout_targs = (RTOp_value_type*)inout_targ_array + 3;
  int i, k;
  for( i = 0; i < *len; ++i, inout_targs += (3 + num_values), in_targs += (3 + num_values) ) {
    for( k = 0; k < num_values; ++k )
      inout_targs[k] += in_targs[k];
  }
}

int RTOp_get_reduct_sum_values_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr )
{
  *reduct_op_func_ptr = external_reduct_op;
  return 0;
}
