/* ///////////////////////////////////////////////////////////////
// RTOp_reduct_max_value.c
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

#define MY_MAX(a,b) a > b ? a : b

#include "RTOp_reduct_max_value.h"

int RTOp_reduct_max_value(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget in_targ_obj, RTOp_ReductTarget inout_targ_obj )
{
  /* inout_dot_prod += in_dot_prod */
  *((RTOp_value_type*)inout_targ_obj)
    = MY_MAX( *((RTOp_value_type*)inout_targ_obj)
          ,*((RTOp_value_type*)in_targ_obj)
      );
  return 0;
}

static void external_reduct_op( void* in_targ_array, void* inout_targ_array
  , int* len, RTOp_Datatype* datatype )
{
  /* inout_dot_prod += in_dot_prod */
  RTOp_value_type /* index past the size members */
    *in_targs    = (RTOp_value_type*)in_targ_array    + 3,
    *inout_targs = (RTOp_value_type*)inout_targ_array + 3;
  int i;
  for( i = 0; i < *len; ++i, inout_targs += 4, in_targs += 4 )
    *inout_targs = MY_MAX(*inout_targs,*in_targs);
}

int RTOp_get_reduct_max_value_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr )
{
  *reduct_op_func_ptr = external_reduct_op;
  return 0;
}
