/* /////////////////////////////////////////
// RTOp_obj_value_vtbl.c
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
#include <malloc.h>

#include "RTOp_obj_value_vtbl.h"
#include "RTOp_obj_free_free.h"

/* Local function definitions */

static int get_obj_type_num_entries(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void* instance_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
  *num_values  = 1;
  *num_indexes = 0;
  *num_chars   = 0;
  return 0;
}

static int obj_create(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void** obj )
{
  const int mem_size = sizeof(RTOp_value_type);
  *obj = malloc( mem_size );
  *((RTOp_value_type*)*obj) = 0.0;
  return 0;
}

static int obj_reinit(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void* obj )
{
  *((RTOp_value_type*)obj) = 0.0;
  return 0;
}

static int extract_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *       instance_data
  ,void *             obj
  ,int                num_values
  ,RTOp_value_type    value_data[]
  ,int                num_indexes
  ,RTOp_index_type    index_data[]
  ,int                num_chars
  ,RTOp_char_type     char_data[]
  )
{
  assert(obj);
  assert( num_values  == 1 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
  value_data[0] = *((RTOp_value_type*)obj);
  return 0;
}

static int load_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *            instance_data
  ,int                     num_values
  ,const RTOp_value_type   value_data[]
  ,int                     num_indexes
  ,const RTOp_index_type   index_data[]
  ,int                     num_chars
  ,const RTOp_char_type    char_data[]
  ,void **                 obj
  )
{
  assert( obj );
  assert( num_values  == 1 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
  if(*obj == NULL)
    obj_create(vtbl,instance_data,obj);
  *((RTOp_value_type*)*obj) = value_data[0];
  return 0;
}

const struct RTOp_obj_type_vtbl_t   RTOp_obj_value_vtbl =
{
   get_obj_type_num_entries
  ,obj_create
  ,obj_reinit
  ,RTOp_obj_free_free
  ,extract_state
  ,load_state
};
