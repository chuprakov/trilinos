// /////////////////////////////////////////
// RTOp_obj_value_value_vtbl.c
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

#include <assert.h>
#include <malloc.h>

#include "RTOp_obj_value_value_vtbl.h"
#include "RTOp_obj_free_free.h"

// Local function definitions

static int get_obj_type_num_entries(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void* instance_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
#ifdef RTOp_DEBUG
  assert(num_values);
  assert(num_indexes);
  assert(num_chars);
#endif
  *num_values  = 2;
  *num_indexes = 0;
  *num_chars   = 0;
  return 0;
}

static int obj_create( const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void** obj )
{
  struct RTOp_value_value_type* vi_obj;
  const int mem_size = sizeof(struct RTOp_value_value_type);
#ifdef RTOp_DEBUG
  assert(obj);
#endif
  *obj = malloc( mem_size );
  vi_obj = (struct RTOp_value_value_type*)*obj;
  vi_obj->value1 = 0.0;
  vi_obj->value2 = 0.0;
  return 0;
}

static int obj_reinit( const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void* obj )
{
  struct RTOp_value_value_type* vi_obj;
#ifdef RTOp_DEBUG
  assert(obj);
#endif
  vi_obj = (struct RTOp_value_value_type*)obj;
  vi_obj->value1 = 0.0;
  vi_obj->value2 = 0.0;
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
  struct RTOp_value_value_type* vi_obj;
#ifdef RTOp_DEBUG
  assert( obj );
  assert( num_values  == 2 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
#endif
  vi_obj = (struct RTOp_value_value_type*)obj;
  value_data[0] = vi_obj->value1;
  value_data[1] = vi_obj->value2;
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
  struct RTOp_value_value_type* vi_obj;
#ifdef RTOp_DEBUG
  assert( obj );
  assert( num_values  == 2 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
#endif
  if(*obj == NULL)
    obj_create(vtbl,instance_data,obj);
  vi_obj = (struct RTOp_value_value_type*)*obj;
  vi_obj->value1 = value_data[0];
  vi_obj->value2 = value_data[1];
  return 0;
}

const struct RTOp_obj_type_vtbl_t   RTOp_obj_value_value_vtbl =
{
   get_obj_type_num_entries
  ,obj_create
  ,obj_reinit
  ,RTOp_obj_free_free
  ,extract_state
  ,load_state
};
