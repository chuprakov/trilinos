/*
// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include <assert.h>

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
