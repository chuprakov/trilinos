// /////////////////////////////////////////////
// RTOp_TOp_ele_wise_sqrt.c

//
// Note: This file was created automatically by 'new_rtop.pl'
//       on 7/2/2002 at 19:1
//

#include <assert.h>
#include <math.h>

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

#include "RTOp_TOp_ele_wise_sqrt.h"
#include "RTOp_obj_null_vtbl.h"  // vtbl for operator object instance data



// Implementation functions for RTOp_RTOp

static int RTOp_TOp_ele_wise_sqrt_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget reduct_obj )
{
  //
  // Declare local variables
  //

  // Vector data
  RTOp_index_type           sub_dim;
  // z0
  RTOp_value_type           *z0_val;
  ptrdiff_t                 z0_val_s;

  // Automatic temporary variables
  register RTOp_index_type  k;

  //
  // Validate the input
  //
  if( num_vecs != 0 || ( num_vecs && vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 || ( num_targ_vecs && targ_vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;


  //
  // Get pointers to data
  //
  sub_dim       = targ_vecs[0].sub_dim;
  // z0
  z0_val        = targ_vecs[0].values;
  z0_val_s      = targ_vecs[0].values_stride;


  //
  // Apply the operator:
  //
  for( k = 0; k < sub_dim; ++k, z0_val += z0_val_s )
    {
    // Element-wise transformation
    assert(*z0_val >= 0);
    (*z0_val) = sqrt((*z0_val));
    }

  return 0; // success?
}



// Name of this transformation operator class
const char RTOp_TOp_ele_wise_sqrt_name[] = "TOp_ele_wise_sqrt";

// Virtual function table
const struct RTOp_RTOp_vtbl_t RTOp_TOp_ele_wise_sqrt_vtbl =
{
  &RTOp_obj_null_vtbl
  ,&RTOp_obj_null_vtbl
  ,NULL
  ,RTOp_TOp_ele_wise_sqrt_apply_op
  ,NULL
  ,NULL
};

// Class specific functions

int RTOp_TOp_ele_wise_sqrt_construct(  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(op);
#endif
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_TOp_ele_wise_sqrt_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return 0;
}

int RTOp_TOp_ele_wise_sqrt_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0;
}



