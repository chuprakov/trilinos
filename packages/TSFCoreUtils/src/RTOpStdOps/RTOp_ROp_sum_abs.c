/* ///////////////////////////////////////////// */
/* RTOp_ROp_sum_abs.c */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 8/12/2003 at 19:31 */
/* */

#include <assert.h>
#include <math.h>

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

#include "RTOp_ROp_sum_abs.h"
#include "RTOp_obj_null_vtbl.h"     /* vtbl for operator object instance data */
#include "RTOp_reduct_sum_value.h"  /* Reduction of intermediate reduction objects */

/* Implementation functions for RTOp_RTOp */

static int RTOp_ROp_sum_abs_apply_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, const int num_vecs, const struct RTOp_SubVector vecs[]
	, const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
	, RTOp_ReductTarget reduct_obj )
{
	/* */
	/* Declare local variables */
	/* */

	/* Access to the reduction object data */
	RTOp_value_type *abs_sum = (RTOp_value_type*)reduct_obj;
	/* Vector data */
	RTOp_index_type           sub_dim;
	/* v0 */
	const RTOp_value_type     *v0_val;
	ptrdiff_t                 v0_val_s;

	/* Automatic temporary variables */
	register RTOp_index_type  k;
	/* Temporary element-wise reduction object */
	RTOp_value_type abs_sum_ith;

	/* */
	/* Validate the input */
	/* */
	if( num_vecs != 1 || ( num_vecs && vecs == NULL ) )
		return RTOp_ERR_INVALID_NUM_VECS;
	if( num_targ_vecs != 0 || ( num_targ_vecs && targ_vecs == NULL ) )
		return RTOp_ERR_INVALID_NUM_TARG_VECS;
	assert(reduct_obj);


	/* */
	/* Get pointers to data */
	/* */
	sub_dim       = vecs[0].sub_dim;
	/* v0 */
	v0_val        = vecs[0].values;
	v0_val_s      = vecs[0].values_stride;


	/* */
	/* Apply the operator: */
	/* */
	for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s )
	{
		/* Element-wise reduction */
		abs_sum_ith = abs((*v0_val));
		/* Reduction of intermediates */
		(*abs_sum) += abs_sum_ith;
	}

	return 0; /* success? */
}




/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_ROp_sum_abs_vtbl =
{
	&RTOp_obj_null_vtbl
	,&RTOp_obj_value_vtbl
	,"ROp_sum_abs"
	,NULL
	,RTOp_ROp_sum_abs_apply_op
	,RTOp_reduct_sum_value
	,RTOp_get_reduct_sum_value_op
};

/* Class specific functions */

int RTOp_ROp_sum_abs_construct(  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
	assert(op);
#endif
	op->obj_data  = NULL;
	op->vtbl      = &RTOp_ROp_sum_abs_vtbl;
	op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
	return 0;
}

int RTOp_ROp_sum_abs_destroy( struct RTOp_RTOp* op )
{
	op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
	op->obj_data  = NULL;
	op->vtbl      = NULL;
	return 0;
}


RTOp_value_type RTOp_ROp_sum_abs_val(RTOp_ReductTarget reduct_obj)
{
	return *((RTOp_value_type*)reduct_obj);
}

