// /////////////////////////////////////////////////////
// RTOp_SparseSubVector.c
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
//

#include <assert.h>

#include "RTOp_SparseSubVector.h"

void RTOp_sparse_sub_vector(
	RTOp_index_type global_offset, RTOp_index_type sub_dim
	,RTOp_index_type sub_nz, const RTOp_value_type values[], ptrdiff_t values_stride
	,const RTOp_index_type indices[], ptrdiff_t indices_stride
	,ptrdiff_t local_offset, int is_sorted
	,struct RTOp_SparseSubVector *sub_vec
	)
{
	// Validate input
#ifdef RTOp_DEBUG
	assert( sub_vec );
	assert(
		( sub_nz && ( values != NULL && indices != NULL && indices_stride != 0 && sub_nz <= sub_dim ) )
		|| !sub_nz || ( sub_nz == sub_dim && indices == NULL )
		);
#endif
	// Set members
	sub_vec->global_offset  = global_offset;
	sub_vec->sub_dim        = sub_dim;
	sub_vec->sub_nz         = sub_nz;
	sub_vec->values         = values;
	sub_vec->values_stride  = values_stride;
	sub_vec->indices        = indices;
	sub_vec->indices_stride = indices_stride;
	sub_vec->local_offset   = local_offset;
	sub_vec->is_sorted      = is_sorted;
}

void RTOp_sparse_sub_vector_null( struct RTOp_SparseSubVector *sub_vec )
{
	sub_vec->global_offset  = 0;
	sub_vec->sub_dim        = 0;
	sub_vec->sub_nz         = 0;
	sub_vec->values         = NULL;
	sub_vec->values_stride  = 0;
	sub_vec->indices        = NULL;
	sub_vec->indices_stride = 0;
	sub_vec->local_offset   = 0;
	sub_vec->is_sorted      = 0;
}

void RTOp_sparse_sub_vector_from_dense(
	const struct RTOp_SubVector     *sub_vec
	,struct RTOp_SparseSubVector    *spc_sub_vec
	)
{
	spc_sub_vec->global_offset  = sub_vec->global_offset;
	spc_sub_vec->sub_dim        = sub_vec->sub_dim;
	spc_sub_vec->sub_nz         = sub_vec->sub_dim;
	spc_sub_vec->values         = sub_vec->values;
	spc_sub_vec->values_stride  = sub_vec->values_stride;
	spc_sub_vec->indices        = NULL;
	spc_sub_vec->indices_stride = 0;
	spc_sub_vec->local_offset   = 0;
	spc_sub_vec->is_sorted      = 0;
}
