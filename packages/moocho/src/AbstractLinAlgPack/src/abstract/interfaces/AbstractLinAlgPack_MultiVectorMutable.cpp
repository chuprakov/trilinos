// ///////////////////////////////////////////////////////////////
// MultiVectorMutable.cpp
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

#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "WorkspacePack.h"

namespace AbstractLinAlgPack {

// Sub-view methods

MultiVectorMutable::multi_vec_mut_ptr_t
MultiVectorMutable::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng)
{
	assert(0); // ToDo: return a MultiVectorMutableSubView object.
	// Note that the MultiVectorMutableSubView class should derive from
	// MultiVectorSubView.
	return ReferenceCountingPack::null;
}

// Collective apply_transformation() methods

void MultiVectorMutable::apply_transformation(
	EApplyBy apply_by, const RTOpPack::RTOp& prim_op
	,const size_t num_multi_vecs,         const MultiVector**   multi_vecs
	,const size_t num_targ_multi_vecs_in, MultiVectorMutable**  targ_multi_vecs_in
	,RTOp_ReductTarget reduct_objs[]
	,const index_type prim_first_ele, const index_type prim_sub_dim, const index_type prim_global_offset
	,const index_type sec_first_ele, const index_type sec_sub_dim
	)
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// ToDo: Validate the input!

	wsp::Workspace<MultiVectorMutable*> targ_multi_vecs(wss,num_targ_multi_vecs_in+1);
	targ_multi_vecs[0] = this; // I am the first!
	for(size_type k = 1; k <= num_targ_multi_vecs_in; ++k)
		targ_multi_vecs[k] = targ_multi_vecs_in[k-1];

	this->apply_op(
		apply_by, prim_op
		,num_multi_vecs,            multi_vecs
		,num_targ_multi_vecs_in+1,  &targ_multi_vecs[0]
		,reduct_objs
		,prim_first_ele, prim_sub_dim, prim_global_offset
		,sec_first_ele, sec_sub_dim
		);
}

void MultiVectorMutable::apply_transformation(
	EApplyBy apply_by, const RTOpPack::RTOp& prim_op, const RTOpPack::RTOp& sec_op
	,const size_t num_multi_vecs,         const MultiVector**   multi_vecs
	,const size_t num_targ_multi_vecs_in, MultiVectorMutable**  targ_multi_vecs_in
	,RTOp_ReductTarget reduct_obj
	,const index_type prim_first_ele, const index_type prim_sub_dim, const index_type prim_global_offset
	,const index_type sec_first_ele, const index_type sec_sub_dim
	)
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// ToDo: Validate the input!

	wsp::Workspace<MultiVectorMutable*> targ_multi_vecs(wss,num_targ_multi_vecs_in+1);
	targ_multi_vecs[0] = this; // I am the first!
	for(size_type k = 1; k <= num_targ_multi_vecs_in; ++k)
		targ_multi_vecs[k] = targ_multi_vecs_in[k-1];

	this->apply_op(
		apply_by, prim_op, sec_op
		,num_multi_vecs,            multi_vecs
		,num_targ_multi_vecs_in+1,  &targ_multi_vecs[0]
		,reduct_obj
		,prim_first_ele, prim_sub_dim, prim_global_offset
		,sec_first_ele, sec_sub_dim
		);
}

// Overriddend form MultiVector

MultiVectorMutable::vec_ptr_t MultiVectorMutable::row(index_type i) const
{
	return const_cast<MultiVectorMutable*>(this)->row(i);
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::col(index_type j) const
{
	return const_cast<MultiVectorMutable*>(this)->col(j);
}

MultiVectorMutable::vec_ptr_t MultiVectorMutable::diag(int k) const
{
	return const_cast<MultiVectorMutable*>(this)->diag(k);
}

MultiVectorMutable::multi_vec_ptr_t
MultiVectorMutable::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	return const_cast<MultiVectorMutable*>(this)->mv_sub_view(row_rng,col_rng);
}

} // end namespace AbstractLinAlgPack
