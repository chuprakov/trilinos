// //////////////////////////////////////////////
// RTOpCppToMPI.cpp
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

#include "RTOpCppToMPI.hpp"
#include "RTOpCppToC.hpp"
#include "RTOpToMPI.h"
#include "WorkspacePack.hpp"
#include "Teuchos_TestForException.hpp"

void RTOpPack::MPI_apply_op(
	MPI_Comm comm, const RTOp& op, int root_rank
	,const int num_vecs, const RTOpPack::SubVector sub_vecs[]
	,const int num_targ_vecs, const RTOpPack::MutableSubVector targ_sub_vecs[]
	,RTOp_ReductTarget reduct_obj
	)
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	int k;

	wsp::Workspace<RTOp_SubVector> c_sub_vecs(wss,num_vecs);
	if(sub_vecs) {
		for( k = 0; k < num_vecs; ++k ) {
			const SubVector& v = sub_vecs[k];
			RTOp_sub_vector(v.globalOffset(),v.subDim(),v.values(),v.stride(),&c_sub_vecs[k]);
		}
	}

	wsp::Workspace<RTOp_MutableSubVector> c_targ_sub_vecs(wss,num_targ_vecs);
	if(targ_sub_vecs) {
		for( k = 0; k < num_targ_vecs; ++k ) {
			const MutableSubVector& v = targ_sub_vecs[k];
			RTOp_mutable_sub_vector(v.globalOffset(),v.subDim(),v.values(),v.stride(),&c_targ_sub_vecs[k]);
		}
	}

	RTOpCppToC op_c(op);
	RTOp_ReductTarget reduct_objs[] = { reduct_obj }; 
	RTOp_MPI_apply_op(
		comm,&op_c.op_c(),root_rank,1
		,num_vecs, num_vecs && sub_vecs ? &c_sub_vecs[0] : NULL
		,num_targ_vecs, num_targ_vecs && targ_sub_vecs ? &c_targ_sub_vecs[0] : NULL
		,reduct_obj != RTOp_REDUCT_OBJ_NULL ? reduct_objs : NULL
		);
}

void RTOpPack::MPI_apply_op(
	MPI_Comm comm, const RTOp& op, int root_rank
	,const int num_cols
	,const int num_multi_vecs, const RTOpPack::SubMultiVector sub_multi_vecs[]
	,const int num_targ_multi_vecs, const RTOpPack::MutableSubMultiVector targ_sub_multi_vecs[]
	,RTOp_ReductTarget reduct_objs[]
	)
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	int k, j, off;

	wsp::Workspace<RTOp_SubVector> c_sub_vecs(wss,num_multi_vecs*num_cols);
	if(sub_multi_vecs) {
		for( off = 0, j = 0; j < num_cols; ++j ) {
			for( k = 0; k < num_multi_vecs; ++k ) {
				const SubMultiVector& mv = sub_multi_vecs[k];
				RTOp_sub_vector(mv.globalOffset(),mv.subDim(),&mv(1,j+1),1,&c_sub_vecs[off++]);
			}
		}
	}

	wsp::Workspace<RTOp_MutableSubVector> c_targ_sub_vecs(wss,num_targ_multi_vecs*num_cols);
	if(targ_sub_multi_vecs) {
		for( off = 0, j = 0; j < num_cols; ++j ) {
			for( k = 0; k < num_targ_multi_vecs; ++k ) {
				const MutableSubMultiVector& mv = targ_sub_multi_vecs[k];
				RTOp_mutable_sub_vector(mv.globalOffset(),mv.subDim(),&mv(1,j+1),1,&c_targ_sub_vecs[off++]);
			}
		}
	}

	RTOpCppToC op_c(op);
	RTOp_MPI_apply_op(
		comm,&op_c.op_c(),root_rank,num_cols
		,num_multi_vecs, num_multi_vecs && sub_multi_vecs ? &c_sub_vecs[0] : NULL
		,num_targ_multi_vecs, num_targ_multi_vecs && targ_sub_multi_vecs ? &c_targ_sub_vecs[0] : NULL
		,reduct_objs
		);

}
