// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
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

// /////////////////////////////////////////////////////////////////////////////
// TSFCore_apply_op_helper.hpp

#ifndef TSFCORE_APPLY_OP_HELPER_HPP
#define TSFCORE_APPLY_OP_HELPER_HPP

#include "TSFCore_apply_op_helper_decl.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"

template<class Scalar>
void TSFCore::apply_op_validate_input(
	const char                      func_name[]
	,const RTOpPack::RTOpT<Scalar>  &op
	,const int                   num_vecs
	,const Vector<Scalar>*          vecs[]
	,const int                   num_targ_vecs
	,Vector<Scalar>*                targ_vecs[]
	,RTOpPack::ReductTarget         *reduct_obj
	,const Index                    first_ele_in
	,const Index                    sub_dim_in
	,const Index                    global_offset_in
	)
{
	int k;
	const VectorSpace<Scalar>
		&space = ( num_vecs ? *vecs[0]->space() : *targ_vecs[0]->space() );
	const Index
		dim = space.dim();
	TEST_FOR_EXCEPTION(
		global_offset_in < 0, std::logic_error
		,func_name << " : Error!  global_offset_in = "
		<<global_offset_in<<" is not valid" );
	TEST_FOR_EXCEPTION(
		first_ele_in > dim, std::logic_error
		,func_name << " : Error!  first_ele_in = "
		<<first_ele_in<<" is not compatible with space.dim() = " << dim );
	TEST_FOR_EXCEPTION(
		sub_dim_in < 0 || (sub_dim_in > 0 && sub_dim_in > dim-(first_ele_in-1)), std::logic_error
		,func_name << " : Error!  first_ele_in = "
		<<first_ele_in<<" and sub_dim_in = "<<sub_dim_in
		<<" is not compatible with space.dim() = " << dim );
	for(k = 0; k < num_vecs; ++k) {
		const bool isCompatible = space.isCompatible(*vecs[k]->space());
		TEST_FOR_EXCEPTION(
			!isCompatible, Exceptions::IncompatibleVectorSpaces
			,func_name << " : Error!  vecs["<<k<<"]"
			" with dimension vecs["<<k<<"].dim() = " << vecs[k]->space()->dim()
			<< " is not compatible with space.dim() = " << dim
			);
	}
	for(k = 0; k < num_targ_vecs; ++k) {
		const bool isCompatible = space.isCompatible(*targ_vecs[k]->space());
		TEST_FOR_EXCEPTION(
			!isCompatible, Exceptions::IncompatibleVectorSpaces
			,func_name << " : Error!  targ_vecs["<<k<<"]"
			" with dimension targ_vecs["<<k<<"].dim() = " << targ_vecs[k]->space()->dim()
			<< " is not compatible with space.dim() = " << dim
			);
	}
}

template<class Scalar>
void TSFCore::apply_op_serial(
	const RTOpPack::RTOpT<Scalar>  &op
	,const int                  num_vecs
	,const Vector<Scalar>*         vecs[]
	,const int                  num_targ_vecs
	,Vector<Scalar>*               targ_vecs[]
	,RTOpPack::ReductTarget        *reduct_obj
	,const Index                   first_ele_in
	,const Index                   sub_dim_in
	,const Index                   global_offset_in
	)
{
 	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
	
	// Dimension of global sub-vector
	const VectorSpace<Scalar>
		&space = ( num_vecs ? *vecs[0]->space() : *targ_vecs[0]->space() );
	const Index
		full_dim       = space.dim(),
		global_sub_dim = sub_dim_in ? sub_dim_in : full_dim - (first_ele_in-1);
	const Range1D
		global_sub_rng = Range1D(first_ele_in,(first_ele_in-1)+global_sub_dim);

	//
	// Get explicit views of the vector elements
	//

	Workspace<RTOpPack::SubVectorT<Scalar> >         local_vecs(wss,num_vecs);
	Workspace<RTOpPack::MutableSubVectorT<Scalar> >  local_targ_vecs(wss,num_targ_vecs);
	int k;
	for(k = 0; k < num_vecs; ++k) {
		RTOpPack::SubVectorT<Scalar> &v = local_vecs[k];
		vecs[k]->getSubVector( global_sub_rng, &v );
		v.setGlobalOffset( global_offset_in );
	}
	for(k = 0; k < num_targ_vecs; ++k) {
		RTOpPack::MutableSubVectorT<Scalar> &v = local_targ_vecs[k];
		targ_vecs[k]->getSubVector( global_sub_rng, &v );
		v.setGlobalOffset( global_offset_in );
	}

	//
	// Apply the reduction/transformation operator on all elements all at once!
	//

	op.apply_op(
		num_vecs,       num_vecs      ? &local_vecs[0]      : NULL
		,num_targ_vecs, num_targ_vecs ? &local_targ_vecs[0] : NULL
		,reduct_obj
		);

	//
	// Free (and commit) the explicit views of the vector elements
	// which should also inform the vectors that they have
	// changed.
	//

	for(k = 0; k < num_vecs; ++k) {
		RTOpPack::SubVectorT<Scalar> &v = local_vecs[k];
		v.setGlobalOffset( global_sub_rng.lbound() - 1);
		vecs[k]->freeSubVector(&v);
	}
	for(k = 0; k < num_targ_vecs; ++k) {
		RTOpPack::MutableSubVectorT<Scalar> &v = local_targ_vecs[k];
		v.setGlobalOffset( global_sub_rng.lbound() - 1);
		targ_vecs[k]->commitSubVector(&v);
	}

}

#endif // TSFCORE_APPLY_OP_HELPER_HPP
