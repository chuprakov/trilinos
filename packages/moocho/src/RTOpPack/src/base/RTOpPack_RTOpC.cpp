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

// ///////////////////////////////
// RTOpPack_RTOpC.cpp

#include "RTOpPack_RTOpC.hpp"

namespace {

} // namespace

namespace RTOpPack {

RTOpC::RTOpC()
  :RTOpT<RTOp_value_type>("RTOpC") // Should be unused since op_name() if overridden here!
{
	op_.vtbl     = NULL;
	op_.obj_data = NULL;
}

RTOpC::~RTOpC()
{
	if(op_.obj_data)
		RTOp_free_op( &op_ );
}

// Overridden from RTOpT

void RTOpC::get_reduct_type_num_entries(
  int*   num_values
  ,int*  num_indexes
  ,int*  num_chars
  ) const
{
	TEST_FOR_EXCEPTION(
    0!=RTOp_get_reduct_type_num_entries(&op_,num_values,num_indexes,num_chars)
		,UnknownError
		,"RTOpC::get_reduct_type_num_entries(...): Error, "
		"RTOp_get_reduct_type_num_entries(...) returned != 0"
    );
}

Teuchos::RefCountPtr<ReductTarget>
RTOpC::reduct_obj_create() const
{
  RTOp_ReductTarget reduct_obj_raw = RTOp_REDUCT_OBJ_NULL;
	TEST_FOR_EXCEPTION(
    0!=RTOp_reduct_obj_create(&op_,&reduct_obj_raw)
		,UnknownError
		,"RTOpC::reduct_obj_create(...): Error, "
		"RTOp_reduct_obj_create(...) returned != 0"
    );
  return Teuchos::rcp(new ReductTargetC(op_,reduct_obj_raw));
}

void RTOpC::reduce_reduct_objs(
  const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
  ) const
{
	TEST_FOR_EXCEPTION(
    0!=RTOp_reduce_reduct_objs( &op_, (*this)(in_reduct_obj), (*this)(*inout_reduct_obj) )
		,UnknownError
		,"RTOpC::reduce_reduct_objs(...): Error, "
		"RTOp_reduce_reduct_objs(...) returned != 0"
    );
}

void RTOpC::reduct_obj_reinit( ReductTarget* reduct_obj ) const
{
	TEST_FOR_EXCEPTION(
    0!=RTOp_reduct_obj_reinit( &op_, (*this)(*reduct_obj) )
		,UnknownError
		,"RTOpC::reduct_obj_reinit(...): Error, "
		"RTOp_reduct_obj_reinit(...) returned != 0"
    );
}

void RTOpC::extract_reduct_obj_state(
  const ReductTarget     &reduct_obj
  ,int                      num_values
  ,primitive_value_type     value_data[]
  ,int                      num_indexes
  ,index_type               index_data[]
  ,int                      num_chars
  ,char_type                char_data[]
  ) const
{
	TEST_FOR_EXCEPTION(
    0!=RTOp_extract_reduct_obj_state(
      &op_, (*this)(reduct_obj)
      ,num_values,  value_data
      ,num_indexes, index_data
      ,num_chars,   char_data
      )
    ,UnknownError
		,"RTOpC::extract_reduct_obj_state(...): Error, "
		"RTOp_extract_reduct_obj_state(...) returned != 0"
    );
}

void RTOpC::load_reduct_obj_state(
  int                            num_values
  ,const primitive_value_type    value_data[]
  ,int                           num_indexes
  ,const index_type              index_data[]
  ,int                           num_chars
  ,const char_type               char_data[]
  ,ReductTarget               *reduct_obj
  ) const
{
	TEST_FOR_EXCEPTION(
    0!=RTOp_load_reduct_obj_state(
      &op_
      ,num_values,  value_data
      ,num_indexes, index_data
      ,num_chars,   char_data
      ,(*this)(*reduct_obj)
      )
    ,UnknownError
		,"RTOpC::load_reduct_obj_state(...): Error, "
		"RTOp_load_reduct_obj_state(...) returned != 0"
    );
}

void RTOpC::get_op_type_num_entries(
  int*  num_values
  ,int* num_indexes
  ,int* num_chars
  ) const
{
	TEST_FOR_EXCEPTION(
    0!=RTOp_get_op_type_num_entries(&op_,num_values,num_indexes,num_chars)
    ,UnknownError
		,"RTOpC::get_op_type_num_entries(...): Error, "
		"RTOp_get_op_type_num_entries(...) returned != 0"
    );
}

void RTOpC::extract_op_state(
  int                             num_values
  ,primitive_value_type           value_data[]
  ,int                            num_indexes
  ,index_type                     index_data[]
  ,int                            num_chars
  ,char_type                      char_data[]
  ) const
{
  TEST_FOR_EXCEPTION(
    0!=RTOp_extract_op_state(
      &op_
      ,num_values,  value_data
      ,num_indexes, index_data
      ,num_chars,   char_data
      )
    ,UnknownError
		,"RTOpC::extract_opt_state(...): Error, "
		"RTOp_extract_opt_state(...) returned != 0"
    );
}

void RTOpC::load_op_state(
  int                           num_values
  ,const primitive_value_type   value_data[]
  ,int                          num_indexes
  ,const index_type             index_data[]
  ,int                          num_chars
  ,const char_type              char_data[]
  )
{
	TEST_FOR_EXCEPTION(
    0!=RTOp_load_op_state(
      num_values,   value_data
      ,num_indexes, index_data
      ,num_chars,   char_data
      ,&op_
      )
    ,UnknownError
		,"RTOpC::load_opt_state(...): Error, "
		"RTOp_load_opt_state(...) returned != 0"
    );
}

bool RTOpC::coord_invariant() const
{
  return false; // We have to assume this to be safe!
}

const char* RTOpC::op_name() const
{
	const char* op_name = NULL;
	TEST_FOR_EXCEPTION(
    0!=RTOp_get_op_name(&op_,&op_name)
    ,UnknownError
		,"RTOpC::get_op_name(...): Error, "
		"RTOp_op_name(...) returned != 0"
    );
	return op_name;
}

void RTOpC::apply_op(
  const int   num_vecs,       const SubVectorT<Scalar>         sub_vecs[]
  ,const int  num_targ_vecs,  const MutableSubVectorT<Scalar>  targ_sub_vecs[]
  ,ReductTarget *_reduct_obj
  ) const
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  RTOp_ReductTarget reduct_obj = RTOp_REDUCT_OBJ_NULL;
  if(_reduct_obj) reduct_obj = (*this)(*_reduct_obj);

	int k;
	Workspace<RTOp_SubVector>        c_sub_vecs(wss,num_vecs,false);
	for( k = 0; k < num_vecs; ++k ) {
		const SubVector& v = sub_vecs[k];
		RTOp_sub_vector(v.globalOffset(),v.subDim(),v.values(),v.stride(),&c_sub_vecs[k]);
	}
	Workspace<RTOp_MutableSubVector>  c_targ_sub_vecs(wss,num_targ_vecs,false);
	for( k = 0; k < num_targ_vecs; ++k ) {
		const MutableSubVector& v = targ_sub_vecs[k];
		RTOp_mutable_sub_vector(v.globalOffset(),v.subDim(),v.values(),v.stride(),&c_targ_sub_vecs[k]);
	}

	const int err = RTOp_apply_op(
		&op_
		,num_vecs,       num_vecs       ? &c_sub_vecs[0]      : (RTOp_SubVector*)NULL
		,num_targ_vecs,  num_targ_vecs  ? &c_targ_sub_vecs[0] : (RTOp_MutableSubVector*)NULL
		,reduct_obj
		);
	TEST_FOR_EXCEPTION(
		err==RTOp_ERR_INVALID_NUM_VECS, InvalidNumVecs
		,"RTOpC::apply_op(...): Error, "
		"RTOp_apply_op(...) returned RTOp_ERR_INVALID_NUM_VECS" );
	TEST_FOR_EXCEPTION(
		err==RTOp_ERR_INVALID_NUM_TARG_VECS, InvalidNumTargVecs
		,"RTOpC::apply_op(...): Error, "
		"RTOp_apply_op(...) returned RTOp_ERR_INVALID_NUM_TARG_VECS" );
	TEST_FOR_EXCEPTION(
		err!=0, UnknownError
		,"RTOpC::apply_op(...): Error, "
		"RTOp_apply_op(...) returned != 0 with unknown meaning" );

}

} // namespace RTOpPack
