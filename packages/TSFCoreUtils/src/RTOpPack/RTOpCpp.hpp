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

// //////////////////////////////////////////////////////////////////////////
// RTOpCpp.hpp

#ifndef RTOPPACK_RTOP_CPP_HPP
#define RTOPPACK_RTOP_CPP_HPP

#include "RTOpCppDecl.hpp"
#include "RTOpCppC.hpp"
#include "RTOp_RTOp_C_Cpp.h"
#include "WorkspacePack.hpp"
#include "Teuchos_TestForException.hpp"

namespace RTOpPack {

// RTOpT

template<class Scalar>
op_create_func_t RTOpT<Scalar>::get_op_create_func() const
{
	return NULL;
}

template<class Scalar>
op_free_func_t RTOpT<Scalar>::get_op_free_func() const
{
	return NULL;
}

template<class Scalar>
void RTOpT<Scalar>::get_op_type_num_entries(
	int*  num_values
	,int* num_indexes
	,int* num_chars
	) const
{
#ifdef RTOp_DEBUG
	assert( num_values );
	assert( num_indexes );
	assert( num_chars );
#endif
	*num_values  = 0;
	*num_indexes = 0;
	*num_chars   = 0;
}

template<class Scalar>
void RTOpT<Scalar>::extract_op_state(
	int               num_values
	,Scalar           value_data[]
	,int              num_indexes
	,RTOp_index_type  index_data[]
	,int              num_chars
	,RTOp_char_type   char_data[]
	) const
{
#ifdef RTOp_DEBUG
	assert( num_values  == 0 );
	assert( num_indexes == 0 );
	assert( num_chars   == 0 );
#endif
}

template<class Scalar>
void RTOpT<Scalar>::load_op_state(
	int                       num_values
	,const Scalar             value_data[]
	,int                      num_indexes
	,const RTOp_index_type    index_data[]
	,int                      num_chars
	,const RTOp_char_type     char_data[]
	)
{
#ifdef RTOp_DEBUG
	assert( num_values  == 0 );
	assert( num_indexes == 0 );
	assert( num_chars   == 0 );
#endif
}

template<class Scalar>
void RTOpT<Scalar>::get_reduct_type_num_entries(
	int*   num_values
	,int*  num_indexes
	,int*  num_chars
	) const
{
#ifdef RTOp_DEBUG
	assert( num_values );
	assert( num_indexes );
	assert( num_chars );
#endif
	*num_values  = 0;
	*num_indexes = 0;
	*num_chars   = 0;
}

template<class Scalar>
void RTOpT<Scalar>::reduct_obj_create( ReductTargetT<Scalar>* reduct_obj ) const
{
	RTOp_ReductTarget reduct_obj_raw;
	this->reduct_obj_create_raw( &reduct_obj_raw );
	reduct_obj->initialize(this,reduct_obj_raw,true); // clean up memory!
}

template<class Scalar>
bool RTOpT<Scalar>::coord_invariant() const
{
	return true;
}

template<class Scalar>
void RTOpT<Scalar>::reduce_reduct_objs(
	RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj
	) const
{
#ifdef RTOp_DEBUG
	assert(in_reduct_obj == RTOp_REDUCT_OBJ_NULL);
	assert(inout_reduct_obj == RTOp_REDUCT_OBJ_NULL);
#endif
}

template<class Scalar>
void RTOpT<Scalar>::get_reduct_op( RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr ) const
{
#ifdef RTOp_DEBUG
	assert( reduct_op_func_ptr );
#endif
	*reduct_op_func_ptr = NULL;
}

template<class Scalar>
void RTOpT<Scalar>::reduct_obj_create_raw( RTOp_ReductTarget* reduct_obj ) const
{
#ifdef RTOp_DEBUG
	assert(reduct_obj);
#endif
	*reduct_obj = RTOp_REDUCT_OBJ_NULL;
}

template<class Scalar>
void RTOpT<Scalar>::reduct_obj_reinit( RTOp_ReductTarget reduct_obj ) const
{
#ifdef RTOp_DEBUG
	assert(reduct_obj==RTOp_REDUCT_OBJ_NULL);
#endif
}

template<class Scalar>
void RTOpT<Scalar>::reduct_obj_free( RTOp_ReductTarget* reduct_obj ) const
{
#ifdef RTOp_DEBUG
	assert(*reduct_obj==RTOp_REDUCT_OBJ_NULL);
#endif
}

template<class Scalar>
void RTOpT<Scalar>::extract_reduct_obj_state(
	const RTOp_ReductTarget   reduct_obj
	,int                      num_values
	,Scalar          value_data[]
	,int                      num_indexes
	,RTOp_index_type          index_data[]
	,int                      num_chars
	,RTOp_char_type           char_data[]
	) const

{
#ifdef RTOp_DEBUG
	assert(reduct_obj==RTOp_REDUCT_OBJ_NULL);
	assert( num_values  == 0 );
	assert( num_indexes == 0 );
	assert( num_chars   == 0 );
#endif
}

template<class Scalar>
void RTOpT<Scalar>::load_reduct_obj_state(
	int                      num_values
	,const Scalar            value_data[]
	,int                     num_indexes
	,const RTOp_index_type   index_data[]
	,int                     num_chars
	,const RTOp_char_type    char_data[]
	,RTOp_ReductTarget       reduct_obj
	) const
{
#ifdef RTOp_DEBUG
	assert( num_values  == 0 );
	assert( num_indexes == 0 );
	assert( num_chars   == 0 );
	assert(reduct_obj==RTOp_REDUCT_OBJ_NULL);
#endif
}

template<class Scalar>
RTOpT<Scalar>& RTOpT<Scalar>::operator=(const RTOpT<Scalar>& op)
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	int num_values = 0, num_indexes = 0, num_chars = 0;
	op.get_op_type_num_entries( &num_values, &num_indexes, &num_chars );
	wsp::Workspace<Scalar>           value_data(wss,num_values);
	wsp::Workspace<RTOp_index_type>  index_data(wss,num_indexes);
	wsp::Workspace<RTOp_char_type>   char_data(wss,num_chars);
	op.extract_op_state(
		num_values,   num_values  ? &value_data[0] : NULL
		,num_indexes, num_indexes ? &index_data[0] : NULL
		,num_chars,   num_chars   ? &char_data[0]  : NULL
		);
	this->load_op_state(
		num_values,   num_values  ? &value_data[0] : NULL
		,num_indexes, num_indexes ? &index_data[0] : NULL
		,num_chars,   num_chars   ? &char_data[0]  : NULL
		);
	return *this;
}

} // end namespace RTOpPack

#endif // RTOPPACK_RTOP_CPP_HPP
