// ////////////////////////////////////////////////////////////////////////
// RTOpCppC.cpp

#include "RTOpCppC.hpp"
#include "Teuchos_TestForException.hpp"

namespace RTOpPack {

// Members for RTOpC

RTOpC::RTOpC( const RTOp_RTOp_vtbl_t* vtbl )
	: vtbl_(vtbl)
{
	op_.vtbl     = vtbl;
	op_.obj_data = NULL;
}

RTOpC::~RTOpC()
{
	if(op_.obj_data)
		RTOp_free_op( &op_ );
}

// Overridden from RTOp

const char* RTOpC::op_name() const
{
	const char* op_name = NULL;
	const int err = RTOp_get_op_name(&op_,&op_name);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::get_op_name(...): Error, "
		"RTOp_op_name(...) returned != 0" );
	return op_name;
}

op_create_func_t RTOpC::get_op_create_func() const
{
#ifdef RTOp_DEBUG
	assert(vtbl_->obj_data_vtbl->obj_create);
#endif
	return vtbl_->obj_data_vtbl->obj_create;
}

op_free_func_t RTOpC::get_op_free_func() const
{
#ifdef RTOp_DEBUG
	assert(vtbl_->obj_data_vtbl->obj_free);
#endif
	return vtbl_->obj_data_vtbl->obj_free;
}

void RTOpC::get_op_type_num_entries(
	int*  num_values
	,int* num_indexes
	,int* num_chars
	) const
{
	const int err = RTOp_get_op_type_num_entries(&op_,num_values,num_indexes,num_chars);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::get_op_type_num_entries(...): Error, "
		"RTOp_get_op_type_num_entries(...) returned != 0" );
}

void RTOpC::extract_op_state(
	int               num_values
	,RTOp_value_type  value_data[]
	,int              num_indexes
	,RTOp_index_type  index_data[]
	,int              num_chars
	,RTOp_char_type   char_data[]
	) const
{
	const int err = RTOp_extract_op_state(
		&op_
		,num_values,  value_data
		,num_indexes, index_data
		,num_chars,   char_data
		);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::extract_opt_state(...): Error, "
		"RTOp_extract_opt_state(...) returned != 0" );
}

void RTOpC::load_op_state(
	int                       num_values
	,const RTOp_value_type    value_data[]
	,int                      num_indexes
	,const RTOp_index_type    index_data[]
	,int                      num_chars
	,const RTOp_char_type     char_data[]
	)
{
	const int err = RTOp_load_op_state(
		num_values,   value_data
		,num_indexes, index_data
		,num_chars,   char_data
		,&op_
		);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::load_opt_state(...): Error, "
		"RTOp_load_opt_state(...) returned != 0" );
}

void RTOpC::get_reduct_type_num_entries(
	int*   num_values
	,int*  num_indexes
	,int*  num_chars
	) const
{
	const int err = RTOp_get_reduct_type_num_entries(&op_,num_values,num_indexes,num_chars);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::get_reduct_type_num_entries(...): Error, "
		"RTOp_get_reduct_type_num_entries(...) returned != 0" );
}

void RTOpC::reduct_obj_create_raw( RTOp_ReductTarget* reduct_obj ) const
{
	const int err = RTOp_reduct_obj_create(&op_,reduct_obj);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::reduct_obj_create(...): Error, "
		"RTOp_reduct_obj_create(...) returned != 0" );
}

void RTOpC::reduct_obj_reinit( RTOp_ReductTarget reduct_obj ) const
{
	const int err = RTOp_reduct_obj_reinit(&op_,reduct_obj);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::reduct_obj_reinit(...): Error, "
		"RTOp_reduct_obj_reinit(...) returned != 0" );
}

void RTOpC::reduct_obj_free( RTOp_ReductTarget* reduct_obj ) const
{
	const int err = RTOp_reduct_obj_free(&op_,reduct_obj);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::reduct_obj_free(...): Error, "
		"RTOp_reduct_obj_free(...) returned != 0" );
}

void RTOpC::extract_reduct_obj_state(
	const RTOp_ReductTarget   reduct_obj
	,int                      num_values
	,RTOp_value_type          value_data[]
	,int                      num_indexes
	,RTOp_index_type          index_data[]
	,int                      num_chars
	,RTOp_char_type           char_data[]
	) const
{
	const int err = RTOp_extract_reduct_obj_state(
		&op_, reduct_obj
		,num_values,  value_data
		,num_indexes, index_data
		,num_chars,   char_data
		);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::extract_reduct_obj_state(...): Error, "
		"RTOp_extract_reduct_obj_state(...) returned != 0" );
}

void RTOpC::load_reduct_obj_state(
	int                       num_values
	,const RTOp_value_type    value_data[]
	,int                      num_indexes
	,const RTOp_index_type    index_data[]
	,int                      num_chars
	,const RTOp_char_type     char_data[]
	,RTOp_ReductTarget        reduct_obj
	) const
{
	const int err = RTOp_load_reduct_obj_state(
		&op_
		,num_values,  value_data
		,num_indexes, index_data
		,num_chars,   char_data
		,reduct_obj
		);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::load_reduct_obj_state(...): Error, "
		"RTOp_load_reduct_obj_state(...) returned != 0" );
}

void RTOpC::apply_op(
	const int   num_vecs,       const SubVector         sub_vecs[]
	,const int  num_targ_vecs,  const MutableSubVector  targ_sub_vecs[]
	,RTOp_ReductTarget reduct_obj
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	int k;
	wsp::Workspace<RTOp_SubVector>        c_sub_vecs(wss,num_vecs);
	for( k = 0; k < num_vecs; ++k ) {
		const SubVector& v = sub_vecs[k];
		RTOp_sub_vector(v.globalOffset(),v.subDim(),v.values(),v.stride(),&c_sub_vecs[k]);
	}
	wsp::Workspace<RTOp_MutableSubVector>  c_targ_sub_vecs(wss,num_targ_vecs);
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
		err == RTOp_ERR_INVALID_NUM_VECS, InvalidNumVecs
		,"RTOpC::apply_op(...): Error, "
		"RTOp_apply_op(...) returned RTOp_ERR_INVALID_NUM_VECS" );
	TEST_FOR_EXCEPTION(
		err == RTOp_ERR_INVALID_NUM_TARG_VECS, InvalidNumTargVecs
		,"RTOpC::apply_op(...): Error, "
		"RTOp_apply_op(...) returned RTOp_ERR_INVALID_NUM_TARG_VECS" );
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::apply_op(...): Error, "
		"RTOp_apply_op(...) returned != 0 with unknown meaning" );
}

void RTOpC::reduce_reduct_objs(
	RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj
	) const
{
	const int err = RTOp_reduce_reduct_objs(&op_,in_reduct_obj,inout_reduct_obj);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::reduce_reduct_objs(...): Error, "
		"RTOp_reduce_reduct_objs(...) returned != 0" );
}

void RTOpC::get_reduct_op( RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr ) const
{
	const int err = RTOp_get_reduct_op(&op_,reduct_op_func_ptr);
	TEST_FOR_EXCEPTION(
		err != 0, UnknownError
		,"RTOpC::get_reduct_op(...): Error, "
		"RTOp_get_reduct_op(...) returned != 0" );
}

} // end namespace RTOpPack
