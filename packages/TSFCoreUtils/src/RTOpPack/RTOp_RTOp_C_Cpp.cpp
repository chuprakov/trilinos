// //////////////////////////////////////////////////
// RTOp_RTOp_C_Cpp.cpp

#include "RTOp_RTOp_C_Cpp.h"
#include "RTOpCpp.hpp"

extern "C" {

//
// Static Implementation functions
//

// Functions for obj_data_vtbl

static int get_op_type_num_entries(
	const struct RTOp_obj_type_vtbl_t* vtbl
	,const void* op_data
	,int* num_values
	,int* num_indexes
	,int* num_chars
	)
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->get_op_type_num_entries(
			num_values, num_indexes, num_chars );
	}
	catch(...) {
		err = -1;
	}
	return err;
}

static int extract_op_state(
	const struct RTOp_obj_type_vtbl_t* vtbl
	,const void *       dummy
	,void *             op_data
	,int                num_values
	,RTOp_value_type    value_data[]
	,int                num_indexes
	,RTOp_index_type    index_data[]
	,int                num_chars
	,RTOp_char_type     char_data[]
	)
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->extract_op_state(
			num_values,   value_data
			,num_indexes, index_data
			,num_chars,   char_data );
	}
	catch(...) {
		err = -1;
	}
	return err;
}

static int load_op_state(
	const struct RTOp_obj_type_vtbl_t* vtbl
	,const void *            dummy
	,int                     num_values
	,const RTOp_value_type   value_data[]
	,int                     num_indexes
	,const RTOp_index_type   index_data[]
	,int                     num_chars
	,const RTOp_char_type    char_data[]
	,void **                 op_data
	)
{
	int err = 0;
	try {
		if(*op_data == NULL)
			vtbl->obj_create(NULL,NULL,op_data);
		reinterpret_cast<RTOpPack::RTOp*>(op_data)->load_op_state(
			num_values,   value_data
			,num_indexes, index_data
			,num_chars,   char_data );
	}
	catch(...) {
		err = -1;
	}
	return err;
	
}

// Functions for reduct_vtbl

static int get_reduct_type_num_entries(
	const struct RTOp_obj_type_vtbl_t* vtbl
	,const void* op_data
	,int* num_values
	,int* num_indexes
	,int* num_chars
	)
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->get_reduct_type_num_entries(
			num_values, num_indexes, num_chars );
	}
	catch(...) {
		err = -1;
	}
	return err;
}

static int reduct_obj_create(
	const struct RTOp_obj_type_vtbl_t* vtbl, const void* op_data, void** reduct_obj )
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->reduct_obj_create_raw(reduct_obj);
	}
	catch(...) {
		err = -1;
	}
	return err;
}

static int reduct_obj_reinit(
	const struct RTOp_obj_type_vtbl_t* vtbl, const void* op_data, void* reduct_obj )
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->reduct_obj_reinit(reduct_obj);
	}
	catch(...) {
		err = -1;
	}
	return err;
}

static int reduct_obj_free(
	const struct RTOp_obj_type_vtbl_t* vtbl, const void* op_data, void** reduct_obj )
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->reduct_obj_free(reduct_obj);
	}
	catch(...) {
		err = -1;
	}
	return err;
}

static int extract_reduct_obj_state(
	const struct RTOp_obj_type_vtbl_t* vtbl
	,const void *       op_data
	,void *             reduct_obj
	,int                num_values
	,RTOp_value_type    value_data[]
	,int                num_indexes
	,RTOp_index_type    index_data[]
	,int                num_chars
	,RTOp_char_type     char_data[]
	)
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->extract_reduct_obj_state(
			const_cast<void*>(reduct_obj)
			,num_values,   value_data
			,num_indexes, index_data
			,num_chars,   char_data );
	}
	catch(...) {
		err = -1;
	}
	return err;
}

static int load_reduct_obj_state(
	const struct RTOp_obj_type_vtbl_t* vtbl
	,const void *            op_data
	,int                     num_values
	,const RTOp_value_type   value_data[]
	,int                     num_indexes
	,const RTOp_index_type   index_data[]
	,int                     num_chars
	,const RTOp_char_type    char_data[]
	,void **                 reduct_obj
	)
{
	int err = 0;
	try {
		assert(*reduct_obj);
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->load_reduct_obj_state(
			num_values,   value_data
			,num_indexes, index_data
			,num_chars,   char_data
			,reinterpret_cast<RTOp_ReductTarget>(*reduct_obj) );
	}
	catch(...) {
		err = -1;
	}
	return err;
}

} // end extern "C"

static const struct RTOp_obj_type_vtbl_t reduct_obj_vtbl =
{
	get_reduct_type_num_entries
	,reduct_obj_create
	,reduct_obj_reinit
	,reduct_obj_free
	,extract_reduct_obj_state
	,load_reduct_obj_state
};

//  other static functions

extern "C" {

static int apply_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* op_data
	, const int num_vecs, const struct RTOp_SubVector sub_vecs[]
	, const int num_targ_vecs, const struct RTOp_MutableSubVector targ_sub_vecs[]
	, RTOp_ReductTarget reduct_obj )
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	int err = 0;
	try {
		int k;
		wsp::Workspace<RTOpPack::SubVector>         cpp_sub_vecs(wss,num_vecs);
		for( k = 0; k < num_vecs; ++k ) {
			const RTOp_SubVector& v = sub_vecs[k];
			cpp_sub_vecs[k].initialize(v.global_offset,v.sub_dim,v.values,v.values_stride);
		}
		wsp::Workspace<RTOpPack::MutableSubVector>  cpp_targ_sub_vecs(wss,num_targ_vecs);
		for( k = 0; k < num_targ_vecs; ++k ) {
			const RTOp_MutableSubVector& v = targ_sub_vecs[k];
			cpp_targ_sub_vecs[k].initialize(v.global_offset,v.sub_dim,v.values,v.values_stride);
		}
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->apply_op(
			num_vecs,        num_vecs      ? &cpp_sub_vecs[0]      : (RTOpPack::SubVector*)NULL
			,num_targ_vecs,  num_targ_vecs ? &cpp_targ_sub_vecs[0] : (RTOpPack::MutableSubVector*)NULL
			,reduct_obj	);
	}
	catch ( const RTOpPack::InvalidNumVecs& ) {
		err = RTOp_ERR_INVALID_NUM_VECS;
	}
	catch( const RTOpPack::InvalidNumTargVecs& ) {
		err = RTOp_ERR_INVALID_NUM_TARG_VECS;
	}
	catch( const RTOpPack::IncompatibleVecs& ) {
		err = RTOp_ERR_INCOMPATIBLE_VECS;
	}
	catch(...) {
		err = RTOp_ERR_INVALID_USAGE;
	}
	return err;
}

static int reduce_reduct_objs(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* op_data
	, RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj )
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->reduce_reduct_objs(
			in_reduct_obj, inout_reduct_obj );
	}
	catch(...) {

		err = -1;
	}
	return err;
}

static int get_reduct_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* op_data
	, RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr )
{
	int err = 0;
	try {
		reinterpret_cast<const RTOpPack::RTOp*>(op_data)->get_reduct_op(reduct_op_func_ptr);
	}
	catch(...) {

		err = -1;
	}
	return err;
}

} // end extern "C"

// Public functions

int RTOp_create_C_Cpp_vtbl( 
	RTOp_op_create_func_ptr_t op_create, RTOp_op_free_func_ptr_t op_free
	,struct RTOp_RTOp_vtbl_t* vtbl )
{
#ifdef RTOp_DEBUG
	assert(op_create);
	assert(op_free);
#endif

	assert( vtbl );
	RTOp_obj_type_vtbl_t
		*obj_data_vtbl = new RTOp_obj_type_vtbl_t;
	assert(obj_data_vtbl);
	obj_data_vtbl->get_obj_type_num_entries    = get_op_type_num_entries;
	obj_data_vtbl->obj_create                  = op_create;
	obj_data_vtbl->obj_reinit                  = NULL;
	obj_data_vtbl->obj_free                    = op_free;   
	obj_data_vtbl->extract_state               = extract_op_state;
	obj_data_vtbl->load_state                  = load_op_state;   
	vtbl->obj_data_vtbl      = obj_data_vtbl;
	vtbl->reduct_vtbl        = &reduct_obj_vtbl;
	vtbl->op_name            = NULL;
	vtbl->reduct_obj_reinit  = NULL;
	vtbl->apply_op           = apply_op;
	vtbl->reduce_reduct_objs = reduce_reduct_objs;
	vtbl->get_reduct_op      = get_reduct_op;
	return 0;
}

int RTOp_free_C_Cpp_vtbl( struct RTOp_RTOp_vtbl_t* vtbl )
{
	assert(vtbl->obj_data_vtbl);
	delete vtbl->obj_data_vtbl;
	vtbl->obj_data_vtbl      = NULL;
	vtbl->reduct_vtbl        = NULL;
	vtbl->reduct_obj_reinit  = NULL;
	vtbl->apply_op           = NULL;
	vtbl->reduce_reduct_objs = NULL;
	vtbl->get_reduct_op      = NULL;
	return 0;
}
