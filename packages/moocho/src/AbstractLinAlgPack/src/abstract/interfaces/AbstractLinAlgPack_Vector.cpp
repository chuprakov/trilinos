// ///////////////////////////////////////////////////////////////////
// VectorWithOp.cpp

#include <assert.h>

#include <limits>
#include <ostream>

#include "AbstractLinAlgPack/include/VectorWithOp.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_dot_prod.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_get_ele.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_norms.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_num_nonzeros.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_get_sub_vector.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "RTOpPack/include/print_sub_vector.h"
#include "Range1D.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace {

// get element operator
static RTOpPack::RTOpC          get_ele_op;
static RTOpPack::ReductTarget   get_ele_targ;
// number nonzros
static RTOpPack::RTOpC          num_nonzeros_op;
static RTOpPack::ReductTarget   num_nonzeros_targ;
// Norm 1
static RTOpPack::RTOpC          norm_1_op;
static RTOpPack::ReductTarget   norm_1_targ;
// Norm 2
static RTOpPack::RTOpC          norm_2_op;
static RTOpPack::ReductTarget   norm_2_targ;
// Norm inf
static RTOpPack::RTOpC          norm_inf_op;
static RTOpPack::ReductTarget   norm_inf_targ;
// dot product operator
static RTOpPack::RTOpC          dot_op;
static RTOpPack::ReductTarget   dot_targ;
// get sub-vector operator
static RTOpPack::RTOpC          get_sub_vector_op;

// Simple class for an object that will initialize the RTOp_Server and operators.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
		// Operator and target for getting a vector element
		if(0!=RTOp_ROp_get_ele_construct( 0, &get_ele_op.op() ))
			assert(0);
		get_ele_op.reduct_obj_create(&get_ele_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_get_ele_name
			   ,&RTOp_ROp_get_ele_vtbl
			   ))
			assert(0);
		// Operator and target for norm 1
		if(0!=RTOp_ROp_num_nonzeros_construct(&num_nonzeros_op.op() ))
			assert(0);
		num_nonzeros_op.reduct_obj_create(&num_nonzeros_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_num_nonzeros_name
			   ,&RTOp_ROp_num_nonzeros_vtbl
			   ))
			assert(0);
		// Operator and target for norm 1
		if(0!=RTOp_ROp_norm_1_construct(&norm_1_op.op() ))
			assert(0);
		norm_1_op.reduct_obj_create(&norm_1_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_norm_1_name
			   ,&RTOp_ROp_norm_1_vtbl
			   ))
			assert(0);
		// Operator and target for norm 1
		if(0!=RTOp_ROp_norm_2_construct(&norm_2_op.op() ))
			assert(0);
		norm_2_op.reduct_obj_create(&norm_2_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_norm_2_name
			   ,&RTOp_ROp_norm_2_vtbl
			   ))
			assert(0);
		// Operator and target for norm 1
		if(0!=RTOp_ROp_norm_inf_construct(&norm_inf_op.op() ))
			assert(0);
		norm_inf_op.reduct_obj_create(&norm_inf_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_norm_inf_name
			   ,&RTOp_ROp_norm_inf_vtbl
			   ))
			assert(0);
		// Dot product operator and target
		if(0!=RTOp_ROp_dot_prod_construct(&dot_op.op()))
			assert(0);
		dot_op.reduct_obj_create(&dot_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_dot_prod_name
			   ,&RTOp_ROp_dot_prod_vtbl
			   ))
			assert(0);
		// Get sub-vector operator
		if(0!=RTOp_ROp_get_sub_vector_construct(1,1,&get_sub_vector_op.op()))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_get_sub_vector_name
			   ,&RTOp_ROp_get_sub_vector_vtbl
			   ))
			assert(0);
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace AbstractLinAlgPack {

VectorWithOp::VectorWithOp()
	:num_nonzeros_(std::numeric_limits<size_type>::max())   // uninitalized
	,norm_1_(-1.0), norm_2_(-1.0), norm_inf_(-1.0)
{}

size_type VectorWithOp::nz() const
{
	if( num_nonzeros_ > this->dim() ) {
		num_nonzeros_targ.reinit();
		this->apply_reduction(num_nonzeros_op,0,NULL,0,NULL,num_nonzeros_targ.obj());
		num_nonzeros_ = RTOp_ROp_num_nonzeros_val(num_nonzeros_targ.obj());
	}
	return num_nonzeros_;
}

std::ostream& VectorWithOp::output(
	std::ostream& out, bool print_dim , bool newline
	,index_type global_offset
	) const
{
	RTOp_SubVector sub_vec;
	RTOp_sub_vector_null(&sub_vec);
	this->get_sub_vector( Range1D(), SPARSE, &sub_vec );
	RTOp_SubVector  sub_vec_print = sub_vec;
	sub_vec_print.global_offset += global_offset;
	RTOpPack::output(out,sub_vec_print,print_dim,newline);
	this->free_sub_vector( &sub_vec );
	return out;
}

value_type VectorWithOp::get_ele(index_type i) const {
	assert(0==RTOp_ROp_get_ele_set_i( i, &get_ele_op.op() ));
	get_ele_targ.reinit();
	this->apply_reduction(get_ele_op,0,NULL,0,NULL,get_ele_targ.obj());
	return RTOp_ROp_get_ele_val(get_ele_targ.obj());
}

value_type VectorWithOp::norm_1() const {
	if( norm_1_ < 0.0 ) {
		norm_1_targ.reinit();
		this->apply_reduction(norm_1_op,0,NULL,0,NULL,norm_1_targ.obj());
		norm_1_ = RTOp_ROp_norm_1_val(norm_1_targ.obj());
	}
	return norm_1_;
}

value_type VectorWithOp::norm_2() const {
	if( norm_2_ < 0.0 ) {
		norm_2_targ.reinit();
		this->apply_reduction(norm_2_op,0,NULL,0,NULL,norm_2_targ.obj());
		norm_2_ = RTOp_ROp_norm_2_val(norm_2_targ.obj());
	}
	return norm_2_;
}

value_type VectorWithOp::norm_inf() const {
	if( norm_inf_ < 0.0 ) {
		norm_inf_targ.reinit();
		this->apply_reduction(norm_inf_op,0,NULL,0,NULL,norm_inf_targ.obj());
		norm_inf_ = RTOp_ROp_norm_inf_val(norm_inf_targ.obj());
	}
	return norm_inf_;
}

VectorWithOp::vec_ptr_t
VectorWithOp::create_sub_view( const Range1D& rng ) const
{
	if( rng.full_range() || ( rng.lbound() == 1 && rng.ubound() == this->dim() ) )
		return vec_ptr_t(this,false); // returned ref counted pointer does not own memory!
	return NULL;
}

void VectorWithOp::get_sub_vector(
	const Range1D& rng_in, ESparseOrDense sparse_or_dense
	, RTOp_SubVector* sub_vec
	) const
{
	const Range1D rng = rng_in.full_range() ? Range1D(1,this->dim()) : rng_in;
#ifdef _DEBUG
	THROW_EXCEPTION(
		this->dim() < rng.ubound(), std::out_of_range
		,"VectorWithOp::get_sub_vector(rng,...): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
		<<"] is not in range = [1,"<<this->dim()<<"]" );
#endif
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
	if( sub_vec->values ) {
		free( (void*)sub_vec->values  );
		free( (void*)sub_vec->indices );
	}
	RTOp_sub_vector_null( sub_vec );
	// Initialize the operator
	if(0!=RTOp_ROp_get_sub_vector_set_range( rng.lbound(), rng.ubound(), &get_sub_vector_op.op() ))
		assert(0);
	// Create the reduction object (another sub_vec)
	RTOp_ReductTarget
		reduct_obj = RTOp_REDUCT_OBJ_NULL;
	get_sub_vector_op.reduct_obj_create_raw(&reduct_obj); // This is really of type RTOp_SubVector!
	// Perform the reduction (get the sub-vector requested)
	this->apply_reduction(get_sub_vector_op,0,NULL,0,NULL,reduct_obj);
	// Set the sub-vector.  Note reduct_obj will go out of scope so the sub_vec parameter will
	// own the memory allocated within get_sub_vector_op.create_reduct_obj_raw(...).  This is okay
	//  since the client is required to call release_sub_vector(...) so release memory!
	RTOp_ROp_get_sub_vector_ESparseOrDense  _sparse_or_dense;
	if( sparse_or_dense == SPARSE )
		_sparse_or_dense = RTOP_ROP_GET_SUB_VECTOR_SPARSE;
	else if ( sparse_or_dense == DENSE )
		_sparse_or_dense = RTOP_ROP_GET_SUB_VECTOR_DENSE;
	else
		assert(0); // Error!

	*sub_vec = RTOp_ROp_get_sub_vector_val(_sparse_or_dense,reduct_obj);
	free(reduct_obj); // Now *sub_vec owns the values[] and indices[] arrays!
}

void VectorWithOp::free_sub_vector( RTOp_SubVector* sub_vec ) const
{
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
	if( sub_vec->values )
		free( (void*)sub_vec->values  );
	if( sub_vec->indices )
		free( (void*)sub_vec->indices );
	RTOp_sub_vector_null( sub_vec );
}

void VectorWithOp::has_changed() const
{
	num_nonzeros_= std::numeric_limits<size_type>::max();  // uninitalized;
	norm_1_ = norm_2_ = norm_inf_ = -1.0;
}

// Overridden from VectorBase

value_type VectorWithOp::inner_product(  const VectorBase& vec ) const
{
	using DynamicCastHelperPack::dyn_cast;

	dot_targ.reinit();
	const int num_vecs = 1;
	const VectorWithOp*
		vec_args[1] = {
			&dyn_cast<const VectorWithOp>(vec)
		};
	this->apply_reduction(dot_op,num_vecs,vec_args,0,NULL,dot_targ.obj());
	return RTOp_ROp_dot_prod_val(dot_targ.obj());
}

} // end namespace AbstractLinAlgPack
