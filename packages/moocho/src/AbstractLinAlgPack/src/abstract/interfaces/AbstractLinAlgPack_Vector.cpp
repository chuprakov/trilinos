// ///////////////////////////////////////////////////////////////////
// VectorWithOp.cpp

#include <assert.h>

#include "AbstractLinAlgPack/include/VectorWithOp.h"
#include "ExampleRTOpsLib/include/RTOp_ROp_dot_prod.h"
#include "ExampleRTOpsLib/include/RTOp_ROp_get_ele.h"
#include "RTOpPack/include/RTOpCppC.h"

namespace {

// get element operator
static RTOpPack::RTOpC          get_ele_op;
static RTOpPack::ReductTarget   get_ele_targ;
// dot product operator
static RTOpPack::RTOpC          dot_op;
static RTOpPack::ReductTarget   dot_targ;

// Simple class for an object that will initialize the RTOp_Server for get_ele operator.
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
		// Dot product operator and target
		if(0!=RTOp_ROp_dot_prod_construct(&dot_op.op()));
		dot_op.reduct_obj_create(&dot_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_dot_prod_name
			   ,&RTOp_ROp_dot_prod_vtbl
			   ))
			assert(0);
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace AbstractLinAlgPack {

RTOp_value_type VectorWithOp::get_ele(RTOp_index_type i) const {
	assert(0==RTOp_ROp_get_ele_set_i( i, &get_ele_op.op() ));
	get_ele_targ.reinit();
	this->apply_reduction(get_ele_op,0,NULL,0,NULL,get_ele_targ.obj());
	return RTOp_ROp_get_ele_val(get_ele_targ.obj());
}

void VectorWithOp::get_sub_vector( const LinAlgPack::Range1D& rng, RTOp_SubVector* sub_vec ) const
{
	assert(0); // ToDo: Implement!
}

void VectorWithOp::release_sub_vector( RTOp_SubVector* sub_vec ) const
{
	assert(0); // Todo: Implement!
}

// Overridden from VectorBase

RTOp_value_type VectorWithOp::inner_product(  const VectorBase& vec ) const
{
	dot_targ.reinit();
	const int num_vecs = 1;
	const VectorWithOp*
		vec_args[1] = { dynamic_cast<const VectorWithOp*>(&vec) };
	if(vec_args[0] == NULL)
		throw IncompatibleVectors(
			"VectorWithOp::inner_product(vec): Error, vec is not of type VectorWithOp!" );
	this->apply_reduction(dot_op,num_vecs,vec_args,0,NULL,dot_targ.obj());
	return RTOp_ROp_dot_prod_val(dot_targ.obj());
}

} // end namespace AbstractLinAlgPack
