// //////////////////////////////////////////////////////////////////////
// VectorWithOpMutable.cpp

#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "ExampleRTOpsLib/include/RTOp_TOp_assign_scalar.h"
#include "ExampleRTOpsLib/include/RTOp_TOp_assign_vectors.h"
#include "ExampleRTOpsLib/include/RTOp_TOp_axpy.h"
#include "ExampleRTOpsLib/include/RTOp_TOp_set_ele.h"
#include "RTOpPack/include/RTOpCppC.h"

namespace {

// vector scalar assignment operator
static RTOpPack::RTOpC          assign_scalar_op;
// vector assignment operator
static RTOpPack::RTOpC          assign_vec_op;
// set element operator
static RTOpPack::RTOpC          set_ele_op;
// axpy operator
static RTOpPack::RTOpC          axpy_op;

// Simple class for an object that will initialize the RTOp_Server for get_ele operator.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
		// Vector scalar assignment operator
		if(0!=RTOp_TOp_assign_scalar_construct( 0.0, &assign_scalar_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_assign_scalar_name
			   ,&RTOp_TOp_assign_scalar_vtbl
			   ))
			assert(0);
		// Vector assignment operator
		if(0!=RTOp_TOp_assign_vectors_construct( &assign_vec_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_assign_vectors_name
			   ,&RTOp_TOp_assign_vectors_vtbl
			   ))
			assert(0);
		// Set element operator
		if(0!=RTOp_TOp_set_ele_construct( 0, 0.0, &set_ele_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_set_ele_name
			   ,&RTOp_TOp_set_ele_vtbl
			   ))
			assert(0);	
		// axpy operator
		if(0!=RTOp_TOp_axpy_construct( 0.0, &axpy_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_axpy_name
			   ,&RTOp_TOp_axpy_vtbl
			   ))
			assert(0);	
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace AbstractLinAlgPack {

// VectorSpace

VectorSpaceBase::vec_ptr_t  VectorSpace::new_member() const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp_implicit_cast<vec_ptr_t::element_type>(create_member());
}

// VectorWithOpMutable

VectorWithOpMutable& VectorWithOpMutable::operator=(RTOp_value_type alpha)
{
	if(0!=RTOp_TOp_assign_scalar_set_alpha( alpha, &assign_scalar_op.op() ))
		assert(0);
	this->apply_transformation(assign_scalar_op,0,NULL,0,NULL,RTOp_REDUCT_OBJ_NULL);
	return *this;
}

VectorWithOpMutable& VectorWithOpMutable::operator=(const VectorWithOpMutable& vec)
{
	const int num_vecs = 1;
	const VectorWithOp*
		vec_args[1] = { &vec };
	this->apply_transformation(assign_vec_op,num_vecs,vec_args,0,NULL,RTOp_REDUCT_OBJ_NULL);
	return *this;
}

void VectorWithOpMutable::set_ele( RTOp_index_type i, RTOp_value_type alpha )
{
	if(0!=RTOp_TOp_set_ele_set_i_alpha( i, alpha, &set_ele_op.op() ))
		assert(0);
	this->apply_transformation(set_ele_op,0,NULL,0,NULL,RTOp_REDUCT_OBJ_NULL);
}

VectorWithOpMutable::vec_mut_ptr_t VectorWithOpMutable::clone() const
{
	vec_mut_ptr_t
		vec = this->space().create_member();
	*vec = *this;
	return vec;
}

void VectorWithOpMutable::set_sub_vector( const RTOp_SubVector& sub_vec )
{
	assert(0); //  ToDo: Implement!
}

// Overridden from VectorWithOp

VectorWithOp::vec_ptr_t
VectorWithOpMutable::create_sub_view( const LinAlgPack::Range1D& rng ) const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp_implicit_cast<vec_ptr_t::element_type>(
		const_cast<VectorWithOpMutable*>(this)->create_sub_view(rng)
		);
}

// Overridden from VectorBaseMutable

const VectorSpaceBase& VectorWithOpMutable::get_space() const
{
	return space();
}

void VectorWithOpMutable::zero()
{
	this->operator=(0.0);
}

void VectorWithOpMutable::axpy( RTOp_value_type alpha, const VectorBase& x )
{
	if( 0!=RTOp_TOp_axpy_set_alpha( alpha, &axpy_op.op() ) )
		assert(0);
	const int num_vecs = 1;
	const VectorWithOp*
		vec_args[1] = { dynamic_cast<const VectorWithOp*>(&x) };
	if( vec_args[0] == NULL )
		throw IncompatibleVectors(
			"VectorWithOp::axpy(alpha,x): Error, x is not of type VectorWithOp!" );
	this->apply_transformation(axpy_op,num_vecs,vec_args,0,NULL,RTOp_REDUCT_OBJ_NULL);
}

} // end namespace AbstractLinAlgPack
