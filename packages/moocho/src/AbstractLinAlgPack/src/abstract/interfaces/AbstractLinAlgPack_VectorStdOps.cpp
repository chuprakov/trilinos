// ////////////////////////////////////////////////////////////////////
// VectorStdOps.cpp

#include <assert.h>

#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_max_near_feas_step.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_num_bounded.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_sum.h"
#include "RTOpStdOpsLib/include/RTOp_TOp_ele_wise_prod.h"
#include "RTOpStdOpsLib/include/RTOp_TOp_force_in_bounds.h"
#include "RTOpStdOpsLib/include/RTOp_TOp_random_vector.h"
#include "RTOpStdOpsLib/include/RTOp_TOp_scale_vector.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "ThrowException.h"

namespace {

// sum
static RTOpPack::RTOpC          sum_op;
static RTOpPack::ReductTarget   sum_targ;
// maximum near feasible step
static RTOpPack::RTOpC          max_near_feas_step_op;
static RTOpPack::ReductTarget   max_near_feas_step_targ;
// number of bounded elements
static RTOpPack::RTOpC          num_bounded_op;
static RTOpPack::ReductTarget   num_bounded_targ;
// scale vector
static RTOpPack::RTOpC          scale_vector_op;
// random vector
static RTOpPack::RTOpC          random_vector_op;
// element-wise product
static RTOpPack::RTOpC          ele_wise_prod_op;
// force in bounds
static RTOpPack::RTOpC          force_in_bounds_op;

// Simple class for an object that will initialize the RTOp_Server.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
		// Operator and target obj for sum
		if(0!=RTOp_ROp_sum_construct(&sum_op.op() ))
			assert(0);
		sum_op.reduct_obj_create(&sum_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_sum_name
			   ,&RTOp_ROp_sum_vtbl
			   ))
			assert(0);
		// Operator and target obj for max_near_feas_step
		if(0!=RTOp_ROp_max_near_feas_step_construct(0.0,&max_near_feas_step_op.op() ))
			assert(0);
		max_near_feas_step_op.reduct_obj_create(&max_near_feas_step_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_max_near_feas_step_name
			   ,&RTOp_ROp_max_near_feas_step_vtbl
			   ))
			assert(0);
		// Operator and target obj for num_bounded
		if(0!=RTOp_ROp_num_bounded_construct(0.0,&num_bounded_op.op() ))
			assert(0);
		num_bounded_op.reduct_obj_create(&num_bounded_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_num_bounded_name
			   ,&RTOp_ROp_num_bounded_vtbl
			   ))
			assert(0);
		// Operator scale_vector
		if(0!=RTOp_TOp_scale_vector_construct( 0.0, &scale_vector_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_scale_vector_name
			   ,&RTOp_TOp_scale_vector_vtbl
			   ))
			assert(0);
		// Operator random_vector
		if(0!=RTOp_TOp_random_vector_construct( 0.0, 0.0, &random_vector_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_random_vector_name
			   ,&RTOp_TOp_random_vector_vtbl
			   ))
			assert(0);
		// Operator ele_wise_prod
		if(0!=RTOp_TOp_ele_wise_prod_construct( 0.0, &ele_wise_prod_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_ele_wise_prod_name
			   ,&RTOp_TOp_ele_wise_prod_vtbl
			   ))
			assert(0);
		// Operator force_in_bounds
		if(0!=RTOp_TOp_force_in_bounds_construct( &force_in_bounds_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_force_in_bounds_name
			   ,&RTOp_TOp_force_in_bounds_vtbl
			   ))
			assert(0);
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

AbstractLinAlgPack::value_type
AbstractLinAlgPack::sum( const VectorWithOp& v_rhs )
{
	sum_targ.reinit();
	v_rhs.apply_reduction(sum_op,0,NULL,0,NULL,sum_targ.obj() );
	return RTOp_ROp_sum_val(sum_targ.obj());
}

std::pair<AbstractLinAlgPack::value_type,AbstractLinAlgPack::value_type>
AbstractLinAlgPack::max_near_feas_step(
	const VectorWithOp& x, const VectorWithOp& d
	,const VectorWithOp& xl, const VectorWithOp& xu
	,value_type max_bnd_viol
	)
{
	const int num_vecs = 3;
	const VectorWithOp*
		vecs[num_vecs] = { &x, &d, &xu };
	assert(0==RTOp_ROp_max_near_feas_step_set_beta( max_bnd_viol, &max_near_feas_step_op.op() ));
	max_near_feas_step_targ.reinit();
	xl.apply_reduction(
		max_near_feas_step_op, num_vecs, vecs, 0, NULL
		,max_near_feas_step_targ.obj() );
	RTOp_ROp_max_near_feas_step_reduct_obj_t
		u = RTOp_ROp_max_near_feas_step_val(max_near_feas_step_targ.obj());
	return std::pair<value_type,value_type>(u.alpha_pos,u.alpha_neg);
}

AbstractLinAlgPack::size_type
AbstractLinAlgPack:: num_bounded(
	const VectorWithOp& xl, const VectorWithOp& xu
	,value_type inf_bound
	)
{
	const int num_vecs = 1;
	const VectorWithOp*
		vecs[num_vecs] = { &xu };
	assert(0==RTOp_ROp_num_bounded_set_inf_bnd( inf_bound, &num_bounded_op.op() ));
	num_bounded_targ.reinit();
	xl.apply_reduction(
		num_bounded_op, num_vecs, vecs, 0, NULL
		,num_bounded_targ.obj() );
	return RTOp_ROp_num_bounded_val(num_bounded_targ.obj());
}

void AbstractLinAlgPack::Vt_S(
	VectorWithOpMutable* v_lhs, const value_type& alpha )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"Vt_S(...), Error");
#endif
	assert(0==RTOp_TOp_scale_vector_set_alpha( alpha, &scale_vector_op.op() ));
	v_lhs->apply_transformation(scale_vector_op,0,NULL,0,NULL,RTOp_REDUCT_OBJ_NULL);
}

void AbstractLinAlgPack::ele_wise_prod(
	const value_type& alpha, const VectorWithOp& v_rhs1, const VectorWithOp& v_rhs2
	, VectorWithOpMutable* v_lhs )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v_lhs==NULL,std::logic_error,"force_in_bounds(...), Error");
#endif
	assert(0==RTOp_TOp_ele_wise_prod_set_alpha(alpha,&ele_wise_prod_op.op()));
	const int num_vecs = 2;
	const VectorWithOp*
		vecs[num_vecs] = { &v_rhs1, &v_rhs2 };
	v_lhs->apply_transformation(
		ele_wise_prod_op, num_vecs, vecs, 0, NULL, RTOp_REDUCT_OBJ_NULL );
}

void AbstractLinAlgPack::seed_random_vector_generator( unsigned int s )
{
	srand(s);
}

void AbstractLinAlgPack::random_vector(
	value_type l, value_type u, VectorWithOpMutable* v )
{
#ifdef _DEBUG
	THROW_EXCEPTION(v==NULL,std::logic_error,"Vt_S(...), Error");
#endif
	assert(0==RTOp_TOp_random_vector_set_bounds( l, u, &random_vector_op.op() ));
	v->apply_transformation(random_vector_op,0,NULL,0,NULL,RTOp_REDUCT_OBJ_NULL);
}

void AbstractLinAlgPack::force_in_bounds(
	const VectorWithOp& xl, const VectorWithOp& xu
	, VectorWithOpMutable* x )
{
#ifdef _DEBUG
	THROW_EXCEPTION(x==NULL,std::logic_error,"force_in_bounds(...), Error");
#endif
	const int num_vecs = 2;
	const VectorWithOp*
		vecs[num_vecs] = { &xl, &xu };
	x->apply_transformation(
		force_in_bounds_op, num_vecs, vecs, 0, NULL, RTOp_REDUCT_OBJ_NULL );
}
