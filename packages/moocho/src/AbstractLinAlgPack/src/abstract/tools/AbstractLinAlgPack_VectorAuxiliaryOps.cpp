// ////////////////////////////////////////////////////////////////////
// VectorAuxiliaryOps.cpp
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

#include <assert.h>

#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_max.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_max_near_feas_step.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_max_rel_step.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_num_bounded.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_log_bound_barrier.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_max_inequ_viol.h"
#include "RTOpStdOpsLib/include/RTOp_TOp_force_in_bounds.h"
#include "RTOpStdOpsLib/include/RTOp_TOp_max_vec_scalar.h"
#include "RTOpStdOpsLib/include/RTOp_TOp_max_abs_vec_scalar.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "ThrowException.h"

namespace {

// maximum near feasible step
static RTOpPack::RTOpC          max_near_feas_step_op;
static RTOpPack::ReductTarget   max_near_feas_step_targ;
// maximum relative step
static RTOpPack::RTOpC          max_rel_step_op;
static RTOpPack::ReductTarget   max_rel_step_targ;
// number of bounded elements
static RTOpPack::RTOpC          num_bounded_op;
static RTOpPack::ReductTarget   num_bounded_targ;
// force in bounds
static RTOpPack::RTOpC          force_in_bounds_op;

// Simple class for an object that will initialize the RTOp_Server.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
		// Operator and target obj for max_near_feas_step
		if(0>RTOp_ROp_max_near_feas_step_construct(0.0,&max_near_feas_step_op.op() ))
			assert(0);
		max_near_feas_step_op.reduct_obj_create(&max_near_feas_step_targ);
		if(0>RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_max_near_feas_step_name
			   ,&RTOp_ROp_max_near_feas_step_vtbl
			   ))
			assert(0);
		// Operator and target obj for max_rel_step
		if(0>RTOp_ROp_max_rel_step_construct(&max_rel_step_op.op() ))
			assert(0);
		max_rel_step_op.reduct_obj_create(&max_rel_step_targ);
		if(0>RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_max_rel_step_name
			   ,&RTOp_ROp_max_rel_step_vtbl
			   ))
			assert(0);
		// Operator and target obj for num_bounded
		if(0>RTOp_ROp_num_bounded_construct(0.0,&num_bounded_op.op() ))
			assert(0);
		num_bounded_op.reduct_obj_create(&num_bounded_targ);
		if(0>RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_num_bounded_name
			   ,&RTOp_ROp_num_bounded_vtbl
			   ))
			assert(0);
		// Operator force_in_bounds
		if(0>RTOp_TOp_force_in_bounds_construct( &force_in_bounds_op.op() ))
			assert(0);
		if(0>RTOp_Server_add_op_name_vtbl(
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
AbstractLinAlgPack::max( const VectorWithOp& v )
{
	RTOpPack::RTOpC          op;
	RTOpPack::ReductTarget   reduct_obj;
	assert(0==RTOp_ROp_max_construct(&op.op()));
	op.reduct_obj_create(&reduct_obj);
	v.apply_reduction( op, 0, NULL, 0, NULL,reduct_obj.obj() );
	return RTOp_ROp_max_val(reduct_obj.obj());
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

AbstractLinAlgPack::value_type
AbstractLinAlgPack::max_rel_step(
	const VectorWithOp& x, const VectorWithOp& d
	)
{
	const int num_vecs = 1;
	const VectorWithOp*
		vecs[num_vecs] = { &d };
	max_rel_step_targ.reinit();
	x.apply_reduction(
		max_rel_step_op, num_vecs, vecs, 0, NULL
		,max_rel_step_targ.obj() );
	return RTOp_ROp_max_rel_step_val(max_rel_step_targ.obj());
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

AbstractLinAlgPack::value_type
AbstractLinAlgPack::log_bound_barrier(
	const VectorWithOp    &x
	,const VectorWithOp   &xl
	,const VectorWithOp   &xu
	)
{
	RTOpPack::RTOpC          op;
	RTOpPack::ReductTarget   reduct_obj;
	assert(0==RTOp_ROp_log_bound_barrier_construct(&op.op()));
	op.reduct_obj_create(&reduct_obj);
	const int num_vecs = 2;
	const VectorWithOp*
		vecs[num_vecs] = { &xl, &xu };
	x.apply_reduction(
		op, num_vecs, vecs, 0, NULL
		,reduct_obj.obj()
		);
	return RTOp_ROp_log_bound_barrier_val(reduct_obj.obj());
}

bool AbstractLinAlgPack::max_inequ_viol(
	const AbstractLinAlgPack::VectorWithOp   &v
	,const AbstractLinAlgPack::VectorWithOp  &vL
	,const AbstractLinAlgPack::VectorWithOp  &vU
	,AbstractLinAlgPack::size_type           *max_viol_i
	,AbstractLinAlgPack::value_type          *max_viol
	,AbstractLinAlgPack::value_type          *v_i
	,int                                     *bnd_type
	,AbstractLinAlgPack::value_type          *vLU_i
	)
{
	RTOpPack::RTOpC          op;
	RTOpPack::ReductTarget   reduct_obj;
	RTOp_ROp_max_inequ_viol_construct(&op.op());
	op.reduct_obj_create(&reduct_obj);
	const int num_vecs = 2;
	const VectorWithOp*
		vecs[num_vecs] = { &vL, &vU };
	v.apply_reduction(
		op, num_vecs, vecs, 0, NULL
		,reduct_obj.obj()
		);
	const RTOp_ROp_max_inequ_viol_reduct_obj_t
		ro = RTOp_ROp_max_inequ_viol_val(reduct_obj.obj());
	*max_viol_i = ro.max_viol_i;
	*max_viol   = ro.max_viol;
	*v_i        = ro.v_i;
	*bnd_type   = ro.bnd_type;
	*vLU_i      = ro.vLU_i;
	return *max_viol_i > 0.0;
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

void AbstractLinAlgPack::max_vec_scalar(
	value_type              min_ele
	,VectorWithOpMutable    *y
	)
{
	RTOpPack::RTOpC op;
	assert(0==RTOp_TOp_max_vec_scalar_construct(min_ele,&op.op()));
	y->apply_transformation( op, 0, NULL, 0, NULL, RTOp_REDUCT_OBJ_NULL );
}

void AbstractLinAlgPack::max_abs_vec_scalar(
	value_type              min_ele
	,VectorWithOpMutable    *y
	)
{
	RTOpPack::RTOpC op;
	assert(0==RTOp_TOp_max_abs_vec_scalar_construct(min_ele,&op.op()));
	y->apply_transformation( op, 0, NULL, 0, NULL, RTOp_REDUCT_OBJ_NULL );
}
