// //////////////////////////////////////////
// ExampleNLPObjGradient.cpp
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

#include <stdexcept>

#include "ExampleNLPObjGradient.h"
#include "ExampleNLPFirstOrderDirectRTOps.h"
#include "AbstractLinAlgPack/include/BasisSystemCompositeStd.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "Range1D.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"
#include "AbstractFactoryStd.h"

namespace {

// Evaluate c(x)
static RTOpPack::RTOpC          explnlp2_c_eval_op;

// Simple class for an object that will initialize the RTOp_Server.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
		// Operator evaluate c(x)
		if(0!=RTOp_TOp_explnlp2_c_eval_construct( &explnlp2_c_eval_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_explnlp2_c_eval_name
			   ,&RTOp_TOp_explnlp2_c_eval_vtbl
			   ))
			assert(0);
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace NLPInterfacePack {

ExampleNLPObjGradient::ExampleNLPObjGradient(
	const VectorSpace::space_ptr_t&  vec_space
	,value_type                      xo
	,bool                            has_bounds
	,bool                            dep_bounded
	)
	:vec_space_(vec_space), vec_space_comp_(NULL,0)
	,initialized_(false), obj_scale_(1.0)
	,has_bounds_(has_bounds), force_xinit_in_bounds_(true), n_(2*vec_space->dim())
{
	namespace rcp = MemMngPack;

	// Assert the size of the NLP
	THROW_EXCEPTION(
		vec_space->dim() <= 0, std::logic_error
		,"ExampleNLPObjGradient::ExampleNLPObjGradient(...) Error!" );

	// Setup the aggregate vector space object
	BasisSystemCompositeStd::initialize_space_x(
		vec_space, vec_space, &var_dep_, &var_indep_, &vec_space_comp_ );

	// Set the initial starting point.
	xinit_ = vec_space_comp_->create_member();
	*xinit_ = xo;

	// Setup the sparse bounds
	//
	// xl(i) = 0.01  \ 
	//                }  for i <: bounded_rng
	// xu(i) = 20    /

	const size_type
		m = n_ / 2;

	xl_ = vec_space_comp_->create_member();
	xu_ = vec_space_comp_->create_member();

	if(has_bounds) {
		const Range1D
			bounded_rng   = ( dep_bounded ? var_dep_   : var_indep_ ),
			unbounded_rng = ( dep_bounded ? var_indep_ : var_dep_   );
		*xl_->sub_view(bounded_rng)   = 0.01;
		*xl_->sub_view(unbounded_rng) = -NLP::infinite_bound();
		*xu_->sub_view(bounded_rng)   = 20.0;
		*xu_->sub_view(unbounded_rng) = +NLP::infinite_bound();
	}
	else {
		*xl_ = -NLP::infinite_bound();
		*xu_ = +NLP::infinite_bound();
	}
}

// Overridden public members from NLP

void ExampleNLPObjGradient::initialize(bool test_setup)
{

#ifndef _WINDOWS
    using NLPInterfacePack::NLPFirstOrderDirect;
#endif

	if( initialized_ ) {
		NLPObjGradient::initialize(test_setup);
		return;
	}

	AbstractLinAlgPack::force_in_bounds( *xl_, *xu_, xinit_.get() );

	NLPObjGradient::initialize(test_setup);

	initialized_ = true;
}

bool ExampleNLPObjGradient::is_initialized() const
{
	return initialized_;
}

size_type ExampleNLPObjGradient::n() const
{
	assert_is_initialized();
	return n_;
}

size_type ExampleNLPObjGradient::m() const
{
	assert_is_initialized();
	return n_ / 2;
}

NLP::vec_space_ptr_t ExampleNLPObjGradient::space_x() const
{
	return vec_space_comp_;
}

NLP::vec_space_ptr_t ExampleNLPObjGradient::space_c() const
{
	return vec_space_;
}

NLP::vec_space_ptr_t ExampleNLPObjGradient::space_h() const
{
	return MemMngPack::null;
}

size_type ExampleNLPObjGradient::num_bounded_x() const
{
	return has_bounds_ ? n_/2 : 0;
}

void ExampleNLPObjGradient::force_xinit_in_bounds(bool force_xinit_in_bounds)
{
	force_xinit_in_bounds_ = force_xinit_in_bounds;
}

bool ExampleNLPObjGradient::force_xinit_in_bounds() const
{
	return force_xinit_in_bounds_;
}

const VectorWithOp& ExampleNLPObjGradient::xinit() const
{
	assert_is_initialized();
	return *xinit_;
}

const VectorWithOp& ExampleNLPObjGradient::xl() const
{
	assert_is_initialized();
	return *xl_;
}

const VectorWithOp& ExampleNLPObjGradient::xu() const
{
	assert_is_initialized();
	return *xu_;
}

value_type ExampleNLPObjGradient::max_var_bounds_viol() const
{
	return std::numeric_limits<value_type>::max(); // No limits on the bounds
}

const VectorWithOp& ExampleNLPObjGradient::hl() const
{
	THROW_EXCEPTION( true, NoBounds, "ExampleNLPObjGradient::hl(), Error, default is for mI() == 0" );
	return xl(); // will never execute.
}

const VectorWithOp& ExampleNLPObjGradient::hu() const
{
	THROW_EXCEPTION( true, NoBounds, "ExampleNLPObjGradient::hl(), Error, default is for mI() == 0" );
	return xu(); // will never execute.
}

void ExampleNLPObjGradient::scale_f( value_type scale_f )
{
	assert_is_initialized();
	obj_scale_ = scale_f;
}

value_type ExampleNLPObjGradient::scale_f() const
{
	assert_is_initialized();
	return obj_scale_;
}

void ExampleNLPObjGradient::report_final_solution(
	const VectorWithOp&    x
	,const VectorWithOp*   lambda
	,const VectorWithOp*   lambdaI
	,const VectorWithOp*   nu
	,bool                  optimal
	) const
{
	assert_is_initialized();

	// Do what you want with the solution (or final values) here.
	// For this example we will just ignore it.
}

Range1D ExampleNLPObjGradient::var_dep() const
{
	return var_dep_;
}

Range1D ExampleNLPObjGradient::var_indep() const
{
	return var_indep_;
}

// Overridden protected members from NLP

void ExampleNLPObjGradient::imp_calc_f(const VectorWithOp& x, bool newx
	, const ZeroOrderInfo& zero_order_info) const
{
	using AbstractLinAlgPack::dot;
	assert_is_initialized();
	f(); // assert f is set
	THROW_EXCEPTION( n() != x.dim(), std::length_error, "ExampleNLPObjGradient::imp_calc_f(...)"  );
	// f(x) = (obj_scale/2) * sum( x(i)^2, for i = 1..n )
	*zero_order_info.f = obj_scale_ / 2.0 * dot(x,x);
}

void ExampleNLPObjGradient::imp_calc_c(const VectorWithOp& x, bool newx
	, const ZeroOrderInfo& zero_order_info) const
{
	assert_is_initialized();
	const size_type n = this->n(), m = n/2;
	THROW_EXCEPTION( n != x.dim(), std::length_error, "ExampleNLPObjGradient::imp_calc_c(...)"  );

	// c(x)(j) = x(j) * (x(m+j) -1) - 10 * x(m+j) = 0, for j = 1...m

	VectorWithOp::vec_ptr_t
		xD= x.sub_view(var_dep()),
		xI = x.sub_view(var_indep());

	const int num_vecs = 2;
	const VectorWithOp*
		vecs[num_vecs] = { xD.get(), xI.get() };
	zero_order_info.c->apply_transformation(
		explnlp2_c_eval_op, num_vecs, vecs, 0, NULL, RTOp_REDUCT_OBJ_NULL );

}

void ExampleNLPObjGradient::imp_calc_h(
	const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const
{
	assert(0); // Should never be called!
}

// Overridden protected members from NLPFirstOrderInfo

void ExampleNLPObjGradient::imp_calc_Gf(const VectorWithOp& x, bool newx
	, const ObjGradInfo& obj_grad_info) const
{
	assert_is_initialized();
	THROW_EXCEPTION( n() != x.dim(), std::length_error, "ExampleNLPObjGradient::imp_calc_Gf(...)"  );
	// Gf = obj_scale * x
	LinAlgOpPack::V_StV(obj_grad_info.Gf,obj_scale_,x);
}

}	// end namespace NLPInterfacePack
