// //////////////////////////////////////////
// ExampleNLPFirstOrderDirect.cpp
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

#include "ExampleNLPFirstOrderDirect.h"
#include "ExampleNLPFirstOrderDirectRTOps.h"
#include "AbstractLinAlgPack/include/BasisSystemCompositeStd.h"
#include "AbstractLinAlgPack/include/MatrixSpaceStd.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "Range1D.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace {

// Evaluate c(x)
static RTOpPack::RTOpC          explnlp2_c_eval_op;
// Calculate py and/or D
static RTOpPack::RTOpC          explnlp2_calc_py_D_op;

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
		// Operator calcualte py and/or D
		if(0!=RTOp_TOp_explnlp2_calc_py_D_construct( 0, &explnlp2_calc_py_D_op.op() ))
			assert(0);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_explnlp2_calc_py_D_name
			   ,&RTOp_TOp_explnlp2_calc_py_D_vtbl
			   ))
			assert(0);
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace NLPInterfacePack {

ExampleNLPFirstOrderDirect::ExampleNLPFirstOrderDirect(
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
		,"ExampleNLPFirstOrderDirect::ExampleNLPFirstOrderDirect(...) Error!" );

	// Setup the aggregate vector space object
	BasisSystemCompositeStd::initialize_space_x(
		vec_space, vec_space, &var_dep_, &var_indep_, &vec_space_comp_ );

	// Create the MatrixSpace object for D
	typedef MatrixSpaceStd<MatrixWithOp,MatrixSymDiagonalStd>  space_D_con_t;
	space_D_ = rcp::rcp_implicit_cast<mat_space_ptr_t::element_type>(
		rcp::ref_count_ptr< const MatrixSpaceStd<MatrixWithOp,MatrixSymDiagonalStd> >(
			new MatrixSpaceStd<MatrixWithOp,MatrixSymDiagonalStd>(vec_space,vec_space)
			)
		);
	
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
			bounded_rng   = ( dep_bounded 
							  ? Range1D(var_dep_.lbound()  ,var_dep_.lbound()-1+m)
							  : Range1D(var_indep_.lbound(),var_dep_.lbound()-1+m) ),
			unbounded_rng = ( dep_bounded
							  ? Range1D(var_indep_.lbound(),var_dep_.lbound()-1+m)
							  : Range1D(var_dep_.lbound()  ,var_dep_.lbound()-1+m) );
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

void ExampleNLPFirstOrderDirect::initialize()
{

#ifndef _WINDOWS
    using NLPInterfacePack::NLPFirstOrderDirect;
#endif

	if( initialized_ ) {
		NLPFirstOrderDirect::initialize();
		return;
	}

	AbstractLinAlgPack::force_in_bounds( *xl_, *xu_, xinit_.get() );

//	xinit_->set_ele(n_/3+1,1.0/0.0); // Uncomment to throw in an invalid value
//	xinit_->set_ele((2*n_)/3+1,0.0/0.0); // Uncomment to throw in an invalid value

	NLPFirstOrderDirect::initialize();

	initialized_ = true;
}

bool ExampleNLPFirstOrderDirect::is_initialized() const
{
	return initialized_;
}

size_type ExampleNLPFirstOrderDirect::n() const
{
	assert_is_initialized();
	return n_;
}

size_type ExampleNLPFirstOrderDirect::m() const
{
	assert_is_initialized();
	return n_ / 2;
}

NLP::vec_space_ptr_t ExampleNLPFirstOrderDirect::space_x() const
{
	return vec_space_comp_;
}

NLP::vec_space_ptr_t ExampleNLPFirstOrderDirect::space_c() const
{
	return vec_space_;
}

size_type ExampleNLPFirstOrderDirect::num_bounded_x() const
{
	return has_bounds_ ? n_/2 : 0;
}

void ExampleNLPFirstOrderDirect::force_xinit_in_bounds(bool force_xinit_in_bounds)
{
	force_xinit_in_bounds_ = force_xinit_in_bounds;
}

bool ExampleNLPFirstOrderDirect::force_xinit_in_bounds() const
{
	return force_xinit_in_bounds_;
}

const VectorWithOp& ExampleNLPFirstOrderDirect::xinit() const
{
	assert_is_initialized();
	return *xinit_;
}

const VectorWithOp& ExampleNLPFirstOrderDirect::xl() const
{
	assert_is_initialized();
	return *xl_;
}

const VectorWithOp& ExampleNLPFirstOrderDirect::xu() const
{
	assert_is_initialized();
	return *xu_;
}

value_type ExampleNLPFirstOrderDirect::max_var_bounds_viol() const
{
	return std::numeric_limits<value_type>::max(); // No limits on the bounds
}

void ExampleNLPFirstOrderDirect::scale_f( value_type scale_f )
{
	assert_is_initialized();
	obj_scale_ = scale_f;
}

value_type ExampleNLPFirstOrderDirect::scale_f() const
{
	assert_is_initialized();
	return obj_scale_;
}

void ExampleNLPFirstOrderDirect::report_final_solution(
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

// Overridden public members from NLPFirstOrderDirect

Range1D ExampleNLPFirstOrderDirect::var_dep() const
{
	return var_dep_;
}

Range1D ExampleNLPFirstOrderDirect::var_indep() const
{
	return var_indep_;
}

const NLPFirstOrderDirect::mat_space_ptr_t&
ExampleNLPFirstOrderDirect::space_D() const
{
	return space_D_;
}

void ExampleNLPFirstOrderDirect::calc_point(
	const VectorWithOp     &x
	,value_type            *f
	,VectorWithOpMutable   *c
	,bool                  recalc_c
	,VectorWithOpMutable   *h
	,VectorWithOpMutable   *Gf
	,VectorWithOpMutable   *py
	,VectorWithOpMutable   *rGf
	,MatrixWithOp          *GcU
	,MatrixWithOp          *Gh
	,MatrixWithOp          *D
	,MatrixWithOp          *V
	,MatrixWithOp          *P
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgOpPack::Vp_MtV;

	assert_is_initialized();

	const size_type
		n = this->n(),
		m = n/2;

	// Validate the input

#ifdef _DEBUG
	THROW_EXCEPTION(
		x.dim() != n, std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error x.dim() = " << x.dim()
		<< " != n = " << n );
	THROW_EXCEPTION(
		c && !this->space_c()->is_compatible(c->space()), std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error c is not compatible" );
	THROW_EXCEPTION(
		h, std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, there are no inequalities h(x)" );
	THROW_EXCEPTION(
		Gf && !this->space_x()->is_compatible(Gf->space()), std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, Gf is not compatible" );
	THROW_EXCEPTION(
		py && !this->space_x()->sub_space(this->var_dep())->is_compatible(py->space()), std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, py is not compatible" );
	THROW_EXCEPTION(
		rGf && !this->space_x()->sub_space(this->var_dep())->is_compatible(rGf->space()), std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, py is not compatible" );
	THROW_EXCEPTION(
		GcU, std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, there are no undecomposed equalities" );
	THROW_EXCEPTION(
		Gh, std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, there are no inequalities h(x)" );
	THROW_EXCEPTION(
		D && !dynamic_cast<MatrixSymDiagonalStd*>(D), std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, D is not compatible" );
	THROW_EXCEPTION(
		V, std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, there are no undecomposed equalities" );
	THROW_EXCEPTION(
		P, std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...), Error, there are no general inequalities h(x)" );
	THROW_EXCEPTION(
		py!=NULL && c==NULL, std::invalid_argument
		,"ExampleNLPFirstOrderDirect::calc_point(...) : "
		"Error, if py!=NULL then c!=NULL must also be true" );
#endif

	// ///////////////////////////////////
	// Compute f(x), c(x) and Gf(x)

	typedef ExampleNLPFirstOrderDirect  this_t;

	// Make temp Gf if needed
	VectorSpace::vec_mut_ptr_t  Gf_ptr;
	if( rGf && !Gf ) {
		Gf_ptr = this->space_x()->create_member();
		Gf = Gf_ptr.get();
	}

	// Make temp D if needed
	mat_space_ptr_t::element_type::mat_ptr_t  D_ptr;
	if( rGf && !D ) {
		D_ptr = this->space_D()->create_member();
		D = D_ptr.get();
	}

	// Remember what references are already set
	value_type             *f_saved  = const_cast<this_t*>(this)->get_f();
	VectorWithOpMutable    *c_saved  = const_cast<this_t*>(this)->get_c();
	VectorWithOpMutable    *Gf_saved = const_cast<this_t*>(this)->get_Gf();
	// Set and compute the quantities
	const_cast<this_t*>(this)->set_f(f);
	const_cast<this_t*>(this)->set_c(c);
	const_cast<this_t*>(this)->set_Gf(Gf);
	if(Gf)
		this->calc_Gf(x,true);
	if(c && recalc_c)
		this->calc_c(x,false);
	if(f)
		this->calc_f(x,false);
	// Reset the remembered references
	const_cast<this_t*>(this)->set_f(f_saved);
	const_cast<this_t*>(this)->set_c(c_saved);
	const_cast<this_t*>(this)->set_Gf(Gf_saved);

	// ////////////////////////////////////////////////////////////////////////
	// Compute py = -inv(C)*c and/or D = -inv(C)*N at the same time
	// 
	//				[ 1/(1-x(m+1))											]
	//				[				1/(1-x(m+2))							]
	//	-inv(C)	= 	[								.						]
	//				[									.					]
	//				[										1/(1-x(m+m))	]
	//
	//
	//				[ x(1) - 10												]
	//				[				x(2) - 10								]
	//	N 		= 	[								.						]
	//				[									.					]
	//				[										x(m) - 10		]

	MatrixSymDiagonalStd
		*D_diag = dynamic_cast<MatrixSymDiagonalStd*>(D);

	VectorWithOp::vec_ptr_t
		xD= x.sub_view(this->var_dep()),
		xI = x.sub_view(this->var_indep());

	int task;
	if     ( py  && !D )  task = 0;
	else if( !py &&  D )  task = 1;
	else if( py  &&  D )  task = 2;
	
	assert(0==RTOp_TOp_explnlp2_calc_py_D_set_task(task,&explnlp2_calc_py_D_op.op()));

	const int                    num_vecs = task < 2 ? 2 : 3;
	const VectorWithOp*          vecs[3] = { NULL, NULL, NULL };
	const int                    num_targ_vecs = task < 2 ? 0 : 1;
	VectorWithOpMutable*         targ_vec0 = NULL;
	VectorWithOpMutable*         targ_vec1 = NULL;

	// targ_vec0 will have apply_transformation(...) called on it.
	if(D) {
		D_diag->init_identity( *vec_space_, 0.0 );
		targ_vec0= &D_diag->diag();
	}
	else if(py)
		targ_vec0 = py;
	else
		assert(0); // Only local error?
	// targ_vec1 will be passed to apply_transformation(...)
	if(py && D)
		targ_vec1 = py;
	
	// vecs[...]
	int k = 0;
	vecs[k] = xD.get(); ++k;
	if(D)  { vecs[k] = xI.get(); ++k; }
	if(py) { vecs[k] = c;        ++k; }

	targ_vec0->apply_transformation(
		explnlp2_calc_py_D_op, num_vecs, vecs, num_targ_vecs, &targ_vec1
		,RTOp_REDUCT_OBJ_NULL );

	// rGf = Gf(var_indep) + D' * Gf(var_dep)
	if(rGf) {
		*rGf = *Gf->sub_view(this->var_indep());
		Vp_MtV( rGf, *D, BLAS_Cpp::trans, *Gf->sub_view(this->var_dep()) );
	}
}

void ExampleNLPFirstOrderDirect::calc_semi_newton_step(
	const VectorWithOp    &x
	,VectorWithOpMutable  *c
	,bool                 recalc_c
	,VectorWithOpMutable  *py
	) const
{
	// In this case just call calc_point(...).
	// In a more specialized application, this could be much cheaper!
	calc_point(x,NULL,c,recalc_c,NULL,NULL,py,NULL,NULL,NULL,NULL,NULL,NULL);
}

// Overridden protected members from NLP

void ExampleNLPFirstOrderDirect::imp_calc_f(const VectorWithOp& x, bool newx
	, const ZeroOrderInfo& zero_order_info) const
{
	using AbstractLinAlgPack::dot;
	assert_is_initialized();
	f(); // assert f is set
	THROW_EXCEPTION( n() != x.dim(), std::length_error, "ExampleNLPFirstOrderDirect::imp_calc_f(...)"  );
	// f(x) = (obj_scale/2) * sum( x(i)^2, for i = 1..n )
	*zero_order_info.f = obj_scale_ / 2.0 * dot(x,x);
}

void ExampleNLPFirstOrderDirect::imp_calc_c(const VectorWithOp& x, bool newx
	, const ZeroOrderInfo& zero_order_info) const
{
	assert_is_initialized();
	const size_type n = this->n(), m = n/2;
	THROW_EXCEPTION( n != x.dim(), std::length_error, "ExampleNLPFirstOrderDirect::imp_calc_c(...)"  );

	// c(x)(j) = x(j) * (x(m+j) -1) - 10 * x(m+j) = 0, for j = 1...m

	VectorWithOp::vec_ptr_t
		xD= x.sub_view(this->var_dep()),
		xI = x.sub_view(this->var_indep());

	const int num_vecs = 2;
	const VectorWithOp*
		vecs[num_vecs] = { xD.get(), xI.get() };
	zero_order_info.c->apply_transformation(
		explnlp2_c_eval_op, num_vecs, vecs, 0, NULL, RTOp_REDUCT_OBJ_NULL );

}

// Overridden protected members from NLPFirstOrderInfo

void ExampleNLPFirstOrderDirect::imp_calc_Gf(const VectorWithOp& x, bool newx
	, const ObjGradInfo& obj_grad_info) const
{
	assert_is_initialized();
	THROW_EXCEPTION( n() != x.dim(), std::length_error, "ExampleNLPFirstOrderDirect::imp_calc_Gf(...)"  );
	// Gf = obj_scale * x
	LinAlgOpPack::V_StV(obj_grad_info.Gf,obj_scale_,x);
}

}	// end namespace NLPInterfacePack
