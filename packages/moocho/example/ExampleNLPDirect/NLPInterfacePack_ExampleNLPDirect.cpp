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
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "Range1D.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"
#include "AbstractFactoryStd.h"

namespace {

// Calculate py and/or D
static RTOpPack::RTOpC          explnlp2_calc_py_D_op;

// Simple class for an object that will initialize the RTOp_Server.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
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
	:ExampleNLPObjGradient(vec_space,xo,has_bounds,dep_bounded)
{
	namespace rcp = MemMngPack;

	// Create the factory object for D
	factory_D_ = rcp::rcp(new MemMngPack::AbstractFactoryStd<MatrixWithOp,MatrixSymDiagonalStd>());
	NLPFirstOrderDirect::set_factories(
		MemMngPack::rcp(
			new MemMngPack::AbstractFactoryStd<MatrixSymWithOp,MatrixSymDiagonalStd>())               // D'*D
		,MemMngPack::rcp(
			new MemMngPack::AbstractFactoryStd<MatrixSymWithOpNonsingular,MatrixSymDiagonalStd>())    // S
		);
}

// Overridden public members from NLP

void ExampleNLPFirstOrderDirect::initialize(bool test_setup)
{

	if( initialized_ ) {
		NLPFirstOrderDirect::initialize(test_setup);
		ExampleNLPObjGradient::initialize(test_setup);
		return;
	}

	NLPFirstOrderDirect::initialize(test_setup);
	ExampleNLPObjGradient::initialize(test_setup);

	initialized_ = true;
}

bool ExampleNLPFirstOrderDirect::is_initialized() const
{
	return initialized_;
}

NLP::vec_space_ptr_t ExampleNLPFirstOrderDirect::space_h() const
{
	return ExampleNLPObjGradient::space_h();
}

const VectorWithOp& ExampleNLPFirstOrderDirect::hl() const
{
	return ExampleNLPObjGradient::hl();
}

const VectorWithOp& ExampleNLPFirstOrderDirect::hu() const
{
	return ExampleNLPObjGradient::hl();
}

// Overridden public members from NLPFirstOrderDirect

Range1D ExampleNLPFirstOrderDirect::var_dep() const
{
	return ExampleNLPObjGradient::var_dep();
}

Range1D ExampleNLPFirstOrderDirect::var_indep() const
{
	return ExampleNLPObjGradient::var_indep();
}

const NLPFirstOrderDirect::mat_fcty_ptr_t
ExampleNLPFirstOrderDirect::factory_D() const
{
	return factory_D_;
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
	mat_fcty_ptr_t::element_type::obj_ptr_t  D_ptr;
	if( rGf && !D ) {
		D_ptr = this->factory_D()->create();
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
		D_diag->init_identity( *this->space_c(), 0.0 );
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

// Overriden protected members for NLP

void ExampleNLPFirstOrderDirect::imp_calc_h(
	const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const
{
	assert(0); // Should never be called!
}

}	// end namespace NLPInterfacePack
