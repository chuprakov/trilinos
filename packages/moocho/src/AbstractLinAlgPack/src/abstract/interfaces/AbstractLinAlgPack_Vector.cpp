// ///////////////////////////////////////////////////////////////////
// Vector.cpp
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

#include <limits>
#include <ostream>

#include "VectorMutable.hpp"
#include "VectorSubView.hpp"
#include "RTOpStdOpsLib/src/RTOp_ROp_dot_prod.h"
#include "RTOpStdOpsLib/src/RTOp_ROp_get_ele.h"
#include "RTOpStdOpsLib/src/RTOp_ROp_norms.h"
#include "RTOpStdOpsLib/src/RTOp_ROp_num_nonzeros.h"
#include "RTOpStdOpsLib/src/RTOp_ROp_get_sub_vector.h"
#include "RTOpPack/src/RTOpCppC.hpp"
#include "RTOpPack/src/print_sub_vector.hpp"
#include "Range1D.hpp"
#include "dynamic_cast_verbose.hpp"
#include "ThrowException.hpp"

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
		if(0>RTOp_ROp_get_ele_construct( 0, &get_ele_op.op() ))
			assert(0);
		get_ele_op.reduct_obj_create(&get_ele_targ);
		// Operator and target for norm 1
		if(0>RTOp_ROp_num_nonzeros_construct(&num_nonzeros_op.op() ))
			assert(0);
		num_nonzeros_op.reduct_obj_create(&num_nonzeros_targ);
		// Operator and target for norm 1
		if(0>RTOp_ROp_norm_1_construct(&norm_1_op.op() ))
			assert(0);
		norm_1_op.reduct_obj_create(&norm_1_targ);
		// Operator and target for norm 1
		if(0>RTOp_ROp_norm_2_construct(&norm_2_op.op() ))
			assert(0);
		norm_2_op.reduct_obj_create(&norm_2_targ);
		// Operator and target for norm 1
		if(0>RTOp_ROp_norm_inf_construct(&norm_inf_op.op() ))
			assert(0);
		norm_inf_op.reduct_obj_create(&norm_inf_targ);
		// Dot product operator and target
		if(0>RTOp_ROp_dot_prod_construct(&dot_op.op()))
			assert(0);
		dot_op.reduct_obj_create(&dot_targ);
		// Get sub-vector operator
		if(0>RTOp_ROp_get_sub_vector_construct(1,1,&get_sub_vector_op.op()))
			assert(0);
	}
}; 

// When the program starts, this object will be created
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace AbstractLinAlgPack {

Vector::Vector()
	:num_nonzeros_(-1), norm_1_(-1.0), norm_2_(-1.0), norm_inf_(-1.0) // uninitalized
{}

index_type Vector::dim() const
{
	return this->space().dim();
}

index_type Vector::nz() const
{
	if( num_nonzeros_ < 0 ) {
		num_nonzeros_targ.reinit();
		const Vector *vecs[1] = { this };
		AbstractLinAlgPack::apply_op(num_nonzeros_op,1,vecs,0,NULL,num_nonzeros_targ.obj());
		num_nonzeros_ = RTOp_ROp_num_nonzeros_val(num_nonzeros_targ.obj());
	}
	return num_nonzeros_;
}

std::ostream& Vector::output(
	std::ostream& out, bool print_dim , bool newline
	,index_type global_offset
	) const
{
	RTOp_SubVector sub_vec;
	RTOp_sub_vector_null(&sub_vec);
	this->get_sub_vector( Range1D(), &sub_vec );
	RTOp_SubVector  sub_vec_print = sub_vec;
	sub_vec_print.global_offset += global_offset;
	RTOpPack::output(out,sub_vec_print,print_dim,newline);
	this->free_sub_vector( &sub_vec );
	return out;
}

VectorMutable::vec_mut_ptr_t Vector::clone() const
{
	vec_mut_ptr_t
		vec = this->space().create_member();
	*vec = *this;
	return vec;
}

value_type Vector::get_ele(index_type i) const {
	assert(0==RTOp_ROp_get_ele_set_i( i, &get_ele_op.op() ));
	get_ele_targ.reinit();
	const Vector *vecs[1] = { this };
	AbstractLinAlgPack::apply_op(
		get_ele_op,1,vecs,0,NULL,get_ele_targ.obj()
		,i,1,i-1 // first_ele, sub_dim, global_offset
		);
	return RTOp_ROp_get_ele_val(get_ele_targ.obj());
}

value_type Vector::norm_1() const {
	if( norm_1_ < 0.0 ) {
		norm_1_targ.reinit();
		const Vector *vecs[1] = { this };
		AbstractLinAlgPack::apply_op(norm_1_op,1,vecs,0,NULL,norm_1_targ.obj());
		norm_1_ = RTOp_ROp_norm_1_val(norm_1_targ.obj());
	}
	return norm_1_;
}

value_type Vector::norm_2() const {
	if( 1 /*norm_2_ < 0.0*/ ) {
		norm_2_targ.reinit();
		const Vector *vecs[1] = { this };
		AbstractLinAlgPack::apply_op(norm_2_op,1,vecs,0,NULL,norm_2_targ.obj());
		norm_2_ = RTOp_ROp_norm_2_val(norm_2_targ.obj());
	}
	return norm_2_;
}

value_type Vector::norm_inf() const {
	if( 1 /*norm_inf_ < 0.0*/ ) {
		norm_inf_targ.reinit();
		const Vector *vecs[1] = { this };
		AbstractLinAlgPack::apply_op(norm_inf_op,1,vecs,0,NULL,norm_inf_targ.obj());
		norm_inf_ = RTOp_ROp_norm_inf_val(norm_inf_targ.obj());
	}
	return norm_inf_;
}

value_type Vector::inner_product(  const Vector& v ) const
{
	return this->space().inner_prod()->inner_prod(*this,v);
}

Vector::vec_ptr_t
Vector::sub_view( const Range1D& rng_in ) const
{
	namespace rcp = MemMngPack;
	const index_type dim = this->dim();
	const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.ubound() > dim, std::out_of_range
		,"Vector::sub_view(rng): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1,"<<dim<<"]" );
#endif	
	if( rng.lbound() == 1 && rng.ubound() == dim )
		return vec_ptr_t( this, false );
	return rcp::rcp(
		new VectorSubView(
			vec_ptr_t( this, false )
			,rng ) );
}

void Vector::get_sub_vector( const Range1D& rng_in, RTOp_SubVector* sub_vec	) const
{
	const Range1D rng = rng_in.full_range() ? Range1D(1,this->dim()) : rng_in;
#ifdef _DEBUG
	THROW_EXCEPTION(
		this->dim() < rng.ubound(), std::out_of_range
		,"Vector::get_sub_vector(rng,...): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
		<<"] is not in range = [1,"<<this->dim()<<"]" );
#endif
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
	if( sub_vec->values ) {
		free( (void*)sub_vec->values  );
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
	const Vector *vecs[1] = { this };
	AbstractLinAlgPack::apply_op(
		get_sub_vector_op,1,vecs,0,NULL,reduct_obj
		,rng.lbound(),rng.size(),rng.lbound()-1 // first_ele, sub_dim, global_offset
		);
	// Set the sub-vector.  Note reduct_obj will go out of scope so the sub_vec parameter will
	// own the memory allocated within get_sub_vector_op.create_reduct_obj_raw(...).  This is okay
	//  since the client is required to call release_sub_vector(...) so release memory!

	*sub_vec = RTOp_ROp_get_sub_vector_val(reduct_obj);
	free(reduct_obj); // Now *sub_vec owns the values[] and indices[] arrays!
}

void Vector::free_sub_vector( RTOp_SubVector* sub_vec ) const
{
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
	if( sub_vec->values )
		free( (void*)sub_vec->values  );
	RTOp_sub_vector_null( sub_vec );
}

void Vector::has_changed() const
{
	num_nonzeros_= -1;  // uninitalized;
	norm_1_ = norm_2_ = norm_inf_ = -1.0;
}

// protected

void Vector::finalize_apply_op(
	const size_t num_targ_vecs, VectorMutable** targ_vecs
	) const
{
	for( int k = 0; k < num_targ_vecs; ++k )
		targ_vecs[k]->has_changed();
}

} // end namespace AbstractLinAlgPack

// nonmember functions

void AbstractLinAlgPack::apply_op(
	const RTOpPack::RTOp       &op
	,const size_t              num_vecs
	,const Vector*             vecs[]
	,const size_t              num_targ_vecs
	,VectorMutable*            targ_vecs[]
	,RTOp_ReductTarget         reduct_obj
	,const index_type          first_ele
	,const index_type          sub_dim
	,const index_type          global_offset
	)
{
	if(num_vecs)
		vecs[0]->apply_op(op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset);
	else if (num_targ_vecs)
		targ_vecs[0]->apply_op(op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset);
}
