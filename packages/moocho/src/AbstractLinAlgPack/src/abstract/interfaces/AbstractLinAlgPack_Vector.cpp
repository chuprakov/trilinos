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
#include "RTOp_ROp_dot_prod.h"
#include "RTOp_ROp_sum.h"
#include "RTOp_ROp_norms.h"
#include "RTOp_ROp_num_nonzeros.h"
#include "RTOp_ROp_get_sub_vector.h"
#include "RTOpCppC.hpp"
#include "print_sub_vector.hpp"
#include "Range1D.hpp"
#include "dynamic_cast_verbose.hpp"
#include "Teuchos_TestForException.hpp"

namespace {

// get element operator
static RTOpPack::RTOpC          sum_op;
static RTOpPack::ReductTarget   sum_targ;
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
		if(0>RTOp_ROp_sum_construct( &sum_op.op() ))
			assert(0);
		sum_op.reduct_obj_create(&sum_targ);
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
	RTOpPack::SubVector sub_vec;
	this->get_sub_vector( Range1D(), &sub_vec );
	sub_vec.initialize( sub_vec.globalOffset() + global_offset, sub_vec.subDim(), sub_vec.values(), sub_vec.stride() );
	RTOpPack::output(out,sub_vec,print_dim,newline);
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
	sum_targ.reinit();
	const Vector *vecs[1] = { this };
	AbstractLinAlgPack::apply_op(
		sum_op,1,vecs,0,NULL,sum_targ.obj()
		,i,1,0 // first_ele, sub_dim, global_offset
		);
	return RTOp_ROp_sum_val(sum_targ.obj());
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
	TEST_FOR_EXCEPTION(
		rng.ubound() > dim, std::out_of_range
		,"Vector::sub_view(rng): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1,"<<dim<<"]" );
#endif	
	if( rng.lbound() == 1 && rng.ubound() == dim )
		return vec_ptr_t( this, false );
	return Teuchos::rcp(
		new VectorSubView(
			vec_ptr_t( this, false )
			,rng ) );
}

void Vector::get_sub_vector( const Range1D& rng_in, RTOpPack::SubVector* sub_vec_inout ) const
{
	const Range1D rng = rng_in.full_range() ? Range1D(1,this->space().dim()) : rng_in;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		this->space().dim() < rng.ubound(), std::out_of_range
		,"Vector::get_sub_vector(rng,...): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
		<<"] is not in range = [1,"<<this->space().dim()<<"]" );
#endif
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
	if( sub_vec_inout->values() ) {
		free( (void*)sub_vec_inout->values()  );
	}
	// Initialize the operator
	RTOpPack::RTOpC get_sub_vector_op;
	if(0>RTOp_ROp_get_sub_vector_construct( rng.lbound(), rng.ubound(),&get_sub_vector_op.op()))
		assert(0);
	// Create the reduction object (another sub_vec)
	RTOp_ReductTarget reduct_obj = RTOp_REDUCT_OBJ_NULL;
	get_sub_vector_op.reduct_obj_create_raw(&reduct_obj); // This is really of type RTOpPack::SubVectorT<Scalar>!
	// Perform the reduction (get the sub-vector requested)
	const size_t  num_vecs = 1;
	const Vector* sub_vecs[num_vecs] = { this };
	AbstractLinAlgPack::apply_op(
		get_sub_vector_op,num_vecs,sub_vecs,0,NULL,reduct_obj
		,rng.lbound(),rng.size(),rng.lbound()-1 // first_ele, sub_dim, global_offset
		);
	// Set the sub-vector.  Note reduct_obj will go out of scope so the sub_vec parameter will
	// own the memory allocated within get_sub_vector_op.create_reduct_obj_raw(...).  This is okay
	// since the client is required to call release_sub_vector(...) so release memory!
	RTOp_SubVector sub_vec = RTOp_ROp_get_sub_vector_val(reduct_obj);
	sub_vec_inout->initialize(sub_vec.global_offset,sub_vec.sub_dim,sub_vec.values,sub_vec.values_stride);
	free(reduct_obj); // Now *sub_vec owns the values[] and indices[] arrays!
}

void Vector::free_sub_vector( RTOpPack::SubVector* sub_vec ) const
{
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
	if( sub_vec->values() )
		free( (void*)sub_vec->values() );
	sub_vec->set_uninitialized();
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
