// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// /////////////////////////////////////////////////////////////
// TSFCoreVector.hpp

#ifndef TSFCORE_VECTOR_HPP
#define TSFCORE_VECTOR_HPP

#include "TSFCoreVectorDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "RTOpPack_ROpGetSubVector.hpp"
#include "RTOpPack_TOpSetSubVector.hpp"
//#include "RTOp_TOp_set_sub_vector.h"
//#include "RTOpPack_RTOpC.hpp"
#include "Teuchos_TestForException.hpp"

// RAB: 3/24/2004: Write direct templated RTOpT<> subclasses for get_sub_vector and set_sub_vector
// operators!

namespace TSFCore {

template<class Scalar>
void Vector<Scalar>::getSubVector( const Range1D& rng_in, RTOpPack::SubVectorT<Scalar>* sub_vec_inout ) const
{
  using DynamicCastHelperPack::dyn_cast;
	const Range1D rng = rng_in.full_range() ? Range1D(1,this->space()->dim()) : rng_in;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		this->space()->dim() < rng.ubound(), std::out_of_range
		,"Vector<Scalar>::getSubVector(rng,...): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
		<<"] is not in range = [1,"<<this->space()->dim()<<"]" );
#endif
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
	if( sub_vec_inout->values() ) {
		free( (void*)sub_vec_inout->values()  );
	}
	// Initialize the operator
	RTOpPack::ROpGetSubVector<Scalar> get_sub_vector_op(rng.lbound(),rng.ubound());
	// Create the reduction object (another sub_vec)
  Teuchos::RefCountPtr<RTOpPack::ReductTarget>
    reduct_obj = get_sub_vector_op.reduct_obj_create(); // This is really of type RTOpPack::SubVectorT<Scalar>!
	// Perform the reduction (get the sub-vector requested)
	const Vector* sub_vecs[] = { this };
	applyOp(
		get_sub_vector_op,1,sub_vecs,0,NULL,&*reduct_obj
		,rng.lbound(),rng.size(),rng.lbound()-1 // first_ele, sub_dim, global_offset
		);
	// Get the sub-vector.  Note reduct_obj will go out of scope so the sub_vec parameter will
	// own the memory allocated within get_sub_vector_op.create_reduct_obj_raw(...).  This is okay
	// since the client is required to call release_sub_vector(...) so release memory!
  RTOpPack::ReductTargetSubVectorT<Scalar> &sub_vec_ro = get_sub_vector_op(*reduct_obj);
  sub_vec_ro.transfer(sub_vec_inout);
}

template<class Scalar>
void Vector<Scalar>::freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
  RTOpPack::ReductTargetSubVectorT<Scalar>::free(sub_vec);
}

template<class Scalar>
void Vector<Scalar>::getSubVector( const Range1D& rng, RTOpPack::MutableSubVectorT<Scalar>* sub_vec_inout )
{
	//
	// Here we get a copy of the data for the sub-vector that the
	// client will modify.  We must later commit these changes to the
	// actual vector when the client calls commitSubVector(...).
	// Note, this implementation is very dependent on the behavior of
	// the default implementation of constant version of
	// Vector<Scalar>::getSubVector(...) and the implementation of
	// Vector<Scalar>::setSubVector(...)!
	//
	RTOpPack::SubVectorT<Scalar> sub_vec;
	Vector<Scalar>::getSubVector( rng, &sub_vec );
	sub_vec_inout->initialize(
		sub_vec.globalOffset(),sub_vec.subDim(),const_cast<Scalar*>(sub_vec.values()),sub_vec.stride());
}

template<class Scalar>
void Vector<Scalar>::commitSubVector( RTOpPack::MutableSubVectorT<Scalar>* sub_vec_inout )
{
	RTOpPack::SparseSubVectorT<Scalar> spc_sub_vec(
		sub_vec_inout->globalOffset(), sub_vec_inout->subDim()
		,sub_vec_inout->values(), sub_vec_inout->stride()
		);
	Vector<Scalar>::setSubVector( spc_sub_vec );           // Commit the changes!
	RTOpPack::SubVectorT<Scalar> sub_vec(*sub_vec_inout);
	Vector<Scalar>::freeSubVector( &sub_vec );             // Free the memory!
	sub_vec_inout->set_uninitialized();                    // Make null as promised!
}

template<class Scalar>
void Vector<Scalar>::setSubVector( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec )
{
/*
	RTOp_SparseSubVector spc_sub_vec;
	RTOp_sparse_sub_vector(
		sub_vec.globalOffset(), sub_vec.subDim(), sub_vec.subNz()
		,sub_vec.values(), sub_vec.valuesStride(), sub_vec.indices(), sub_vec.indicesStride()
		,sub_vec.localOffset(), sub_vec.isSorted()
		,&spc_sub_vec
		);
	RTOpPack::RTOpC  set_sub_vector_op;
  TEST_FOR_EXCEPTION(
    0!=RTOp_TOp_set_sub_vector_construct( &spc_sub_vec, &set_sub_vector_op.op() )
    ,std::invalid_argument,"Error!"
    );
*/
  RTOpPack::TOpSetSubVector<Scalar> set_sub_vector_op(sub_vec);
	Vector* targ_vecs[1] = { this };
	this->applyOp(
		set_sub_vector_op,0,NULL,1,targ_vecs,NULL
		,sub_vec.globalOffset()+1,sub_vec.subDim(),sub_vec.globalOffset() // first_ele, sub_dim, global_offset
		);
}

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_HPP
