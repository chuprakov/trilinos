// /////////////////////////////////////////////////////////////
// TSFCoreVector.hpp

#ifndef TSFCORE_VECTOR_HPP
#define TSFCORE_VECTOR_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreVectorDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "RTOp_ROp_get_sub_vector.h"
#include "RTOp_TOp_set_sub_vector.h"
#include "RTOpCppC.hpp"
#include "ThrowException.hpp"

namespace TSFCore {

template<class Scalar>
void Vector<Scalar>::getSubVector( const Range1D& rng_in, RTOpPack::SubVectorT<Scalar>* sub_vec_inout ) const
{
	const Range1D rng = rng_in.full_range() ? Range1D(1,this->space()->dim()) : rng_in;
#ifdef _DEBUG
	THROW_EXCEPTION(
		this->space()->dim() < rng.ubound(), std::out_of_range
		,"Vector<Scalar>::getSubVector(rng,...): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
		<<"] is not in range = [1,"<<this->space()->dim()<<"]" );
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
	applyOp(
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

template<class Scalar>
void Vector<Scalar>::freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
	// Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
	if( sub_vec->values() )
		free( (void*)sub_vec->values() );
	sub_vec->set_uninitialized();
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
	RTOp_SparseSubVector spc_sub_vec;
	RTOp_sparse_sub_vector(
		sub_vec.globalOffset(), sub_vec.subDim(), sub_vec.subNz()
		,sub_vec.values(), sub_vec.valuesStride(), sub_vec.indices(), sub_vec.indicesStride()
		,sub_vec.localOffset(), sub_vec.isSorted()
		,&spc_sub_vec
		);
	RTOpPack::RTOpC  set_sub_vector_op;
	if(0>RTOp_TOp_set_sub_vector_construct( &spc_sub_vec, &set_sub_vector_op.op() ))
		assert(0);
	Vector* targ_vecs[1] = { this };
	this->applyOp(
		set_sub_vector_op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL
		,sub_vec.globalOffset()+1,sub_vec.subDim(),sub_vec.globalOffset() // first_ele, sub_dim, global_offset
		);
}

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_HPP
