// /////////////////////////////////////////////////////////////////
// SerialVector.hpp

#ifndef TSFCORE_VECTOR_SERIAL_HPP
#define TSFCORE_VECTOR_SERIAL_HPP

#include <stdexcept>

#include "TSFCoreSerialVectorDecl.hpp"
#include "ThrowException.hpp"

namespace TSFCore {

template<class Scalar>
SerialVector<Scalar>::~SerialVector()
{
	free_mem();
}

template<class Scalar>
SerialVector<Scalar>::SerialVector(
	int dim
	)
	:v_(NULL),vs_(0),ownsMem_(false),space_serial_(dim)
{
	this->initialize(dim);
}

template<class Scalar>
SerialVector<Scalar>::SerialVector(
	Scalar  v[]
	,int    vs
	,int    dim
	,bool   ownsMem
	)
	:v_(NULL),vs_(0),ownsMem_(false),space_serial_(dim)
{
	this->initialize(v,vs,dim,ownsMem);
}

template<class Scalar>
void SerialVector<Scalar>::initialize(
	int dim
	)
{
	namespace mmp = MemMngPack;
	free_mem();
	this->initialize(
		new Scalar[dim]
		,1
		,dim
		,true
		);
	space_serial_.initialize(dim);
}

template<class Scalar>
void SerialVector<Scalar>::initialize(
	Scalar  v[]
	,int    vs
	,int    dim
	,bool   ownsMem
	)
{
	free_mem();
	v_       = v;
	vs_      = vs;
	dim_     = dim;
	ownsMem_ = ownsMem;
}

// Overridden from Vector

template<class Scalar>
MemMngPack::ref_count_ptr< const VectorSpace<Scalar> >
SerialVector<Scalar>::space() const
{
	return MemMngPack::rcp(&space_serial_,false);
}

template<class Scalar>
void SerialVector<Scalar>::getSubVector( const Range1D& rng_in, RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
	const Index      this_dim = dim_;
	const Range1D    rng      = RangePack::full_range(rng_in,1,this_dim);
	THROW_EXCEPTION(
		rng.ubound() > this_dim, std::out_of_range
		,"SerialVector<Scalar>::getSubVector(...) : Error, "
		"rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1," << this_dim << "]!" );
	sub_vec->initialize(
		rng.lbound()-1            // globalOffset
		,rng.size()               // subDim
		,v_+vs_*(rng.lbound()-1)  // values
		,vs_                      // stride
		);
}

template<class Scalar>
void SerialVector<Scalar>::freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
	sub_vec->set_uninitialized();  // Nothing to deallocate!
}

void SerialVector<Scalar>::getSubVector( const Range1D& rng_in, RTOpPack::MutableSubVectorT<Scalar>* sub_vec )
{
	const Index     this_dim = dim_;
	const Range1D   rng = RangePack::full_range(rng_in,1,this_dim);
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.ubound() > this_dim, std::out_of_range
		,"SerialVector<Scalar>::getSubVector(...) : Error, "
		"rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1," << this_dim << "]!" );
#endif
	sub_vec->initialize(
		rng.lbound()-1           // globalOffset
		,rng.size()              // subDim
		,v_+vs_*(rng.lbound()-1) // values
		,vs_                     // stride
		);
}

template<class Scalar>
void SerialVector<Scalar>::commitSubVector( RTOpPack::MutableSubVectorT<Scalar>* sub_vec )
{
	sub_vec->set_uninitialized();  // Nothing to deallocate!
}

template<class Scalar>
void SerialVector<Scalar>::setSubVector( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec )
{
	Vector<Scalar>::setSubVector(sub_vec);  // This implementation is okay
}

// private

template<class Scalar>
void SerialVector<Scalar>::free_mem()
{
	if(ownsMem_) delete [] v_;
	v_ = NULL;
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SERIAL_HPP
