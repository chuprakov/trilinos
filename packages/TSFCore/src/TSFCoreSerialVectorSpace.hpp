// ////////////////////////////////////////////////////////////////////////////
// SerialVectorSpace.hpp

#ifndef TSFCORE_VECTOR_SPACE_SERIAL_HPP
#define TSFCORE_VECTOR_SPACE_SERIAL_HPP

#include <assert.h>

#include "TSFCoreSerialVectorSpaceDecl.hpp"
#include "TSFCoreSerialVector.hpp"
#include "ThrowException.hpp"

namespace TSFCore {

template<class Scalar>
SerialVectorSpace<Scalar>::SerialVectorSpace( int dim )
{
	initialize(dim);
}

template<class Scalar>
void SerialVectorSpace<Scalar>::initialize( int dim )
{
	dim_ = dim;
}

// Overridden from VectorSpace

template<class Scalar>
Index SerialVectorSpace<Scalar>::dim() const
{
	return dim_;
}

template<class Scalar>
MemMngPack::ref_count_ptr<Vector<Scalar> >
SerialVectorSpace<Scalar>::createMember() const
{
	return MemMngPack::rcp(new SerialVector<Scalar>(dim_));
}

template<class Scalar>
MemMngPack::ref_count_ptr< const VectorSpace<Scalar> >
SerialVectorSpace<Scalar>::clone() const
{
	return MemMngPack::rcp(new SerialVectorSpace<Scalar>(*this));
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_SERIAL_HPP
