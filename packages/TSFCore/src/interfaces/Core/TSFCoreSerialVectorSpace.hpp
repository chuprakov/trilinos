// ////////////////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorSpace.hpp

#ifndef TSFCORE_VECTOR_SPACE_SERIAL_HPP
#define TSFCORE_VECTOR_SPACE_SERIAL_HPP

#include "TSFCoreSerialVectorSpaceDecl.hpp"
#include "TSFCoreSerialVectorSpaceBase.hpp"
#include "TSFCoreSerialVector.hpp"
#include "Teuchos_TestForException.hpp"

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
Teuchos::RefCountPtr<Vector<Scalar> >
SerialVectorSpace<Scalar>::createMember() const
{
	return Teuchos::rcp(new SerialVector<Scalar>(dim_));
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
SerialVectorSpace<Scalar>::clone() const
{
	return Teuchos::rcp(new SerialVectorSpace<Scalar>(*this));
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_SERIAL_HPP
