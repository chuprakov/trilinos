// //////////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorSpaceFactory.hpp

#ifndef TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_HPP
#define TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_HPP

#include "TSFCoreSerialVectorSpaceFactoryDecl.hpp"
#include "TSFCoreSerialVectorSpace.hpp"

namespace TSFCore {

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpace<Scalar> >
SerialVectorSpaceFactory<Scalar>::createVecSpc(int dim) const
{
	return Teuchos::rcp(new SerialVectorSpace<Scalar>(dim));
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_HPP
