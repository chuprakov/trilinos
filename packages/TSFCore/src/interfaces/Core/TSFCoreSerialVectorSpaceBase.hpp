// ////////////////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorSpaceBase.hpp

#ifndef TSFCORE_SERIAL_VECTOR_SPACE_BASE_HPP
#define TSFCORE_SERIAL_VECTOR_SPACE_BASE_HPP

#include "TSFCoreSerialVectorSpaceBaseDecl.hpp"

namespace TSFCore {

template<class Scalar>
bool SerialVectorSpaceBase<Scalar>::isInCore() const
{
	return true;
}

template<class Scalar>
bool SerialVectorSpaceBase<Scalar>::isCompatible(const VectorSpace<Scalar>& aVecSpc ) const
{
	return this->dim() == aVecSpc.dim();
}

} // end namespace TSFCore

#endif // TSFCORE_SERIAL_VECTOR_SPACE_BASE_HPP
