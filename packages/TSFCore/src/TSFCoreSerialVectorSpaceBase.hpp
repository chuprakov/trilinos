// ////////////////////////////////////////////////////////////////////////////
// SerialVectorSpaceBase.cpp

#include <assert.h>

#include "TSFCoreSerialVectorSpaceBaseDecl.hpp"

namespace TSFCore {

template<class Scalar>
bool SerialVectorSpaceBase<Scalar>::isCompatible(const VectorSpace<Scalar>& aVecSpc ) const
{
	return this->dim() == aVecSpc.dim();
}

} // end namespace TSFCore
