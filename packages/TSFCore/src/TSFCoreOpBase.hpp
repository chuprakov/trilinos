// //////////////////////////////////////////////////////////////////////////////////
// OpBase.hpp

#ifndef TSFCORE_OP_BASE_HPP
#define TSFCORE_OP_BASE_HPP

#include <assert.h>

#include "TSFCoreOpBaseDecl.hpp"

namespace TSFCore {

// Virtual functions with default implemenations

template<class Scalar>
bool OpBase<Scalar>::opSupported(ETransp) const
{
	return true;
}
}	// end namespace TSFCore

#endif // TSFCORE_OP_BASE_HPP
