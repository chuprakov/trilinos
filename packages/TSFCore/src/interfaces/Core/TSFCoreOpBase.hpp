// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreOpBase.hpp

#ifndef TSFCORE_OP_BASE_HPP
#define TSFCORE_OP_BASE_HPP

#include <assert.h>

#include "TSFCoreOpBaseDecl.hpp"

namespace TSFCore {

template<class Scalar>
OpBase<Scalar>::~OpBase()
{}

template<class Scalar>
bool OpBase<Scalar>::opSupported(ETransp) const
{
	return true;
}

}	// end namespace TSFCore

#endif // TSFCORE_OP_BASE_HPP
