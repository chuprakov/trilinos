// //////////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorSpaceFactory.hpp

#ifndef TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_HPP
#define TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreSerialVectorSpaceFactoryDecl.hpp"
#include "TSFCoreSerialVectorSpace.hpp"

namespace TSFCore {

template<class Scalar>
MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >
SerialVectorSpaceFactory<Scalar>::createVecSpc(int dim) const
{
	return MemMngPack::rcp(new SerialVectorSpace<Scalar>(dim));
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_HPP
