// ///////////////////////////////////////////////////////////////
// TSFCoreSerialVectorSpaceFactoryDecl.hpp

#ifndef TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_DECL_HPP
#define TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_DECL_HPP

#include "TSFCoreVectorSpaceFactory.hpp"

namespace TSFCore {

///
/** Implementation of a vector-space factory for serial vector spaces.
 */
template<class Scalar>
class SerialVectorSpaceFactory : public VectorSpaceFactory<Scalar> {
public:

	///
	SerialVectorSpaceFactory() {};

	/** @name Overridden from VectorSpaceFactory */
	//@{

	///
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > createVecSpc(int dim) const;

	//@}
	
}; // end class SerialVectorSpaceFactory

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_SPACE_FACTORY_SERIAL_DECL_HPP
