// ////////////////////////////////////////////////////////////////////////////
// SerialVectorSpaceBaseDecl.hpp

#ifndef TSFCORE_VECTOR_SPACE_SERIAL_BASE_DECL_HPP
#define TSFCORE_VECTOR_SPACE_SERIAL_BASE_DECL_HPP

#include "TSFCoreVectorSpaceDecl.hpp"

namespace TSFCore {

///
/** <tt>%VectorSpace</tt> node subclass for serial vectors and multi-vectors.
 *
 * All a concrete subclass must do is to override the <tt>createMember()</tt>
 * method.
 */
template<class Scalar>
class SerialVectorSpaceBase : public VectorSpace<Scalar> {
public:

	/** @name Overridden from VectorSpece */
	//@{

	///
	/** Returns true if <tt>vecSpc.dim() == this->dim()</tt>.
	 *
	 * The assumption here is that <tt>Vector::getSubVector()</tt>,
	 * <tt>Vector::freeSubVector()</tt> and <tt>Vector::commitSubVector()</tt>
	 * can be used to implement all of the methods on an SMP machine in an
	 * efficient manner.
	 */
 	bool isCompatible(const VectorSpace<Scalar>& vecSpc) const;

	//@}

}; // end class SerialVectorSpaceBase

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_SERIAL_BASE_DECL_HPP
