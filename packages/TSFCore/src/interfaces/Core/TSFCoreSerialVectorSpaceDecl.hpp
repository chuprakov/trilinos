// ////////////////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorSpaceDecl.hpp

#ifndef TSFCORE_VECTOR_SPACE_SERIAL_DECL_HPP
#define TSFCORE_VECTOR_SPACE_SERIAL_DECL_HPP

#include "TSFCoreSerialVectorSpaceBase.hpp"

namespace TSFCore {

///
/** <tt>%VectorSpace</tt> Subclass for serial vectors and multi-vectors.
 *
 * The default copy constructor and assignment operators are allowed
 * since they have the correct semantics.
 */
template<class Scalar>
class SerialVectorSpace : public SerialVectorSpaceBase<Scalar> {
public:

	/** @name Constructors / initializers */
	//@{

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	SerialVectorSpace( int dim = 0 );

	///
	/** Initialize given the dimension of the vector space.
	 *
	 * @param  dim   [in] The dimension of the vector space.
	 */
	void initialize( int dim );

	//@}

	/** @name Overridden from VectorSpece */
	//@{

	/// Returns 0 if uninitialized
	Index dim() const;
	/// Returns a <tt>SerialVector</tt> object.
	MemMngPack::ref_count_ptr<Vector<Scalar> > createMember() const;
	///
	MemMngPack::ref_count_ptr< const VectorSpace<Scalar> > clone() const;

	//@}

private:

	int   dim_;

}; // end class SerialVectorSpace

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SPACE_SERIAL_DECL_HPP
