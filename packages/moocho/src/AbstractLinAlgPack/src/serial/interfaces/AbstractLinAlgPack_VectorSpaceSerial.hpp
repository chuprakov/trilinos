// ////////////////////////////////////////////////////////////////////////////
// VectorSpaceSerial.h

#ifndef VECTOR_SPACE_SERIAL_H
#define VECTOR_SPACE_SERIAL_H

#include "SparseLinAlgPack/include/SparseLinAlgPackTypes.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"

namespace SparseLinAlgPack {

///
/** Subclass for vector space objects that create <tt>LinAlgPack::Vector</tt> vector objects.
 *
 * The default constructor, copy constructor and assignment operators
 * are allowed since they have the correct semantics.
 */
class VectorSpaceSerial : public VectorSpace {
public:

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorSpaceSerial( size_type dim );

	///
	/** Initialize given the dimension of the vector space.
	 *
	 * @param  dim   [in] The dimension of the vector space.
	 */
	void initialize( size_type dim );

	/** @name Overridden from VectorSpece */
	//@{

	///
	/** Returns true if <tt>vec_space.dim() == this->dim()</tt>.
	 *
	 * The assumption here is that since <tt>VectorWithOp::get_sub_vector()</tt>
	 * and <tt>VectorWithOpMutable::get_sub_vector()</tt> can be used to implement
	 * all of the methods on an SMP machine.
	 */
 	bool is_compatible(const VectorSpace& vec_space) const;
	/// Returns 0 if uninitialized
	index_type dim() const;
	///
	space_ptr_t clone() const;
	///
	vec_mut_ptr_t create_member() const;

	//@}

private:

	size_type     dim_;

}; // end class VectorSpaceSerial

} // end namespace SparseLinAlgPack

#endif // VECTOR_SPACE_SERIAL_H
