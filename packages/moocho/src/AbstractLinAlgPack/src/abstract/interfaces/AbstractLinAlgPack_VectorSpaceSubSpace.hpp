// /////////////////////////////////////////////////////////////
// VectorSpaceSubSpace.h

#ifndef VECTOR_SPACE_SUB_SPACE_H
#define VECTOR_SPACE_SUB_SPACE_H

#include "VectorSpace.h"
#include "Range1D.h"

namespace AbstractLinAlgPack {

///
/** Concrete subclass for a default sub-space of a vector.
 *
 * The default constructor, copy constructor and assignment operator
 * functions are allowed and have the correct behavior.
 *
 * ToDo: Finish Documentation!
 */
class VectorSpaceSubSpace : public virtual VectorSpace {
public:

	///
	/** Calls #this->initialize(...)#.
	 */
	VectorSpaceSubSpace( const space_ptr_t& space, const Range1D& rng );

	///
	/** Initialize.
	 *
	 * Constructs a sub-space of the vector space this = space.sub_space(rng).
	 *
	 * @param  space  [in] The original full vector space.
	 * @param  rng    [in] The range of element that #this# vector sub-space will represent.
	 */
	void initialize( const space_ptr_t& space, const Range1D& rng );

	///
	const space_ptr_t& space() const;

	///
	const Range1D& rng() const;

	/// Validate rng
	void validate_range( const Range1D& rng ) const;

	/** @name Overridden from VectorSpaceBase */
	//@{

	///
	bool is_compatible(const VectorSpaceBase& ) const;

	//@}

	/** @name Overridden form VectorSpace */
	//@{

	///
	index_type dim() const;
	///
	vec_mut_ptr_t create_member() const;
	///
	space_ptr_t sub_space(const Range1D& rng) const;

	//@}

private:

	space_ptr_t     space_;   // If space_.get() == NULL, then uninitalized (dim == 0)
	Range1D         rng_;     // The range of elements from this space to represent!

}; // end class VectorSpaceSubSpace

// //////////////////////////////
// Inline members

inline
const VectorSpace::space_ptr_t& VectorSpaceSubSpace::space() const
{
	return space_;
}

inline
const Range1D& VectorSpaceSubSpace::rng() const
{
	return rng_;
}

#ifndef _DEBUG
inline
void VectorSpaceSubSpace::validate_range(const Range1D&) const
{}
#endif

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_SUB_SPACE_H
