// /////////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorSpaceSubSpace.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef VECTOR_SPACE_SUB_SPACE_H
#define VECTOR_SPACE_SUB_SPACE_H

#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Thyra_Range1D.hpp"

namespace AbstractLinAlgPack {

///
/** Concrete subclass for a default sub-space of a vector.
 *
 * There is not much to this subclass.  It basically implements all
 * of its methods based on the external VectorSpace interface to
 * implement is_compatible() and sub_space() and and relys
 * on a default subclass VectorMutableSubView to implement
 * create_member(). 
 *
 * The default constructor, copy constructor and assignment operator
 * functions are allowed and have the correct behavior.
 */
class VectorSpaceSubSpace : public virtual VectorSpace {
public:

	///
	/** Constructs to uninitialized.
	 *
	 * Postconditions: see \c set_uninitialized().
	 */
	VectorSpaceSubSpace();

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorSpaceSubSpace( const space_ptr_t& full_space, const Range1D& rng );

	///
	/** Initialize.
	 *
	 * Constructs a sub-space of the vector space this = space.sub_space(rng).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>full_space.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>).
	 * <li> [rng.full_range() == false</tt>] <tt>rng.lbound() <= full_space->dim()</tt> (throw <tt>std::out_of_range</tt>).
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [&& rng.full_range() == true</tt>] <tt>this->dim() == full_space->dim()</tt>
	 * <li> [&& rng.full_range() == false</tt>] <tt>this->dim() == rng.size()</tt>
	 * </ul>
	 *
	 * @param  full_space  [in] The original full vector space (must be <tt>full_space.get() != NULL</tt>).
	 * @param  rng         [in] The range of element that <tt>this</tt> vector sub-space will represent.
	 */
	void initialize( const space_ptr_t& full_space, const Range1D& rng );

	///
	/** Set uninitilized.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->dim() == 0</tt>
	 * <li> <tt>this->create_member().get() == NULL</tt>
	 * </ul>
	 */
	void set_uninitialized();

	///
	const space_ptr_t& full_space() const;

	///
	const Range1D& rng() const;

	/// Validate rng
	void validate_range( const Range1D& rng ) const;

	/** @name Overridden from VectorSpace */
	//@{

	///
	bool is_compatible(const VectorSpace& ) const;
	///
	bool is_in_core() const;
	///
	index_type dim() const;
	///
	vec_mut_ptr_t create_member() const;
	///
	space_ptr_t clone() const;
	///
	space_ptr_t sub_space(const Range1D& rng) const;

	//@}

private:

	space_ptr_t     full_space_;   ///< If space_.get() == NULL, then uninitalized (dim == 0)
	Range1D         rng_;          ///< The range of elements from this space to represent!

}; // end class VectorSpaceSubSpace

// //////////////////////////////
// Inline members

inline
VectorSpaceSubSpace::VectorSpaceSubSpace()
	: rng_(Range1D::Invalid)
{}

inline
const VectorSpace::space_ptr_t& VectorSpaceSubSpace::full_space() const
{
	return full_space_;
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
