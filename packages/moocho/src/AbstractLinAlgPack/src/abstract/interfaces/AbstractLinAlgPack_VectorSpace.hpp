// //////////////////////////////////////////
// VectorSpace.h

#ifndef VECTOR_SPACE_H
#define VECTOR_SPACE_H

#include "VectorSpaceBase.h"
#include "Range1D.h"

namespace AbstractLinAlgPack {

///
/** Abstract interface for objects that represent a space for mutable coordinate vectors.
 *
 * This class is primary an "Abstract Factory" interface.  A <tt>VectorSpace</tt> object can exist
 * independent from any individual <tt>VectorWithOpMutable</tt> object.  Or, a <tt>VectorSpace</tt> object
 * can have a lifetime that is dependent on a single <tt>VectorWithOpMutable</tt> object.  The
 * same interface can serve both roles.  
 */
class VectorSpace : public virtual VectorSpaceBase {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<const VectorSpace>    space_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<VectorWithOpMutable>  vec_mut_ptr_t;

	///
	/** Return the dimmension of the vector space.
	 */
	virtual index_type dim() const = 0;

	///
	/** Create a member of the vector space.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return->dim() == this->dim()</tt>
	 * <li> <tt>return->get_ele(i) == 0.0</tt>, for <tt>i = 1...this->dim()</tt>
	 * <li> <tt>return->space().is_compatible(*this) == true</tt>
	 * </ul>
	 *
	 * @return  Returns a smart reference counted pointer to a dynamically
	 * allocated vector object.  After construction the values returnd by 
	 * <tt>return->get_ele(i)</tt> are unspecified (uninitialized).  This allows for
	 * faster execution times.  Note that <tt>&return->space()</tt> does not have to
	 * be equal to <tt>this</tt>.
	 */
	virtual vec_mut_ptr_t create_member() const = 0;

	///
	/** Create a transient sub-space of the current vector space.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>rng.ubound() <= this->dim()</tt> (<tt>throw std::out_of_range</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return->dim() == rng->size()</tt>
	 * </ul>
	 *
	 * @param  rng  [in] The range of the elements to extract a vector sub-space.
	 *
	 * @return  Returns a smart reference counted pointer to a dynamically
	 * allocated vector space object.  Note that the vector object returned
	 * by <tt>this->sub_space(rng).create_member()</tt> should be exactly equivalent
	 * to the vector returned by
	 * <tt>this->create_member()->sub_view(rng)->space()->create_member()</tt>.
	 * It is allowed for the implementation to return <tt>return->get() == NULL</tt>
	 * for arbitrary values of <tt>rng</tt>.  Only some <tt>rng</tt> ranges may be allowed
	 * but they will be appropriate for the application at hand.  However, a
	 * very good implementation should be able to accomidate any valid <tt>rng</tt>
	 * that meets the basic preconditions.  The default implementation uses
	 * the subclass VectorSpaceSubSpace to represent any arbitrary
	 * sub-space but this can be very inefficient if the sub-space is
	 * very small compared this this full vector space.
	 *
	 * Note that if two vector space objects <tt>X</tt> and <tt>Y</tt> are compatible (i.e.
	 * <tt>X.is_compatible(Y) == true</tt>, then it is also expected that
	 * <tt>X.sub_space(rng)->is_compatible(*Y.sub_space(rng))</tt> will also be true.
	 * However, in general it can not be expected that
	 * <tt>X.sub_space(rng1)->is_compatible(*X.sub_space(rng2))</tt>, with
	 * <tt>rng1.size() == rng2.size()</tt> true.  For serial vectors, it may
	 * be but for parallel vectors it will most certainly not be.  Therefore, in
	 * general, don't assume that arbitrary subsets of the vector spaces will be
	 * compatible, even if the sizes of these subspaces are the same.
	 */
	virtual space_ptr_t sub_space(const Range1D& rng) const;

	/// Inlined to call <tt>this->sub_space(Range1D(il,iu))</tt>.
	space_ptr_t sub_space( const index_type il, const index_type iu ) const;
	
	/** @name Overrriden from VectorSpaceBase */
	//@{
	
	///
	/** This implementaion just calls <tt>create_member()</tt> and then converts
	 * the smart pointer.
	 */
	vec_ptr_t new_member() const;

	//@}

}; // end class VectorSpace

// ////////////////////////////////////////////////
// Inline members

inline
VectorSpace::space_ptr_t
VectorSpace::sub_space( const index_type il, const index_type iu ) const
{
	return this->sub_space(Range1D(il,iu));
}

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_SPACE_H
