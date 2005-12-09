// ////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorSpaceSerial.hpp
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

#ifndef VECTOR_SPACE_SERIAL_H
#define VECTOR_SPACE_SERIAL_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

///
/** Subclass for serial vector space objects that create <tt>VectorMutableDense</tt>
 * vector and <tt>MultiVectorMutableDense</tt> multi-vector objects.
 *
 * The default constructor, copy constructor and assignment operators
 * are allowed since they have the correct semantics.
 */
class VectorSpaceSerial
	: public AbstractLinAlgPack::VectorSpace
{
public:

	/** @name Constructors / initializers */
	//@{

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	VectorSpaceSerial( size_type dim = 0 );

	///
	/** Initialize given the dimension of the vector space.
	 *
	 * @param  dim   [in] The dimension of the vector space.
	 */
	void initialize( size_type dim );

	//@}

	/** @name Overridden from VectorSpece */
	//@{

	///
	/** Returns true if <tt>vec_space.dim() == this->dim()</tt>.
	 *
	 * The assumption here is that <tt>Vector::get_sub_vector()</tt>
	 * and <tt>VectorMutable::get_sub_vector()</tt> can be used to implement
	 * all of the methods on an SMP machine in an efficient manner.
	 */
 	bool is_compatible(const VectorSpace& vec_space) const;
	/// Returns true
	bool is_in_core() const;
	/// Returns 0 if uninitialized
	index_type dim() const;
	/// Returns a <tt>VectorSpaceFactorySerial</tt> object
	space_fcty_ptr_t small_vec_spc_fcty() const;
	///
	space_ptr_t clone() const;
	/// Returns a <tt>VectorMutableDense</tt> object.
	vec_mut_ptr_t create_member() const;
	/// Returns a <tt>MultiVectorMutableDense</tt> object.
	multi_vec_mut_ptr_t create_members(size_type num_vecs) const;
	///
	space_ptr_t sub_space(const Range1D& rng) const;
	///
	space_ptr_t space(
		const GenPermMatrixSlice  &P
		,BLAS_Cpp::Transp         P_trans
		) const;
	//@}

private:

	size_type     dim_;

}; // end class VectorSpaceSerial

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_SERIAL_H
