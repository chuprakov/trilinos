// ////////////////////////////////////////////////////////////////////////////
// VectorSpaceSerial.h
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

#include "SparseLinAlgPack/include/SparseLinAlgPackTypes.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"

namespace SparseLinAlgPack {

///
/** Subclass for serial vector space objects that create <tt>VectorWithOpMutableDense</tt>
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
	 * The assumption here is that <tt>VectorWithOp::get_sub_vector()</tt>
	 * and <tt>VectorWithOpMutable::get_sub_vector()</tt> can be used to implement
	 * all of the methods on an SMP machine in an efficient manner.
	 */
 	bool is_compatible(const VectorSpace& vec_space) const;
	/// Returns 0 if uninitialized
	index_type dim() const;
	///
	space_ptr_t clone() const;
	/// Returns a <tt>VectorWithOpMutableDense</tt> object.
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

} // end namespace SparseLinAlgPack

#endif // VECTOR_SPACE_SERIAL_H
