// /////////////////////////////////////////////////////////////////
// MatrixWithOpGetGMSMutable.h
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

#ifndef MATRIX_WITH_OP_GET_GMS_MUTABLE_H
#define MATRIX_WITH_OP_GET_GMS_MUTABLE_H

#include "SparseLinAlgPack/include/MatrixWithOpGetGMS.h"

namespace SparseLinAlgPack {

///
/** Abstract interface that allows the extraction of a non-const <tt>GenMatrixSlice</tt>
 * view of an abstract matrix.
 *
 * This interface is ment to be used by <tt>MatrixWithOp</tt> objects
 * that store all of their matrix elements in the local address space or can easily
 * access all of the elements from this process and can modify the elements in their
 * data structures.
 *
 * Subclasses that store a Fortran compatible dense dense matrix can implement
 * these methods without any dynamic memory allocations.  There is no default
 * implementation for these methods so subclasses that derive from this interface
 * must implement these methods.
 *
 * These methods should never be called directly.  Instead, use the helper
 * class type <tt>MatrixDenseMutableEncap</tt>.
 */
class MatrixWithOpGetGMSMutable : virtual public MatrixWithOpGetGMS {
public:

	///
	using MatrixWithOpGetGMS::get_gms_view;

	///
	/** Get a representation of the abstract matrixr in the form <tt>LinAlgPack::GenMatrixSlice</tt>.
	 *
	 * @return On ouput, \c return will be initialized to point to storage to the dense matrix elements.
	 * The output from this function <tt>gms_view = this->get_gms_view()</tt> must be passed to
	 * <tt>this->commit_gms_view(gms)</tt> to commit and free any memory that may have been allocated
	 * and to ensure the that underlying abstract matrix object has been updated.
	 * After <tt>this->commit_gms_view(gms_view)</tt> is called, \c gms_view must not be used any longer!
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.rows() == this->rows()</tt>
	 * <li> <tt>return.cols() == this->cols()</tt>
	 * </ul>
	 *
	 * Warning!  If a subclass overrides this method, it must also override \c commit_gms_view().
	 */
	virtual GenMatrixSlice get_gms_view() = 0;

	///
	/** Commit changes to a view of a dense matrix initialized from <tt>this->get_gms_view()</tt>.
	 *
	 * @param  gms_view
	 *              [in/out] On input, \c gms_view must have been initialized from \c this->get_gms_view().
	 *              On output, \c gms_view will become invalid and must not be used.
	 *
	 * Preconditions:<ul>
	 * <li> \c gms_view must have been initialized by \c this->get_gms_view)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> \c this is guaranteed to be updated.
	 * <li> \c gms_view becomes invalid and must not be used any longer!
	 * </ul>
	 */
	virtual void commit_gms_view(GenMatrixSlice* gms_view) = 0;

}; // end class MatrixWithOpGetGMSMutable

///
/** Helper class type that simplifies the usage of the <tt>MatrixWithOpGetGMSMutable</tt> interface for clients.
 *
 * This takes care of worrying about if the <tt>MatrixWithOpGetGMSMutable</tt> interface is supported or not
 * and remembering to free the <tt>GenMatrixSlice</tt> view properly.
 *
 * This class is only to be used on the stack as an automatic variable.  For example, to extract a
 * <tt>GenMatrixSlice</tt> view of an abstract vector and use it to set the matrix to a scalar
 * one could write a function like:
 \code
 void assign( const value_type alpha, MatrixWithOpGetGMSMutable* mat_inout ) {
     MatrixDenseMutableEncap  gms_inout(*mat_inout);
	 gms_inout() = alpha;
 }
 \endcode
 * In the above code, if the underlying <tt>MatrixWithOpGetGMSMutable</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method
 * <tt>MatrixWithOpGetGMSMutable::get_gms_view()</tt> then the above code will only have a constant
 * time overhead.
 */
class MatrixDenseMutableEncap {
public:

	///
	/** Construct a <tt>GenMatrixSlice</tt> view from a <tt>MatrixWithOpGetGMSMutable</tt> object.
	 */
	MatrixDenseMutableEncap( MatrixWithOpGetGMSMutable*  mat_get );
	///
	/** Construct a <tt>GenMatrixSlice</tt> view from a <tt>MatrixWithOp</tt> object.
	 *
	 * If <tt>dynamic_cast<MatrixWithOpGetGMSMutable*>(mat) == NULL</tt> then a ???
	 * exception is thrown.
	 */
	MatrixDenseMutableEncap( MatrixWithOp* mat );
	/// Frees the <tt>GenMatrixSlice</tt> view and commits the changes.
	~MatrixDenseMutableEncap();
	/// Returns a non-const view of the <tt>GenMatrixSlice</tt> view.
	GenMatrixSlice operator()();
	/// Returns a const view of the <tt>GenMatrixSlice</tt> view.
	const GenMatrixSlice operator()() const;

private:

	MatrixWithOpGetGMSMutable     *mat_get_;
	GenMatrixSlice                gms_view_;
	MatrixDenseMutableEncap();                                          // Not defined and not to be called!
	MatrixDenseMutableEncap(const MatrixDenseMutableEncap&);            // ""
	MatrixDenseMutableEncap& operator=(const MatrixDenseMutableEncap&); // ""

}; // end class MatrixDenseMutableEncap

// ///////////////////////////////////////////
// Inline members

// MatrixDenseMutableEncap

inline
MatrixDenseMutableEncap::MatrixDenseMutableEncap( MatrixWithOpGetGMSMutable*  mat_get )
	:mat_get_(mat_get)
	,gms_view_(mat_get_->get_gms_view())
{}

inline
MatrixDenseMutableEncap::MatrixDenseMutableEncap( MatrixWithOp* mat )
	:mat_get_(&DynamicCastHelperPack::dyn_cast<MatrixWithOpGetGMSMutable>(*mat))
	,gms_view_(mat_get_->get_gms_view())
{}

inline
MatrixDenseMutableEncap::~MatrixDenseMutableEncap()
{
	mat_get_->commit_gms_view(&gms_view_);
}

inline
GenMatrixSlice MatrixDenseMutableEncap::operator()()
{
	return gms_view_;
}

inline
const GenMatrixSlice MatrixDenseMutableEncap::operator()() const
{
	return gms_view_;
}

} // end namespace SparseLinAlgPack

#endif // MATRIX_WITH_OP_GET_GMS_MUTABLE_H
