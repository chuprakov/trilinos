// /////////////////////////////////////////////////////////////////
// MatrixSymWithOpGetGMSSymMutable.hpp
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

#ifndef MATRIX_SYM_WITH_OP_GET_GMS_SYM_MUTABLE_H
#define MATRIX_SYM_WITH_OP_GET_GMS_SYM_MUTABLE_H

#include "SparseLinAlgPack/src/MatrixSymWithOpGetGMSSym.hpp"

namespace SparseLinAlgPack {

///
/** Abstract interface that allows the extraction of a non-const <tt>LinAlgPack::sym_gms</tt>
 * view of a symmetry abstract matrix.
 *
 * This interface is ment to be used by <tt>MatrixSymWithOp</tt> objects
 * that store all of their matrix elements in the local address space or can easily
 * access all of the elements from this process and can modify the elements in their
 * data structures.
 *
 * Subclasses that store a BLAS compatible dense symmetric matrix can implement
 * these methods without any dynamic memory allocations.  There is no default
 * implementation for these methods so subclasses that derive from this interface
 * must implement these methods.
 *
 * These methods should never be called directly.  Instead, use the helper
 * class type <tt>MatrixDenseSymMutableEncap</tt>.
 */
class MatrixSymWithOpGetGMSSymMutable : virtual public MatrixSymWithOpGetGMSSym {
public:

	///
	/** Get a non-const view of the symmetric abstract matrix in the form <tt>LinAlgPack::LinAlgPack::sym_gms</tt>.
	 *
	 * @return On ouput, \c return will be initialized to point to storage to the dense matrix elements.
	 * The output from this function <tt>sym_gms_view = this->get_sym_gms_view()</tt> must be passed to
	 * <tt>this->commit_sym_gms_view(gms)</tt> to free any memory that may have been allocated and to ensure
	 * the that underlying abstract matrix object has been updated.
	 * After <tt>this->commit_sym_gms_view(sym_gms_view)</tt> is called, \c sym_gms_view must not be used any longer!
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.rows() == this->rows()</tt>
	 * <li> <tt>return.cols() == this->cols()</tt>
	 * </ul>
	 *
	 * Warning!  If a subclass overrides this method, it must also override \c commit_sym_gms_view().
	 */
	virtual LinAlgPack::sym_gms get_sym_gms_view() = 0;

	///
	/** Free a view of a dense matrix initialized from <tt>get_sym_gms_view()>/tt>.
	 *
	 * @param  sym_gms_view
	 *              [in/out] On input, \c sym_gms_view must have been initialized from \c this->get_sym_gms_view().
	 *              On output, \c sym_gms_view will become invalid and must not be used.
	 *
	 * Preconditions:<ul>
	 * <li> \c sym_gms_view must have been initialized by \c this->get_sym_gms_view)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> \c this is guaranteed to be updated.
	 * <li> \c sym_gms_view becomes invalid and must not be used any longer!
	 * </ul>
	 */
	virtual void commit_sym_gms_view(LinAlgPack::sym_gms* sym_gms_view) = 0;

}; // end class MatrixSymWithOpGetGMSSymMutable

///
/** Helper class type that simplifies the usage of the <tt>MatrixSymWithOpGetGMSSymMutable</tt> interface for clients.
 *
 * This takes care of worrying about if the <tt>MatrixSymWithOpGetGMSSymMutable</tt> interface is supported or not
 * and remembering to free the <tt>LinAlgPack::sym_gms</tt> view properly.
 *
 * This class is only to be used on the stack as an automatic variable.  For example, to extract a
 * <tt>LinAlgPack::sym_gms</tt> view of an abstract vector and use it to set the matrix to a scalar
 * one could write a function like:
 \code
 void assign( const value_type alpha, MatrixSymWithOpGetGMSSymMutable* mat_inout ) {
     MatrixDenseSymMutableEncap  gms_inout(*mat_inout);
	 gms_inout() = alpha;
 }
 \endcode
 * In the above code, if the underlying <tt>MatrixSymWithOpGetGMSSymMutable</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method
 * <tt>MatrixSymWithOpGetGMSSymMutable::get_sym_gms_view()</tt> then the above code will only have a constant
 * time overhead.
 */
class MatrixDenseSymMutableEncap {
public:

	///
	/** Construct a <tt>LinAlgPack::sym_gms</tt> view from a <tt>MatrixSymWithOpGetGMSSymMutable</tt> object.
	 */
	MatrixDenseSymMutableEncap( MatrixSymWithOpGetGMSSymMutable*  mat_get );
	///
	/** Construct a <tt>LinAlgPack::sym_gms</tt> view from a <tt>MatrixSymWithOp</tt> object.
	 *
	 * If <tt>dynamic_cast<MatrixSymWithOpGetGMSSymMutable*>(mat) == NULL</tt> then a ???
	 * exception is thrown.
	 */
	MatrixDenseSymMutableEncap( MatrixSymWithOp* mat );
	/// Frees the <tt>LinAlgPack::sym_gms</tt> view and commits the changes.
	~MatrixDenseSymMutableEncap();
	/// Returns a non-const view of the <tt>LinAlgPack::sym_gms</tt> view.
	LinAlgPack::sym_gms operator()();
	/// Returns a const view of the <tt>LinAlgPack::sym_gms</tt> view.
	const LinAlgPack::sym_gms operator()() const;

private:

	MatrixSymWithOpGetGMSSymMutable    *mat_get_;
	LinAlgPack::sym_gms                sym_gms_view_;
	MatrixDenseSymMutableEncap();                                             // Not defined and not to be called!
	MatrixDenseSymMutableEncap(const MatrixDenseSymMutableEncap&);            // ""
	MatrixDenseSymMutableEncap& operator=(const MatrixDenseSymMutableEncap&); // ""

}; // end class MatrixDenseSymMutableEncap

// ///////////////////////////////////////////
// Inline members

// MatrixDenseSymMutableEncap

inline
MatrixDenseSymMutableEncap::MatrixDenseSymMutableEncap( MatrixSymWithOpGetGMSSymMutable*  mat_get )
	:mat_get_(mat_get)
	,sym_gms_view_(mat_get_->get_sym_gms_view())
{}

inline
MatrixDenseSymMutableEncap::MatrixDenseSymMutableEncap( MatrixSymWithOp* mat )
	:mat_get_(&DynamicCastHelperPack::dyn_cast<MatrixSymWithOpGetGMSSymMutable>(*mat))
	,sym_gms_view_(mat_get_->get_sym_gms_view())
{}

inline
MatrixDenseSymMutableEncap::~MatrixDenseSymMutableEncap()
{
	mat_get_->commit_sym_gms_view(&sym_gms_view_);
}

inline
LinAlgPack::sym_gms MatrixDenseSymMutableEncap::operator()()
{
	return sym_gms_view_;
}

inline
const LinAlgPack::sym_gms MatrixDenseSymMutableEncap::operator()() const
{
	return sym_gms_view_;
}

} // end namespace SparseLinAlgPack

#endif // MATRIX_SYM_WITH_OP_GET_GMS_SYM_MUTABLE_H
