// /////////////////////////////////////////////////////////////////
// MatrixWithOpGetGMSTri.hpp
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

#ifndef MATRIX_WITH_OP_GET_GMS_TRI_H
#define MATRIX_WITH_OP_GET_GMS_TRI_H

#include "SparseLinAlgPackTypes.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOp.hpp"
#include "LinAlgPack/src/GenMatrixAsTriSym.hpp"
#include "dynamic_cast_verbose.hpp"

namespace SparseLinAlgPack {

///
/** Mix-in interface that allows the extraction of a const <tt>LinAlgPack::tri_gms</tt>
 * view of an non-singular abstract matrix.
 *
 * This interface is ment to be used by <tt>MatrixWithOp</tt> objects
 * that store all of their matrix elements in the local address space or can easily
 * access all of the elements from this process.
 *
 * Subclasses that store a BLAS compatible triangular dense matrix can implement
 * these methods without any dynamic memory allocations.  There is no default
 * implementation for these methods so subclasses that derive from this interface
 * must implement these methods.
 *
 * These methods should never be called directly.  Instead, use the helper
 * class type <tt>MatrixDenseTriEncap</tt>.
 */
class MatrixWithOpGetGMSTri
	: virtual public AbstractLinAlgPack::MatrixWithOp // doxygen needs full name
{
public:

	///
	/** Get a const view of the symmetric abstract matrix in the form <tt>LinAlgPack::tri_gms</tt>.
	 *
	 * @return On ouput, \c return will be initialized to point to storage to the symmetric dense
	 *  matrix elements.
	 * The output from this function <tt>tri_gms_view = this->get_tri_gms_view()</tt> must be passed to
	 * <tt>this->free_tri_gms_view(gms)</tt> to free any memory that may have been allocated.
	 * After <tt>this->free_gms_view(gms_view)</tt> is called, \c gms_view must not be used any longer!
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.rows() == this->rows()</tt>
	 * <li> <tt>return.cols() == this->cols()</tt>
	 * </ul>
	 *
	 * Warning!  If a subclass overrides this method, it must also override \c free_tri_gms_view().
	 */
	virtual const LinAlgPack::tri_gms get_tri_gms_view() const = 0;

	///
	/** Free a view of a symmetric dense matrix initialized from <tt>get_tri_gms_view()>/tt>.
	 *
	 * @param  tri_gms_view
	 *              [in/out] On input, \c tri_gms_view must have been initialized from \c this->get_tri_gms_view().
	 *              On output, \c tri_gms_view will become invalid and must not be used.
	 *
	 * Preconditions:<ul>
	 * <li> \c tri_gms_view must have been initialized by \c this->get_tri_gms_view)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> \c tri_gms_view becomes invalid and must not be used any longer!
	 * </ul>
	 */
	virtual void free_tri_gms_view(const LinAlgPack::tri_gms* tri_gms_view) const = 0;

}; // end class MatrixWithOpGetGMSTri

///
/** Helper class type that simplifies the usage of the <tt>MatrixWithOpGetGMSTri</tt> interface for clients.
 *
 * This takes care of worrying about if the <tt>MatrixWithOpGetGMSTri</tt> interface is supported or not
 * and remembering to free the <tt>LinAlgPack::tri_gms</tt> view properly.
 *
 * This class is only to be used on the stack as an automatic variable.  For example, to extract a
 * <tt>LinAlgPack::tri_gms</tt> view of an abstract vector and use it to call another function
 * one could write a function like:
 \code
 void call_func(const MatrixWithOpGetGMSTri& mat_in ) {
     func( MatrixDenseTriEncap(mat_in)() );
 }
 \endcode
 * In the above code, if the underlying <tt>MatrixWithOpGetGMSTri</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method
 * <tt>MatrixWithOpGetGMSTri::get_tri_gms_view()</tt> then the above code will only have a constant
 * time overhead.
 */
class MatrixDenseTriEncap {
public:

	///
	/** Construct a <tt>LinAlgPack::tri_gms</tt> view from a <tt>MatrixWithOpGetGMSTri</tt> object.
	 */
	MatrixDenseTriEncap( const MatrixWithOpGetGMSTri&  mat_get );
	///
	/** Construct a <tt>LinAlgPack::tri_gms</tt> view from a <tt>MatrixWithOp</tt> object.
	 *
	 * If <tt>dynamic_cast<const MatrixWithOpGetGMSTri*>(&mat) == NULL</tt> then a ???
	 * exception is thrown.
	 */
	MatrixDenseTriEncap( const MatrixWithOp& mat );
	/// Frees the <tt>LinAlgPack::tri_gms</tt> view.
	~MatrixDenseTriEncap();
	/// Returns a constant view of the <tt>LinAlgPack::tri_gms</tt> view.
	const LinAlgPack::tri_gms operator()() const;

private:

	const MatrixWithOpGetGMSTri     &mat_get_;
	const LinAlgPack::tri_gms       tri_gms_view_;
	MatrixDenseTriEncap();                                      // Not defined and not to be called!
	MatrixDenseTriEncap(const MatrixDenseTriEncap&);               // ""
	MatrixDenseTriEncap& operator=(const MatrixDenseTriEncap&);    // ""

}; // end class MatrixDenseTriEncap

// ///////////////////////////////////////////
// Inline members

// MatrixDenseTriEncap

inline
MatrixDenseTriEncap::MatrixDenseTriEncap( const MatrixWithOpGetGMSTri&  mat_get )
	:mat_get_(mat_get)
	,tri_gms_view_(mat_get_.get_tri_gms_view())
{}

inline
MatrixDenseTriEncap::MatrixDenseTriEncap( const MatrixWithOp& mat )
	:mat_get_(DynamicCastHelperPack::dyn_cast<const MatrixWithOpGetGMSTri>(mat))
	,tri_gms_view_(mat_get_.get_tri_gms_view())
{}

inline
MatrixDenseTriEncap::~MatrixDenseTriEncap()
{
	mat_get_.free_tri_gms_view(&tri_gms_view_);
}

inline
const LinAlgPack::tri_gms MatrixDenseTriEncap::operator()() const
{
	return tri_gms_view_;
}

} // end namespace SparseLinAlgPack

#endif // MATRIX_WITH_OP_GET_GMS_TRI_H
