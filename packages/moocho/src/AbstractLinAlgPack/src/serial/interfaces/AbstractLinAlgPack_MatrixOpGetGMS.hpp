// /////////////////////////////////////////////////////////////////
// MatrixWithOpGetGMS.hpp
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

#ifndef MATRIX_WITH_OP_GET_GMS_H
#define MATRIX_WITH_OP_GET_GMS_H

#include "SparseLinAlgPackTypes.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOp.hpp"
#include "DenseLinAlgPack/src/DMatrixClass.hpp"
#include "dynamic_cast_verbose.hpp"

namespace SparseLinAlgPack {

///
/** Abstract interface that allows the extraction of a const <tt>DMatrixSlice</tt>
 * view of an abstract matrix.
 *
 * This interface is ment to be used by <tt>MatrixWithOp</tt> objects
 * that store all of their matrix elements in the local address space or can easily
 * access all of the elements from this process.
 *
 * Subclasses that store a Fortran compatible dense dense matrix can implement
 * these methods without any dynamic memory allocations.  There is no default
 * implementation for these methods so subclasses that derive from this interface
 * must implement these methods.
 *
 * These methods should never be called directly.  Instead, use the helper
 * class type <tt>MatrixDenseEncap</tt>.
 */
class MatrixWithOpGetGMS 
	: virtual public AbstractLinAlgPack::MatrixWithOp // doxygen needs full path
{
public:

	///
	/** Get a const view of the abstract matrix in the form <tt>DenseLinAlgPack::DMatrixSlice</tt>.
	 *
	 * @return On ouput, \c return will be initialized to point to storage to the dense matrix elements.
	 * The output from this function <tt>gms_view = this->get_gms_view()</tt> must be passed to
	 * <tt>this->free_gms_view(gms)</tt> to free any memory that may have been allocated.
	 * After <tt>this->free_gms_view(gms_view)</tt> is called, \c gms_view must not be used any longer!
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.rows() == this->rows()</tt>
	 * <li> <tt>return.cols() == this->cols()</tt>
	 * </ul>
	 *
	 * Warning!  If a subclass overrides this method, it must also override \c free_gms_view().
	 */
	virtual const DMatrixSlice get_gms_view() const = 0;

	///
	/** Free a view of a dense matrix initialized from <tt>get_gms_view()>/tt>.
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
	 * <li> \c gms_view becomes invalid and must not be used any longer!
	 * </ul>
	 */
	virtual void free_gms_view(const DMatrixSlice* gms_view) const = 0;

}; // end class MatrixWithOpGetGMS

///
/** Helper class type that simplifies the usage of the <tt>MatrixWithOpGetGMS</tt> interface for clients.
 *
 * This takes care of worrying about if the <tt>MatrixWithOpGetGMS</tt> interface is supported or not
 * and remembering to free the <tt>DMatrixSlice</tt> view properly.
 *
 * This class is only to be used on the stack as an automatic variable.  For example, to extract a
 * <tt>DMatrixSlice</tt> view of an abstract vector and use it to copy to a <tt>DMatrix</tt>
 * object you could write a function like:
 \code
 void copy(const MatrixWithOpGetGMS& mat_in, GenMatrixClass* gms_out ) {
     MatrixDenseEncap  gms_in(mat_in);
	 *gms_out = gms_in();
 }
 \endcode
 * In the above code, if the underlying <tt>MatrixWithOpGetGMS</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method
 * <tt>MatrixWithOpGetGMS::get_gms_view()</tt> then the above code will only have a constant
 * time overhead.
 */
class MatrixDenseEncap {
public:

	///
	/** Construct a <tt>DMatrixSlice</tt> view from a <tt>MatrixWithOpGetGMS</tt> object.
	 */
	MatrixDenseEncap( const MatrixWithOpGetGMS&  mat_get );
	///
	/** Construct a <tt>DMatrixSlice</tt> view from a <tt>MatrixWithOp</tt> object.
	 *
	 * If <tt>dynamic_cast<const MatrixWithOpGetGMS*>(&mat) == NULL</tt> then a ???
	 * exception is thrown.
	 */
	MatrixDenseEncap( const MatrixWithOp& mat );
	/// Frees the <tt>DMatrixSlice</tt> view.
	~MatrixDenseEncap();
	/// Returns a constant view of the <tt>DMatrixSlice</tt> view.
	const DMatrixSlice operator()() const;

private:

	const MatrixWithOpGetGMS     &mat_get_;
	const DMatrixSlice         gms_view_;
	MatrixDenseEncap();                                      // Not defined and not to be called!
	MatrixDenseEncap(const MatrixDenseEncap&);               // ""
	MatrixDenseEncap& operator=(const MatrixDenseEncap&);    // ""

}; // end class MatrixDenseEncap

// ///////////////////////////////////////////
// Inline members

// MatrixDenseEncap

inline
MatrixDenseEncap::MatrixDenseEncap( const MatrixWithOpGetGMS&  mat_get )
	:mat_get_(mat_get)
	,gms_view_(mat_get_.get_gms_view())
{}

inline
MatrixDenseEncap::MatrixDenseEncap( const MatrixWithOp& mat )
	:mat_get_(DynamicCastHelperPack::dyn_cast<const MatrixWithOpGetGMS>(mat))
	,gms_view_(mat_get_.get_gms_view())
{}

inline
MatrixDenseEncap::~MatrixDenseEncap()
{
	mat_get_.free_gms_view(&gms_view_);
}

inline
const DMatrixSlice MatrixDenseEncap::operator()() const
{
	return gms_view_;
}

} // end namespace SparseLinAlgPack

#endif // MATRIX_WITH_OP_GET_GMS_H
