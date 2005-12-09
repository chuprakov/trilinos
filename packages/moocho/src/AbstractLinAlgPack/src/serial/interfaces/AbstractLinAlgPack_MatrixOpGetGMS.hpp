// /////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixOpGetGMS.hpp
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

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract interface that allows the extraction of a const <tt>DMatrixSlice</tt>
 * view of an abstract matrix.
 *
 * This interface is ment to be used by <tt>MatrixOp</tt> objects
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
class MatrixOpGetGMS 
	: virtual public AbstractLinAlgPack::MatrixOp // doxygen needs full path
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

}; // end class MatrixOpGetGMS

///
/** Helper class type that simplifies the usage of the <tt>MatrixOpGetGMS</tt> interface for clients.
 *
 * This takes care of worrying about if the <tt>MatrixOpGetGMS</tt> interface is supported or not
 * and remembering to free the <tt>DMatrixSlice</tt> view properly.
 *
 * This class is only to be used on the stack as an automatic variable.  For example, to extract a
 * <tt>DMatrixSlice</tt> view of an abstract vector and use it to copy to a <tt>DMatrix</tt>
 * object you could write a function like:
 \code
 void copy(const MatrixOpGetGMS& mat_in, GenMatrixClass* gms_out ) {
     MatrixDenseEncap  gms_in(mat_in);
	 *gms_out = gms_in();
 }
 \endcode
 * In the above code, if the underlying <tt>MatrixOpGetGMS</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method
 * <tt>MatrixOpGetGMS::get_gms_view()</tt> then the above code will only have a constant
 * time overhead.
 */
class MatrixDenseEncap {
public:

	///
	/** Construct a <tt>DMatrixSlice</tt> view from a <tt>MatrixOpGetGMS</tt> object.
	 */
	MatrixDenseEncap( const MatrixOpGetGMS&  mat_get );
	///
	/** Construct a <tt>DMatrixSlice</tt> view from a <tt>MatrixOp</tt> object.
	 *
	 * If <tt>dynamic_cast<const MatrixOpGetGMS*>(&mat) == NULL</tt> then a ???
	 * exception is thrown.
	 */
	MatrixDenseEncap( const MatrixOp& mat );
	/// Frees the <tt>DMatrixSlice</tt> view.
	~MatrixDenseEncap();
	/// Returns a constant view of the <tt>DMatrixSlice</tt> view.
	const DMatrixSlice operator()() const;

private:

	const MatrixOpGetGMS     &mat_get_;
	const DMatrixSlice         gms_view_;
	MatrixDenseEncap();                                      // Not defined and not to be called!
	MatrixDenseEncap(const MatrixDenseEncap&);               // ""
	MatrixDenseEncap& operator=(const MatrixDenseEncap&);    // ""

}; // end class MatrixDenseEncap

// ///////////////////////////////////////////
// Inline members

// MatrixDenseEncap

inline
MatrixDenseEncap::MatrixDenseEncap( const MatrixOpGetGMS&  mat_get )
	:mat_get_(mat_get)
	,gms_view_(mat_get_.get_gms_view())
{}

inline
MatrixDenseEncap::MatrixDenseEncap( const MatrixOp& mat )
	:mat_get_(Teuchos::dyn_cast<const MatrixOpGetGMS>(mat))
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

} // end namespace AbstractLinAlgPack

#endif // MATRIX_WITH_OP_GET_GMS_H
