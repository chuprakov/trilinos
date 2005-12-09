// ////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MatrixIdentConcatStd.hpp
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

#ifndef MATRIX_IDENT_CONCAT_STD_H
#define MATRIX_IDENT_CONCAT_STD_H

#include "ConstrainedOptPack_MatrixIdentConcat.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Thyra_Range1D.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace ConstrainedOptPack {

///
/** Concrete implementation class for a matrix vertically concatonated with an identity matrix.
 *
 * Represents an interface for a matrix that represents:
 \verbatim
 
 M = [ alpha*op(D) ]    (TOP)
     [     I       ]

 or

 M = [     I       ]
     [ alpha*op(D) ]   (BOTTOM)
 \endverbatim
 *
 * This subclass allows a client to set the representation matrix \c D.
 */
class MatrixIdentConcatStd : public MatrixIdentConcat {
public:

	/** @name Public types */
	//@{
	///
	enum ETopBottom { TOP, BOTTOM };
	///
	typedef Teuchos::RefCountPtr<const MatrixOp> D_ptr_t;
	//@}

	/** @name Constructors/initializers. */
	//@{
	///
	/** Constructs to uninitialized.
	 */
	MatrixIdentConcatStd();
	///
	/** Setup with a matrix object.
	 *
	 * @param  top_or_bottom
	 *                 [in] If <tt>TOP</tt> then <tt>M = [ alpha*op(D); I ]</tt> and if <tt>BOTTOM</tt> then
	 *                  <tt>M = [ I; alpha*op(D) ]</tt>
	 * @param  alpha   [in] Scalar that multiplies \c op(D)
	 * @param  D_ptr   [in] Smart pointer to a matrix object that represents \c D.  The matrix object pointed to must
	 *                 not be altered until \c this object is no longer in use or <tt>this->set_uninitialized()</tt>
	 *                 has been called.
	 * @param  D_trans [in] Determines if <tt>op(D) = D</tt> (\c no_trans#) or <tt>op(D) = D'</tt> (\c trans).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>D.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>space_cols->dim() == rows(op(D)) + cols(op(D))</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>space_rows->dim() == cols(op(D))</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>space_cols->sub_space(D_rng)->is_compatible(op(D).space_cols())</tt>
	 *      (throw <tt>std::invalid_argument</tt>) See \c D_rng defined below
	 * <li> <tt>space_rows->is_compatible(op(D).space_rows())</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->D_ptr().get() == D_ptr.get()</tt>
	 * <li> <tt>&this->D()          == this->D_ptr().get()</tt>
	 * <li> <tt>this->D_trans()     == D_trans</tt>
	 * <li> <tt>this->alpha()       == alpha</tt>
	 * <li> <tt>this->rows()        == rows(op(D)) + cols(op(D))</tt>
	 * <li> <tt>this->cols()        == cols(op(D))</tt>
	 * <li> <tt>&this->space_cols() == space_cols.get()</tt>
	 * <li> <tt>&this->space_rows() == space_rows.get()</tt>
	 * <li> [<tt>top_or_bottom == TOP</tt>]    <tt>this->D_rng() = [1,rows(op(D))]</tt>
	 * <li> [<tt>top_or_bottom == TOP</tt>]    <tt>this->I_rng() = [rows(op(D))+1,rows(op(D))+cols(op(D))]</tt>
	 * <li> [<tt>top_or_bottom == BOTTOM</tt>] <tt>this->D_rng() = [cols(op(D))+1,rows(op(D))+cols(op(D))]</tt>
	 * <li> [<tt>top_or_bottom == BOTTOM</tt>] <tt>this->I_rng() = [1,cols(op(D))]</tt>
	 * </ul>
	 */
	virtual void initialize(
		const VectorSpace::space_ptr_t&    space_cols
		,const VectorSpace::space_ptr_t&   space_rows
		,ETopBottom                        top_or_bottom
		,value_type                        alpha
		,const D_ptr_t                     &D_ptr
		,BLAS_Cpp::Transp                  D_trans
		);
	///
	/** Set the matrix to uninitialized.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->space_cols()</tt> throws an exception
	 * <li> <tt>this->space_rows()</tt> throws an exception
	 * <li> <tt>this->D_ptr().get() == NULL</tt>
	 * <li> <tt>&this->D()</tt> throws an exception
	 * <li> <tt>this->D_trans() == no_trans</tt>
	 * <li> <tt>this->alpha() == 0.0</tt>
	 * <li> <tt>this->rows() == 0</tt>
	 * <li> <tt>this->cols() == 0</tt>
	 * <li> [<tt>top_or_bottom == TOP</tt>]    <tt>this->D_rng() = [1,rows(op(D))]</tt>
	 * <li> [<tt>top_or_bottom == TOP</tt>]    <tt>this->I_rng() = [rows(op(D))+1,rows(op(D))+cols(op(D))]</tt>
	 * <li> [<tt>top_or_bottom == BOTTOM</tt>] <tt>this->D_rng() = [cols(op(D))+1,rows(op(D))+cols(op(D))]</tt>
	 * <li> [<tt>top_or_bottom == BOTTOM</tt>] <tt>this->I_rng() = [1,cols(op(D))]</tt>
	 * </ul>
	 */
	virtual void set_uninitialized();
	///
	/** Return the smart reference counted point to the \c D matrix.
	 *
	 * If the matrix object \c D is owned exclusively by \c this matrix object
	 * then <tt>this->D_ptr().count() == 1</tt> on return.
	 */
	virtual const D_ptr_t& D_ptr() const;
	//@}

	/** @name Overridden form MatrixIdentConcat */
	//@{
	///
	Range1D D_rng() const;
	///
	Range1D I_rng() const;
	///
	value_type alpha() const;
	///
	const MatrixOp& D() const;
	///
	BLAS_Cpp::Transp D_trans() const;
	//@}

	/** @name Overridden from MatrixOp */
	//@{
	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;
	///
	/** The default just performs a shallow copy and just copies
	 * the underlying smart reference counted pointer.  If other
	 * behavior is desired then this method must be overridden.
	 */
	MatrixOp& operator=(const MatrixOp& m);
	//@}

private:

#ifdef DOXYGEN_COMPILE
	AbstractLinAlgPack::VectorSpace    *space_cols;
	AbstractLinAlgPack::VectorSpace    *space_rows;
	AbstractLinAlgPack::MatrixOp   *D;
	RangePack::Range1D                 D_rng;
	RangePack::Range1D                 I_rng;
#else
	VectorSpace::space_ptr_t  space_cols_;
	VectorSpace::space_ptr_t  space_rows_;
	D_ptr_t           D_ptr_;
	Range1D           D_rng_;
	Range1D           I_rng_;
#endif
	value_type        alpha_;
	BLAS_Cpp::Transp  D_trans_;
 
	//
	void assert_initialized() const;

}; // end class MatrixIndentConcatStd

} // end namespace ConstrainedOptPack

#endif // MATRIX_IDENT_CONCAT_STD_H
