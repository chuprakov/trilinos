// ////////////////////////////////////////////////////////////////////////
// MatrixIdentConcatStd.h
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

#include "MatrixIdentConcat.h"
#include "LinAlgPack/include/Range1D.h"
#include "Misc/include/ref_count_ptr.h"

namespace ConstrainedOptimizationPack {

///
/** Concrete implementation class for a matrix vertically concatonated with an identity matrix.
 *
 * Represents an interface for a matrix that represents:
 * #M = [ alpha*op(D); I ] (TOP)#
 * or
 * #M = [ I; alpha*op(D) ] (BOTTOM)#
 */
class MatrixIdentConcatStd : public MatrixIdentConcat {
public:

	/** @name Setup and representation access. */
	//@{
	///
	enum ETopBottom { TOP, BOTTOM };
	///
	typedef ReferenceCountingPack::ref_count_ptr<const MatrixWithOp> D_ptr_t;
	///
	/** Setup with a matrix object.
	 *
	 * @param  top_or_bottom
	 *                 [in] If #TOP# then #M = [ alpha*op(D); I ]# and if #BOTTOM# then #M = [ I; alpha*op(D) ]#
	 * @param  alpha   [in] Scalar that multiplies op(D).
	 * @param  D_ptr   [in] Smart pointer to a matrix object that represents #D#.  The matrix object pointed to must
	 *                 not be altered until #this# object is no longer in use or #this->set_uninitialized() has been called.
	 *                 Here #D.get()# must point to a valid matrix object.  If #D.get() == NULL# then an #std::invalid_argument#
	 *                 exception will be thrown.
	 * @param  D_trans [in] Determines if #op(D) = D# (#no_trans#) or #D'# (#trans#).
	 */
	virtual void initialize(ETopBottom top_or_bottom, value_type alpha, const D_ptr_t& D_ptr, BLAS_Cpp::Transp D_trans);
	///
	/** Set the matrix to uninitialized.
	 */
	virtual void set_uninitialized();
	///
	virtual const D_ptr_t& D_ptr() const;
	//@}

	/** @name Overridden form MatrixIdentConcat. */
	//@{
	///
	Range1D D_rng() const;
	///
	Range1D I_rng() const;
	///
	value_type alpha() const;
	///
	const MatrixWithOp& D() const;
	///
	BLAS_Cpp::Transp D_trans() const;
	//@}

	/** @name Overridden from MatrixWithOp. */
	//@{
	///
	/** The default just performs a shallow copy and just copies
	 * the underlying smart reference counted pointer.  If other
	 * behavior is desired then this method must be overridden.
	 */
	MatrixWithOp& operator=(const MatrixWithOp& m);
	//@}

private:
	value_type        alpha_;
	D_ptr_t           D_ptr_;
	BLAS_Cpp::Transp  D_trans_;
	Range1D           D_rng_;
	Range1D           I_rng_;
 
}; // end class MatrixIndentConcatStd

} // end namespace ConstrainedOptimizationPack

#endif // MATRIX_IDENT_CONCAT_STD_H
