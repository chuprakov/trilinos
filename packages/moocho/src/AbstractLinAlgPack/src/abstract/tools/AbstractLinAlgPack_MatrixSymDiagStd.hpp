// ///////////////////////////////////////////
// MatrixSymDiagonalStd.h
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

#ifndef MATRIX_SYM_DIAGONAL_STD_H
#define MATRIX_SYM_DIAGONAL_STD_H

#include "MatrixSymInitDiagonal.h"
#include "MatrixSymWithOpNonsingular.h"
#include "VectorSpace.h"

namespace AbstractLinAlgPack {

///
/** Simple diagonal matrix class.
 *
 * ToDo: Implement clone_mswons() and deal with this->diag_ptr().count() > 1
 * by cloning vector if told to.  This allows lazy evaluation of the clone_mswons()
 * method.
 */
class MatrixSymDiagonalStd
	: public virtual MatrixSymInitDiagonal
	, public virtual MatrixSymWithOpNonsingular
{
public:

	/** @name Constructors/initalizers */
	//@{

	/// Calls <tt>this->initialize()</tt>.
	MatrixSymDiagonalStd(
		const VectorSpace::vec_mut_ptr_t& diag   = MemMngPack::null
		,bool                             unique = true
		);

	///
	/** Initialize given the diagonal vector (or no vector at all).
	 *
	 * @param  diag   [in] Vector to be used for the diagonal.  If <tt>diag.get() == NULL</tt>
	 *                then \c this will be uninitialized.
	 * @param  unique [in] Determines if the underlying \c diag vector is guaranteed to be
	 *                unique and not shared.
	 */
	void initialize(
		const VectorSpace::vec_mut_ptr_t& diag
		,bool                             unique = true
		);

	//@}

	/** @name Access */
	//@{

	///
	/** Give non-const access to the diagonal vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->diag_ptr().get() != NULL</tt> (throw <tt>???</tt>)
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	VectorWithOpMutable& diag();
	///
	const VectorWithOp& diag() const;
	///
	const VectorSpace::vec_mut_ptr_t& diag_ptr() const;
	///
	bool unique() const;

	//@}

	/** @name Overridden from MatrixBase */
	//@{

	/// Returns 0 if not initalized (this->diag() == NULL).
	size_type rows() const;
	///
	size_type nz() const;

	//@}

	/** @name Overridden from MatrixWithOp */
	//@{

	///
	const VectorSpace& space_rows() const;
	///
	const VectorSpace& space_cols() const;
	///
	MatrixWithOp& operator=(const MatrixWithOp& mwo_rhs);
	///
	/** Add to a mutable matrix lhs.
	 *
	 * Preconditions:<ul>
	 * <li> #dynamic_cast<MultiVectorMutable*>(m_lhs) != NULL#.
	 * </ul>
	 */
	bool Mp_StM(MatrixWithOp* g_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const;
	///
	void Vp_StMtV(VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorWithOp& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;

	//@}

	/** Overridden from MatrixWithOpNonsingular */
	//@{

	///
	void V_InvMtV(VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorWithOp& v_rhs2) const;
	///
	void V_InvMtV(VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;

	//@}

	/** @name Overridden from MatrixSymInitDiagonal */
	//@{

	///
	void init_identity( const VectorSpace& space_diag, value_type alpha );
	///
	void init_diagonal( const VectorWithOp& diag );

	//@}

private:

	VectorSpace::vec_mut_ptr_t     diag_;
	bool                           unique_;

	void copy_unique();

}; // end class MatrixSymDiagonalStd

// ////////////////////////////////////////
// Inline members

inline
bool MatrixSymDiagonalStd::unique() const
{
	return unique_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SYM_DIAGONAL_STD_H
