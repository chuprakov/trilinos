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
#include "MatrixSymWithOpFactorized.h"
#include "VectorSpace.h"

namespace AbstractLinAlgPack {

///
/** Simple diagonal matrix class.
 */
class MatrixSymDiagonalStd
	: public virtual MatrixSymInitDiagonal
	, public virtual MatrixSymWithOpFactorized
{
public:

	/// Initialized with the diagonal vector NULL.
	MatrixSymDiagonalStd();

	/// Give access to the diagonal VectorSlice
	VectorWithOpMutable& diag();
	///
	const VectorWithOp& diag() const;

	/** @name Overridden from MatrixBase */
	//@{

	///
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
	/** Add to a mutable matrix lhs.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #dynamic_cast<MatrixWithOpMutable*>(m_lhs) != NULL# (throw ???).
	 * \end{itemize}
	 */
	void Mp_StM(MatrixWithOp* g_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const;
	///
	void Vp_StMtV(VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorWithOp& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;

	//@}

	/** Overridden from MatrixWithOpFactorized */
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

}; // end class MatrixSymInitDiagonal

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SYM_DIAGONAL_STD_H
