// //////////////////////////////////////////////////////////////////
// MatrixZero.h
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

#ifndef ALAP_MATRIX_ZERO_H
#define ALAP_MATRIX_ZERO_H

#include "MatrixWithOp.h"
#include "VectorSpace.h"

namespace AbstractLinAlgPack {

///
/** Implementation of a matrix with all zeros.
 *
 * This may seem like a silly class but it is helpful in some circumstances.
 * This class needs to be updated whenever methods are added or removed
 * from \c MatrixWithOp.
 */
class MatrixZero : public MatrixWithOp {
public:

	/** @name Constructors/initializers */
	//@{

	/// Calls <tt>this->initalize()</tt>
	MatrixZero(
		const VectorSpace::space_ptr_t&     space_cols = VectorSpace::space_ptr_t(NULL)
		,const VectorSpace::space_ptr_t&    space_rows = VectorSpace::space_ptr_t(NULL)
		);

	///
	/** Initialize (or initialize) given the columns and rows vector spaces.
	 *
	 * Preconditions:<ul>
	 * <li> [<tt>space_cols.get() == NULL</tt>] <tt><tt>space_rows.get() == NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> [<tt>space_cols.get() != NULL</tt>] <tt><tt>space_rows.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>space_cols.get() == NULL</tt>] <tt>this->rows() == 0</tt>
	 * <li> [<tt>space_cols.get() == NULL</tt>] <tt>this->cols() == 0</tt>
	 * <li> [<tt>space_cols.get() != NULL</tt>] <tt>&this->space_cols() == space_cols.get()</tt>
	 * <li> [<tt>space_cols.get() != NULL</tt>] <tt>&this->space_rows() == space_rows.get()</tt>
	 * </ul>
	 */
	void initialize(
		const VectorSpace::space_ptr_t&    space_cols
		,const VectorSpace::space_ptr_t&   space_rows
		);

	//@}

	/** @name Overridden from MatrixBase */
	//@{
	///
	size_type rows() const;
	///
	size_type cols() const;
	///
	size_type nz() const;
	//@}

	/** @name Overridden from MatrixWithOp */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;
	///
	void zero_out();
	///
	void Mt_S(value_type alpha);
	///
	MatrixWithOp& operator=(const MatrixWithOp& M);
	///
	std::ostream& output(std::ostream& out) const;
	///
	bool Mp_StM(
		MatrixWithOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs
		) const;
	///
	bool Mp_StMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		) const;
	///
	bool Mp_StPtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		,BLAS_Cpp::Transp M_trans
		) const;
	///
	bool Mp_StPtMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const VectorWithOp& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_rhs2_trans
		,const VectorWithOp& v_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_rhs2_trans
		,const SpVectorSlice& sv_rhs3, value_type beta) const;
	///
	value_type transVtMtV(
		const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2
		,const VectorWithOp& v_rhs3) const;
	///
	value_type transVtMtV(
		const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;
	///
	void syr2k(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		,const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
		,const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
		,value_type beta, MatrixSymWithOp* symwo_lhs ) const;
	///
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		,value_type beta ) const;
	///
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2
		,value_type beta ) const;
	///
	bool syrk(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		,value_type beta, MatrixSymWithOp* sym_lhs ) const;

	//@}

private:

	VectorSpace::space_ptr_t  space_cols_;
	VectorSpace::space_ptr_t  space_rows_;

	//
	void assert_initialized() const;

}; // end class MatrixZero

} // end namespace AbstractLinAlgPack

#endif // ALAP_MATRIX_ZERO_H
