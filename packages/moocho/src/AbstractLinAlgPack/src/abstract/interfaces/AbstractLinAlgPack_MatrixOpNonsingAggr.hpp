// /////////////////////////////////////////////////////////////////////
// MatrixWithOpNonsingularAggr.h
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

#ifndef MATRIX_WITH_OP_NONSINGULAR_AGGR_H
#define MATRIX_WITH_OP_NONSINGULAR_AGGR_H

#include "MatrixWithOpNonsingular.h"

namespace AbstractLinAlgPack {

///
/** Aggregate matrix class pulling together a \c MatrixWithOp object and a
 * \c MatrixNonsingular object into a unified matrix object.
 *
 * ToDo: Finish documentation!
 */
class MatrixWithOpNonsingularAggr
	: virtual public MatrixWithOpNonsingular
{
public:

	/** @name Public types */
	//@{

	///
	typedef MemMngPack::ref_count_ptr<const MatrixWithOp>        mwo_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const MatrixNonsingular>   mns_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	/// Construct to uninitialized
	MatrixWithOpNonsingularAggr();

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	MatrixWithOpNonsingularAggr(
		const mwo_ptr_t       &mwo
		,BLAS_Cpp::Transp     mwo_trans
		,const mns_ptr_t      &mns
		,BLAS_Cpp::Transp     mns_trans
		);

	///
	/** Initialize.
	 *
	 * ToDo: Finish documentation.
	 */
	void initialize(
		const mwo_ptr_t       &mwo
		,BLAS_Cpp::Transp     mwo_trans
		,const mns_ptr_t      &mns
		,BLAS_Cpp::Transp     mns_trans
		);

	///
	/** Set uninitialized.
	 *
	 * ToDo: Finish documentation.
	 */
	void set_uninitialized();

	//@}
	
	/** @name Access */
	//@{

	///
	const mwo_ptr_t& mwo() const;
	///
	BLAS_Cpp::Transp mwo_trans() const;
	///
	const mns_ptr_t& mns() const;
	///
	BLAS_Cpp::Transp mns_trans() const;

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
	MatrixWithOp::mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
	///
	MatrixWithOp& operator=(const MatrixWithOp& M);
	///
	std::ostream& output(std::ostream& out) const;
	///
	bool Mp_StM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs) const;
	///
	bool Mp_StMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp M_trans
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		) const;
	///
	bool Mp_StPtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		, BLAS_Cpp::Transp M_trans
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
		, const VectorWithOp& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorWithOp& v_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;
	///
	value_type transVtMtV(
		const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const VectorWithOp& v_rhs3) const;
	///
	value_type transVtMtV(
		const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;
	///
	void syr2k(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
		, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
		, value_type beta, MatrixSymWithOp* symwo_lhs ) const;
	///
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;
	///
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;
	///
	void syrk(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, value_type beta, MatrixSymWithOp* sym_lhs ) const;
	
	//@}

	/** @name Overridden from MatrixNonsingular */
	//@{

	///
	void V_InvMtV(
		VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const VectorWithOp& v_rhs2) const;
	///
	void V_InvMtV(
		VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;
	///
	value_type transVtInvMtV(
		const VectorWithOp& v_rhs1
		,BLAS_Cpp::Transp trans_rhs2, const VectorWithOp& v_rhs3) const;
	///
	value_type transVtInvMtV(
		const SpVectorSlice& sv_rhs1
		,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;
	///
	void M_StInvMtM(
		MatrixWithOp* m_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		) const;
	///
	void M_StMtInvM(
		MatrixWithOp* m_lhs, value_type alpha
		,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2
		) const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	MatrixWithOp         *mwo;
	MatrixNonsingular    *mns;
#else
	mwo_ptr_t            mwo_;
	BLAS_Cpp::Transp     mwo_trans_;
	mns_ptr_t            mns_;
	BLAS_Cpp::Transp     mns_trans_;
#endif

}; // end class MatrixWithOpNonsingularAggr

// ////////////////////////////////////
// Inline members

inline
const MatrixWithOpNonsingularAggr::mwo_ptr_t&
MatrixWithOpNonsingularAggr::mwo() const
{
	return mwo_;
}

inline
BLAS_Cpp::Transp MatrixWithOpNonsingularAggr::mwo_trans() const
{
	return mwo_trans_;
}

inline
const MatrixWithOpNonsingularAggr::mns_ptr_t&
MatrixWithOpNonsingularAggr::mns() const
{
	return mns_;
}

inline
BLAS_Cpp::Transp MatrixWithOpNonsingularAggr::mns_trans() const
{
	return mns_trans_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_WITH_OP_NONSINGULAR_AGGR_H
