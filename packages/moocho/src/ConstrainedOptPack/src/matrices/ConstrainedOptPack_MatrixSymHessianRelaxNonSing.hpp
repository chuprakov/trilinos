// /////////////////////////////////////////////////////////////////////
// MatrixSymHessianRelaxNonSing.h
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

#ifndef MATRIX_SYM_HESSIAN_RELAX_NON_SING_H
#define MATRIX_SYM_HESSIAN_RELAX_NON_SING_H

#include "ConstrainedOptimizationPackTypes.h"
#include "AbstractLinAlgPack/src/MatrixSymDiagonalStd.h"
#include "ref_count_ptr.h"

namespace ConstrainedOptimizationPack {

///
/** Matrix class for non-singular Hessian matrix augmented with a terms for
 * "Big M" relaxation variables.
 *
 * The matrix that is formed is:
 \begin{verbatim}

 H = [ G    ]
     [    M ]
 \end{verbatim}
 * where <tt>M</tt> is a diagonal matrix made up of entries <tt>M_diag</tt>.
 */
class MatrixSymHessianRelaxNonSing
	: public AbstractLinAlgPack::MatrixSymWithOpNonsingular
{
public:

	///
	typedef MemMngPack::ref_count_ptr<const MatrixSymWithOpNonsingular>  G_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<VectorWithOpMutable>               vec_mut_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const VectorSpace>                 space_ptr_t;
	
	/** @name Constructors/initializers */
	//@{

	///
	/** Construct uninitialized.
	 *
	 * Postconditions:
	 * <ul>
	 * <li> <tt>this->rows() == 0</tt>
	 * <li> <tt>&this->G_ptr().get() == NULL</tt>
	 * <li> <tt>&this->M_diag_ptr().get() == NULL</tt>
	 * </ul>
	 */
	MatrixSymHessianRelaxNonSing();

	///
	/** Constructor (calls \c initialize()).
	 */
	MatrixSymHessianRelaxNonSing(
		const G_ptr_t         &G_ptr
		,const vec_mut_ptr_t  &M_diag_ptr
		,const space_ptr_t    &space = MemMngPack::null
		);

	///
	/** Initialize the Hessian and the relaxation terms.
	 *
	 * Preconditions:
	 * <ul>
	 * <li> <tt>G_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>G_ptr->rows() > 0</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>M_diag_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>M_diag_ptr->dim() > 0</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:
	 * <ul>
	 * <li> <tt>this->rows() == G_ptr->rows() + M_diag->dim()</tt>
	 * <li> <tt>&this->G() == G_ptr.get()</tt>
	 * <li> <tt>&this->M().diag() == M_diag_ptr.get()</tt>
	 * </ul>
	 *
	 * @param  G_ptr   [in] Smart pointer to matrix that this object will represent.  The underlying
	 *                 matrix object <tt>*G_ptr.get()</tt> should not be modified without calling \c this->initialize()
	 *                 again.
	 * @param  M_diag_ptr
	 *                 [in] Smart pointer to the diagonal for <tt>M</tt>.  All of the elements in this vector must
	 *                 be nonzero!  This vector is used to initalize the diagonal matrix <tt>M</tt> and this vector should
	 *                 not be modified externally without calling \c this->initialize() again.
	 * @param  space   [in] Smart pointer to the space used for <tt>this->space_cols()</tt> and <tt>this->space_rows()</tt>.
	 *                 If <tt>space.get() == NULL</tt> then <tt>VectorSpaceCompoisteStd</tt> will be used which is
	 *                 constructed from the spaces <tt>G_ptr->space_cols()</tt> and <tt>M_diag_ptr->space()</tt>.
	 */
	void initialize(
		const G_ptr_t         &G_ptr
		,const vec_mut_ptr_t  &M_diag_ptr
		,const space_ptr_t    &space = MemMngPack::null
		);

	///
	/** Set uninitalized.
	 *
	 */
	void set_uninitialized();
	
	///
	const G_ptr_t& G_ptr() const;

	///
	const vec_mut_ptr_t& M_diag_ptr() const;

	///
	const MatrixSymWithOpNonsingular& G() const;

	///
	const SparseLinAlgPack::MatrixSymDiagonalStd& M() const;

	//@}
	
	/** @name Overridden from MatrixWithOp */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	bool Mp_StM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const VectorWithOp& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* v_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_rhs2_trans
		,const VectorWithOp& v_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* v_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_rhs2_trans
		,const SpVectorSlice& sv_rhs3, value_type beta) const;

	//@}

	/** @name Overridden form MatrixSymWithOp */
	//@{

	///
	void Mp_StPtMtP(
		MatrixSymWithOp* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
		,value_type beta
		) const;

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
		,const SpVectorSlice& sv_rhs2) const;

	//@}

private:
	
	// ///////////////////////////////
	// Private data members

	space_ptr_t              vec_space_;
	G_ptr_t                  G_ptr_;
	MatrixSymDiagonalStd     M_;

	// ///////////////////////////////
	// Private member functions

	void assert_initialized() const;

};

// ////////////////////////////////////
// Inline members

inline
const MatrixSymHessianRelaxNonSing::G_ptr_t&
MatrixSymHessianRelaxNonSing::G_ptr() const
{
	return G_ptr_;
}

inline
const MatrixSymHessianRelaxNonSing::vec_mut_ptr_t&
MatrixSymHessianRelaxNonSing::M_diag_ptr() const
{
	return M_.diag_ptr();
}

inline
const MatrixSymWithOpNonsingular&
MatrixSymHessianRelaxNonSing::G() const
{
	assert_initialized();
	return *G_ptr_;
}

inline
const MatrixSymDiagonalStd&
MatrixSymHessianRelaxNonSing::M() const
{
	assert_initialized();
	return M_;
}
	
} // end namespace ConstrainedOptimizationPack

#endif // MATRIX_SYM_HESSIAN_RELAX_NON_SING_H
