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
#include "SparseLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "Misc/include/ref_count_ptr.h"

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
 * where #M# is a diagonal matrix made up of entries M_diag
 *
 */
class MatrixSymHessianRelaxNonSing
	: public MatrixSymWithOpFactorized
{
public:

	///
	typedef MemMngPack::ref_count_ptr<const MatrixSymWithOpFactorized>  G_ptr_t;
	
	///
	/** Construct
	 *
	 * Calls initialize(...).
	 */
	MatrixSymHessianRelaxNonSing(
		const G_ptr_t       &G_ptr  = G_ptr_t(NULL)
		,const VectorSlice  &M_diag = VectorSlice()
		);

	///
	/** Initialize the Hessian and the relaxation terms.
	 *
	 * Preconditions:
	 * \begin{itemize}
	 * \item [#G_ptr.get() != NULL#] #M_diag.size() >= 1# (throw #std::invalid_argument#)
	 * \item [#G_ptr.get() != NULL#] #G_ptr->rows() >= 1# (throw #std::invalid_argument#)
	 * \item [#G_ptr.get() == NULL#] #M_diag.size() == 0# (throw #std::invalid_argument#)
	 * \end{itemize}
	 *
	 * Postconditions:
	 * \begin{itemize}
	 * \item [#G_ptr.get() != NULL#] #this->rows() == G_ptr->rows() + M_diag.size()#
	 * \item [#G_ptr.get() != NULL#] #&this->G() == G_ptr.get()#
	 * \item [#G_ptr.get() != NULL#] #this->M().diag() == M_diag#
	 * \item [#G_ptr.get() == NULL#] #this->rows() == 0#
	 * \end{itemize}
	 *
	 * @param  G_ptr   [in] Smart pointer to matrix that this object will represent.  The underlying
	 *                 matrix object #*G_ptr.get()# should not be modified without calling initialize(...)
	 *                 again.  It is allowed for #G_ptr.get() == NULL# in which case #this# will become
	 *                 uninitalized (i.e. #this->rows() == 0#).
	 * @param  M_diag  [in] Diagonal for #M#.  All of the elements in this vector must be nonzero!
	 *                 The elements of this vector are copied so no worries.
	 */
	void initialize(
		const G_ptr_t       &G_ptr
		,const VectorSlice  &M_diag = VectorSlice()
		);
	
	///
	const G_ptr_t& G_ptr() const;

	///
	const MatrixSymWithOpFactorized& G() const;

	///
	const SparseLinAlgPack::MatrixSymDiagonalStd& M() const;
	
	// /////////////////////////////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	// /////////////////////////////////////////////////////
	// Overridden from MatrixWithOp

	///
	void Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorSlice& vs_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;

	// /////////////////////////////////////////////////////
	// Overridden form MatrixSymWithOp

	void Mp_StPtMtP( sym_gms* sym_lhs, value_type alpha
		, EMatRhsPlaceHolder dummy_place_holder
		, const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
		, value_type beta ) const;

	// /////////////////////////////////////////////////////
	// Overridden from MatrixWithOpFactorized

	///
	void V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;
	///
	void V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;

private:
	
	// ///////////////////////////////
	// Private data members

	G_ptr_t                                 G_ptr_;
	SparseLinAlgPack::MatrixSymDiagonalStd  M_;

	// ///////////////////////////////
	// Private member functions

	void assert_initialized() const;

};

} // end namespace ConstrainedOptimizationPack

#endif // MATRIX_SYM_HESSIAN_RELAX_NON_SING_H
