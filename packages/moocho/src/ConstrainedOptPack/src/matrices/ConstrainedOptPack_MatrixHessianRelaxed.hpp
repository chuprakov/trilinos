// /////////////////////////////////////
// ConstrainedOptPack_MatrixHessianRelaxed.hpp
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

#ifndef MATRIX_HESSIAN_RELAXED_H
#define MATRIX_HESSIAN_RELAXED_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymOp.hpp"

namespace ConstrainedOptPack {

///
/** Represents a symmetric Hessian matrix with a relaxation variable
  * added.
  *
  * This class is used to represent the matrix:
  \begin{verbatim}
	     [ H       ]
	G =  [    bigM ]
  \end{verbatim}
  *
  */
class MatrixHessianRelaxed : public MatrixSymOp {
public:

	/// Construct to uninitialized
	MatrixHessianRelaxed();

	///
	/** Initialize.
	  *
	  * ToDo: Finish documentation!
	  *
	  */
	void initialize(
		  const MatrixSymOp	&H
		, value_type			bigM
		);

	// ///////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	// //////////////////////////////
	// Overridden from MatrixOp

	///
	void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const DVectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const DVectorSlice& vs_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;
	///
	value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const ;

private:
	size_type				n_;	// size of H
	const MatrixSymOp	*H_;
	value_type				bigM_;

};	// end class MatrixHessianRelaxed

}	// end namespace ConstrainedOptPack

#endif 	// MATRIX_HESSIAN_RELAXED_H
