// //////////////////////////////////////////////////////////
// MatrixSymWithOp.h
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

#ifndef MATRIX_SYM_WITH_OP_H
#define MATRIX_SYM_WITH_OP_H

#include "MatrixWithOp.h"

namespace AbstractLinAlgPack {

///
/** Interface adding operations specific for a symmetric matrix {abstract}.
 *
 *
 * Clients should use the \ref MatrixSymWithOp_funcs_grp "provided non-member functions"
 * to call the methods and not the methods themselves.
 *
 * ToDo: Finish documentation!
 */
class MatrixSymWithOp : public virtual MatrixWithOp {
public:

	///
	using MatrixWithOp::Mp_StPtMtP;

	///
	enum EMatRhsPlaceHolder { DUMMY_ARG };

	/// Returns size equal to rows() (need not be overridden by subclass again)
	size_type cols() const;

	///
	/** sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
	  *
	  * The default operation is based on Vp_StMtV(...) and assumes
	  * that the matrix is symmetric.  Of course, a more efficient implementation
	  * is often needed and the sublcass would like to override this.
	  */
	virtual void Mp_StPtMtP(
		MatrixSymWithOp* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
		,value_type beta = 1.0
		) const;

	///
	/** sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs).
	  *
	  * The default operation is based on Vp_StMtV(...) and assumes
	  * that the matrix is symmetric.  Of course, a more efficient implementation
	  * is often needed and the sublcass would like to override this.
	  */
	virtual void Mp_StMtMtM(
		MatrixSymWithOp* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const MatrixWithOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
		,value_type beta = 1.0
		) const;

};	// end class MatrixSymWithOp

/** \defgroup MatrixSymWithOp_funcs_grp Inline nonmeber functions for MatrixSymWithOp to call methods.
  */
//@{

inline
/// sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
void Mp_StPtMtP(
	MatrixSymWithOp* sym_lhs, value_type alpha
	,MatrixSymWithOp::EMatRhsPlaceHolder dummy_place_holder
	,const MatrixSymWithOp& M
	,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
	,value_type beta = 1.0
	)
{
	M.Mp_StPtMtP(sym_lhs,alpha,dummy_place_holder,gpms_rhs,gpms_rhs_trans,beta);
}

inline
/// sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs) + beta * sym_lhs
void Mp_StMtMtM(
	MatrixSymWithOp* sym_lhs, value_type alpha
	,MatrixSymWithOp::EMatRhsPlaceHolder dummy_place_holder
	,const MatrixSymWithOp& M
	,const MatrixWithOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
	,value_type beta = 1.0
	)
{
	M.Mp_StMtMtM(sym_lhs,alpha,dummy_place_holder,mwo_rhs,mwo_rhs_trans,beta);
}

//@}

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_SYM_WITH_OP_H
