// ///////////////////////////////////////////////////////////////////////////
// MatrixSymNonsingular.h
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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_SYM_NONSINGULAR_H
#define ABSTRACT_LINALG_PACK_MATRIX_SYM_NONSINGULAR_H

#include "MatrixNonsingular.h"

namespace AbstractLinAlgPack {

///
/** Abstract base class for all polymorphic symmetrix nonsingular matrices that
  * can be used to solve for linear systems relatively efficently.
  */
class MatrixSymNonsingular
	: public virtual MatrixNonsingular
{
public:

	///
	enum EMatrixDummyArg { DUMMY_ARG };

	/** @name Level-3 */
	//@{

	///
	/** symwo_lhs = alpha * op(mwo) * inv(M) * op(mwo)'.
	 *
	 * The default implementation is based on the operation M_StInvMtM(...)
	 * assuming that this #M# is a symmetric matrix.  For an efficient implementation
	 * (for this = L*L' for instance) the subclass may want to override this function.
	 */
	virtual void M_StMtInvMtM(
		MatrixSymWithOp* symwo_lhs, value_type alpha
		,const MatrixWithOp& mwo, BLAS_Cpp::Transp mwo_trans
		,EMatrixDummyArg
		) const;

	//@}

};	// end class MatrixSymNonsingular

// //////////////////////////////////////////////////////////
// Inline nonmember helper function.

inline
/// sym_gms_lhs = alpha * op(mwo) * inv(mswof) * op(mwo)'
void M_StMtInvMtM(
	MatrixSymWithOp* sym_gms_lhs, value_type alpha
	, const MatrixWithOp& mwo
	, BLAS_Cpp::Transp mwo_trans, const MatrixSymNonsingular& mswof
	, MatrixSymNonsingular::EMatrixDummyArg mwo_rhs )
{
	mswof.M_StMtInvMtM(sym_gms_lhs,alpha,mwo,mwo_trans,mwo_rhs);
}

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_SYM_NONSINGULAR_H
