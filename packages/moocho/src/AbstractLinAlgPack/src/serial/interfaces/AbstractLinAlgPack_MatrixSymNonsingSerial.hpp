// ///////////////////////////////////////////////////////////////////////////
// MatrixSymNonsingularSerial.h
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

#ifndef SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H
#define SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H

#include "MatrixNonsingularSerial.h"
#include "AbstractLinAlgPack/include/MatrixSymNonsingular.h"

namespace SparseLinAlgPack {

///
/** Abstract base class for all serial polymorphic symmetrix nonsingular matrices that
 * can be used to solve for linear systems relatively efficiently.
 *
 * The methods of this interface should not be called directly but instead through
 * the \ref MatrixSymNonsingularSerial_funcs "provided nonmember functions".
 */
class MatrixSymNonsingularSerial
	: virtual public MatrixNonsingularSerial
	, virtual public AbstractLinAlgPack::MatrixSymNonsingular // doxygen needs full name
{
public:

	///
	using MatrixSymNonsingular::M_StMtInvMtM;

	/** @name Level-3 */
	//@{

	///
	/** sym_gms_lhs = alpha * op(mwo) * inv(M) * op(mwo)'.
	  *
	  * The default implementation is based on the operation M_StInvMtM(...)
	  * assuming that this \c M is a symmetric matrix.  For an efficient implementation
	  * (for this = L*L' for instance) the subclass may want to override this function.
	  */
	virtual void M_StMtInvMtM(
		sym_gms* sym_gms_lhs, value_type alpha
		,const MatrixWithOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
		,EMatrixDummyArg
		) const;

	//@}

	/** @name Overridden from MatrixSymNonsingular */
	//@{

	void M_StMtInvMtM(
		MatrixSymWithOp* sym_lhs, value_type alpha
		,const MatrixWithOp& mwo, BLAS_Cpp::Transp mwo_trans
		,EMatrixDummyArg
		) const;

	//@}

};	// end class MatrixSymNonsingularSerial

/** \defgroup MatrixSymNonsingularSerial_funcs MatrixSymNonsingularSerial nonmember inline functions.
 *
 * These nonmember functions allow operations to be called on \c MatrixSymNonsingularSerial objects
 * in similar manner to those in \c LinAlgPack.
 */
//@{

inline
/// sym_gms_lhs = alpha * op(mwo) * inv(mswof) * op(mwo)'
void M_StMtInvMtM(
	sym_gms* sym_gms_lhs, value_type alpha
	,const MatrixWithOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
	,const MatrixSymNonsingularSerial& mswons
	,MatrixSymNonsingularSerial::EMatrixDummyArg mwo_rhs
	 )
{
	mswons.M_StMtInvMtM(sym_gms_lhs,alpha,mwo,mwo_trans,mwo_rhs);
}

//@}

} // end namespace SparseLinAlgPack

#endif	// SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H
