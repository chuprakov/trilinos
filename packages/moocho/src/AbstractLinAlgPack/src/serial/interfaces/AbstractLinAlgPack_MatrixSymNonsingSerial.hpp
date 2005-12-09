// ///////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymNonsingSerial.hpp
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

#include "AbstractLinAlgPack_MatrixNonsingSerial.hpp"
#include "AbstractLinAlgPack_MatrixSymNonsing.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract base class for all serial polymorphic symmetrix nonsingular matrices that
 * can be used to solve for linear systems relatively efficiently.
 *
 * The methods of this interface should not be called directly but instead through
 * the \ref MatrixSymNonsingularSerial_funcs "provided nonmember functions".
 */
class MatrixSymNonsingSerial
	: virtual public MatrixNonsingSerial
	, virtual public AbstractLinAlgPack::MatrixSymNonsing // doxygen needs full name
{
public:

	///
	using MatrixSymNonsing::M_StMtInvMtM;

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
		DMatrixSliceSym* sym_gms_lhs, value_type alpha
		,const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
		,EMatrixDummyArg
		) const;

	//@}

	/** @name Overridden from MatrixSymNonsing */
	//@{

	void M_StMtInvMtM(
		MatrixSymOp* sym_lhs, value_type alpha
		,const MatrixOp& mwo, BLAS_Cpp::Transp mwo_trans
		,EMatrixDummyArg
		) const;

	//@}

};	// end class MatrixSymNonsingSerial

/** \defgroup MatrixSymNonsingularSerial_funcs MatrixSymNonsingSerial nonmember inline functions.
 *
 * These nonmember functions allow operations to be called on \c MatrixSymNonsingSerial objects
 * in similar manner to those in \c DenseLinAlgPack.
 */
//@{

inline
/// sym_gms_lhs = alpha * op(mwo) * inv(mswof) * op(mwo)'
void M_StMtInvMtM(
	DMatrixSliceSym* sym_gms_lhs, value_type alpha
	,const MatrixOpSerial& mwo, BLAS_Cpp::Transp mwo_trans
	,const MatrixSymNonsingSerial& mswons
	,MatrixSymNonsingSerial::EMatrixDummyArg mwo_rhs
	 )
{
	mswons.M_StMtInvMtM(sym_gms_lhs,alpha,mwo,mwo_trans,mwo_rhs);
}

//@}

} // end namespace AbstractLinAlgPack

#endif	// SLAP_MATRIX_SYM_NONSINGULAR_SERIAL_H
