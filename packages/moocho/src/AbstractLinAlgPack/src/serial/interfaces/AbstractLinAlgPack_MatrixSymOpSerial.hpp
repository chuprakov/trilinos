// //////////////////////////////////////////////////////////
// MatrixSymWithOpSerial.h
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

#ifndef SLAP_MATRIX_SYM_WITH_OP_SERIAL_H
#define SLAP_MATRIX_SYM_WITH_OP_SERIAL_H

#include "MatrixWithOpSerial.h"
#include "AbstractLinAlgPack/include/MatrixSymWithOp.h"

namespace SparseLinAlgPack {

///
/** Abstract base class for all <tt>AbstractLinAlgPack::MatrixSymWithOp</tt> objects
 * implemented in shared memory space.
 *
 * This base class does a mapping from fully abstract linear algebra to shared memory
 * linear algebra.
 *
 * These methods should not be called directly but instead should be called through
 * the line \ref MatrixSymWithOpSerial_funcs "non-member functions" that are provided.
 */
class MatrixSymWithOpSerial
	: virtual public MatrixWithOpSerial
	, virtual public AbstractLinAlgPack::MatrixSymWithOp // doxygen needs full name
{
public:

	///
	using MatrixSymWithOp::Mp_StPtMtP;

	///
	/** sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
	  *
	  * The default operation is based on <tt>this->Vp_StMtV(...)</tt> and assumes
	  * that the matrix is symmetric.  Of course, a more efficient implementation
	  * is often needed and the sublcass would like to override this.
	  */
	virtual void Mp_StPtMtP(
		sym_gms* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
		,value_type beta
		) const;

	///
	/** sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs).
	  *
	  * The default operation is based on <tt>this->Vp_StMtV(...)</tt> and assumes
	  * that the matrix is symmetric.  Of course, a more efficient implementation
	  * is often needed and the sublcass would like to override this.
	  */
	virtual void Mp_StMtMtM(
		sym_gms* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const MatrixWithOpSerial& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
		,value_type beta
		) const;

	/** @name Overridden from MatrixSymWithOp */
	//@{

	/// Must be overridden to call <tt>MatrixWithOpSerial::space_rows()</tt>
	const VectorSpace& space_rows() const;

	/// symwo_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
 	void Mp_StPtMtP(
		MatrixSymWithOp* symwo_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
		,value_type beta
		) const;

	/// symwo_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs).
	void Mp_StMtMtM(
		MatrixSymWithOp* symwo_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const MatrixWithOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
		,value_type beta
		) const;

	//@}

};	// end class MatrixSymWithOpSerial

/** \defgroup MatrixSymWithOpSerial_funcs Inline nonmeber functions for <tt>MatrixSymWithOpSerial</tt>.
  */
//@{

inline
/// sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
void Mp_StPtMtP(
	sym_gms* sym_lhs, value_type alpha
	,MatrixSymWithOpSerial::EMatRhsPlaceHolder dummy_place_holder
	,const MatrixSymWithOpSerial& M
	,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
	,value_type beta = 1.0
	)
{
	M.Mp_StPtMtP(sym_lhs,alpha,dummy_place_holder,gpms_rhs,gpms_rhs_trans,beta);
}

inline
/// sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs) + beta * sym_lhs
void Mp_StMtMtM(
	sym_gms* sym_lhs, value_type alpha
	,MatrixSymWithOpSerial::EMatRhsPlaceHolder dummy_place_holder
	,const MatrixSymWithOpSerial& M
	,const MatrixWithOpSerial& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
	,value_type beta = 1.0
	)
{
	M.Mp_StMtMtM(sym_lhs,alpha,dummy_place_holder,mwo_rhs,mwo_rhs_trans,beta);
}

//@}

}	// end namespace SparseLinAlgPack 

#endif	// SLAP_MATRIX_SYM_WITH_OP_SERIAL_H
