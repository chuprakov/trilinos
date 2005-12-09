// //////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymOp.hpp
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

#include "AbstractLinAlgPack_MatrixOp.hpp"

namespace AbstractLinAlgPack {

///
/** Interface adding operations specific for a symmetric matrix {abstract}.
 *
 * This interface defines two addition methods to those found in \c MatrixOp:
 *
 * <tt>sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs</tt><br>
 * <tt>sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs) + beta * sym_lhs</tt><br>
 *
 * The reason that these methods could not be defined in the \c MatrixOp interface
 * is that the lhs matrix matrix argument \c sym_lhs is only guaranteed to be
 * symmetric if the rhs matrix argument \c M (which is \c this matrix) is guaranteed
 * to be symmetric.  Since a \c MatrixOp matrix object may be unsymmetric (as
 * well as rectangular), it can not implement this operation, only a symmetric
 * matrix can.
 *
 * Clients should use the \ref MatrixSymWithOp_funcs_grp "provided non-member functions"
 * to call the methods and not the methods themselves.
 */
class MatrixSymOp : public virtual MatrixOp {
public:

	///
	using MatrixOp::Mp_StPtMtP;

	/** @name Public types */
	//@{

#ifndef DOXYGEN_COMPILE
	///
	typedef Teuchos::RefCountPtr<const MatrixSymOp>    mat_mswo_ptr_t;
	///
	typedef Teuchos::RefCountPtr<MatrixSymOp>          mat_mswo_mut_ptr_t;
#endif
	///
	enum EMatRhsPlaceHolder { DUMMY_ARG };

	//@}

	/** @Friends */
	//@{

	///
	friend
	void Mp_StPtMtP(
		MatrixSymOp* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const MatrixSymOp& M
		,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
		,value_type beta
		);
	///
	friend
	void Mp_StMtMtM(
		MatrixSymOp* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const MatrixSymOp& M
		,const MatrixOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
		,value_type beta
		);

	//@}

	/** @name Clone */
	//@{

	///
	/** Clone the non-const matrix object (if supported).
	 *
	 * The default implementation returns NULL which is perfectly acceptable.
	 * A matrix object is not required to return a non-NULL value but almost
	 * every good matrix implementation will.
	 */
	virtual mat_mswo_mut_ptr_t clone_mswo();

	///
	/** Clone the const matrix object (if supported).
	 *
	 * The behavior of this method is the same as for the non-const version
	 * above except it returns a smart pointer to a const matrix object.
	 */
	virtual mat_mswo_ptr_t clone_mswo() const;

	//@}

protected:

	/** @name Level-1 BLAS */
	//@{

	///
	/** sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
	  *
	  * The default operation is based on Vp_StMtV(...) and assumes
	  * that the matrix is symmetric.  Of course, a more efficient implementation
	  * is often needed and the sublcass would like to override this.
	  */
	virtual void Mp_StPtMtP(
		MatrixSymOp* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
		,value_type beta
		) const;

	//@}

	/** @name Level-3 BLAS */
	//@{

	///
	/** sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs) + beta * sym_lhs.
	  *
	  * The default operation is based on \c Vp_StMtV() and assumes
	  * that the matrix is symmetric.  Of course, a more efficient implementation
	  * is often needed and the sublcass would like to override this.
	  */
	virtual void Mp_StMtMtM(
		MatrixSymOp* sym_lhs, value_type alpha
		,EMatRhsPlaceHolder dummy_place_holder
		,const MatrixOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
		,value_type beta
		) const;

	//@}

public:

	/** Overridden from MatrixOp */
	//@{
	/// Returns <tt>this->rows()</tt>
	size_type cols() const;
	// Returns <tt>this->space_cols()</tt>
	const VectorSpace& space_rows() const;
	/// Returns <tt>this->clone_mswo()</tt>.
	mat_mut_ptr_t clone();
	/// Returns <tt>this->clone_mswo()</tt>.
	mat_ptr_t clone() const;
	//@}

	/// Calls operator=(MatrixOp&)
	virtual MatrixSymOp& operator=(const MatrixSymOp& M)
	{ static_cast<MatrixOp*>(this)->operator=(M); return *this; }

};	// end class MatrixSymOp

/** \defgroup MatrixSymWithOp_funcs_grp Inline nonmeber functions for MatrixSymOp to call methods.
  *
  * These allow nonmember functions to act like virtual functions.
  */
//@{

inline
/// sym_lhs = alpha * op(gpms_rhs') * M * op(gpms_rhs) + beta * sym_lhs.
void Mp_StPtMtP(
	MatrixSymOp* sym_lhs, value_type alpha
	,MatrixSymOp::EMatRhsPlaceHolder dummy_place_holder
	,const MatrixSymOp& M
	,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
	,value_type beta = 1.0
	)
{
	M.Mp_StPtMtP(sym_lhs,alpha,dummy_place_holder,gpms_rhs,gpms_rhs_trans,beta);
}

inline
/// sym_lhs = alpha * op(mwo_rhs') * M * op(mwo_rhs) + beta * sym_lhs
void Mp_StMtMtM(
	MatrixSymOp* sym_lhs, value_type alpha
	,MatrixSymOp::EMatRhsPlaceHolder dummy_place_holder
	,const MatrixSymOp& M
	,const MatrixOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
	,value_type beta = 1.0
	)
{
	M.Mp_StMtMtM(sym_lhs,alpha,dummy_place_holder,mwo_rhs,mwo_rhs_trans,beta);
}

//@}

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_SYM_WITH_OP_H
