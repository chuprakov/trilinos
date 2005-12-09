// ///////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymNonsing.hpp
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

#include "AbstractLinAlgPack_MatrixNonsing.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract base class for all polymorphic symmetrix nonsingular matrices that
 * can be used to solve for linear systems relatively efficently.
 *
 * This interface defines a single addition method to those found in \c MatrixNonsing:
 *
 * <tt>symwo_lhs = alpha * op(mwo) * inv(M) * op(mwo)'</tt><br>
 *
 * The reason that this method could not be defined in the \c MatrixNonsing interface
 * is that the lhs matrix matrix argument \c symwo_lhs is only guaranteed to be
 * symmetric if the rhs matrix argument \c M (which is \c this matrix) is guaranteed
 * to be symmetric.  Since a \c MatrixNonsing matrix object may be unsymmetric, it
 * can not implement this operation, only a symmetric nonsingular matrix can.
 *
 * Any symmetric nonsingular matrix abstraction that can be used to solve for nonlinear
 * systems should also be able to support the \c MatrixSymOp interface.
 * Therefore, this interface is more of an implementation artifact than
 * an a legitimate domain abstraction.  However, some symmetric linear solvers that
 * can implement this interface, can not easily implement the <tt>%MatrixSymOp</tt>
 * interface and therefore this interface is justified.  A general client should never
 * use this interface directly.  Instead, the combined interface \c MatrixSymOpNonsing
 * should be used with fully formed symmetric matrix abstractions.
 *
 * Clients should use the \ref MatrixSymNonsingular_func_grp "provided non-member functions"
 * to call the methods and not the methods themselves.
 */
class MatrixSymNonsing
	: public virtual MatrixNonsing
{
public:

	/** @name Public types */
	//@{

#ifndef DOXYGEN_COMPILE
	///
	typedef Teuchos::RefCountPtr<const MatrixSymNonsing>    mat_msns_ptr_t;
	///
	typedef Teuchos::RefCountPtr<MatrixSymNonsing>          mat_msns_mut_ptr_t;
#endif
	///
	enum EMatrixDummyArg { DUMMY_ARG };

	//@}

	/** @name Friends */
	//@{

	///
	friend
	void M_StMtInvMtM(
		MatrixSymOp* sym_gms_lhs, value_type alpha
		,const MatrixOp& mwo
		,BLAS_Cpp::Transp mwo_trans, const MatrixSymNonsing& mswof
		,EMatrixDummyArg mwo_rhs
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
	virtual mat_msns_mut_ptr_t clone_msns();

	///
	/** Clone the const matrix object (if supported).
	 *
	 * The behavior of this method is the same as for the non-const version
	 * above except it returns a smart pointer to a const matrix object.
	 */
	virtual mat_msns_ptr_t clone_msns() const;

	//@}

protected:

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
		MatrixSymOp* symwo_lhs, value_type alpha
		,const MatrixOp& mwo, BLAS_Cpp::Transp mwo_trans
		,EMatrixDummyArg
		) const;

	//@}

public:

	/** Overridden from MatrixNonsing */
	//@{
	/// Returns <tt>this->clone_msns()</tt>.
	mat_mns_mut_ptr_t clone_mns();
	/// Returns <tt>this->clone_msns()</tt>.
	mat_mns_ptr_t clone_mns() const;
	//@}

};	// end class MatrixSymNonsing

// ////////////////////////////////////////////////////////////////////////////////////////////////
/** \defgroup MatrixSymNonsingular_func_grp MatrixSymNonsing non-member functions that call virtual functions.
  *
  * These allow nonmember functions to act like virtual functions.
  */
//@{

inline
/// sym_gms_lhs = alpha * op(mwo) * inv(mswof) * op(mwo)'
void M_StMtInvMtM(
	MatrixSymOp* sym_gms_lhs, value_type alpha
	,const MatrixOp& mwo
	,BLAS_Cpp::Transp mwo_trans, const MatrixSymNonsing& mswof
	,MatrixSymNonsing::EMatrixDummyArg mwo_rhs
	)
{
	mswof.M_StMtInvMtM(sym_gms_lhs,alpha,mwo,mwo_trans,mwo_rhs);
}

//@}

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_SYM_NONSINGULAR_H
