// ////////////////////////////////////////////////////////////////////
// MatrixNonsingular.h
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
//
#ifndef ABSTRACT_LINALG_PACK_MATRIX_NONSINGULAR_H
#define ABSTRACT_LINALG_PACK_MATRIX_NONSINGULAR_H

#include "MatrixBase.h"
#include "ref_count_ptr.h"

namespace AbstractLinAlgPack {

///
/** Abstract base class for all nonsingular polymorphic matrices that can solve
 * for linear system with but it may not be convienent to compute matrix vector
 * products {abstract}.
 * 
 * The operations supported are:
 *
 * Level-2 BLAS
 *
 * <tt>v_lhs	= inv(op(M_rhs1)) * vs_rhs2</tt><br>
 * <tt>v_lhs	= inv(op(M_rhs1)) * sv_rhs2</tt><br>
 * <tt>result	= v_rhs1' * inv(op(M_rhs2)) * v_rhs3</tt><br>
 * <tt>result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3</tt><br>
 * 
 * Level-3 BLAS
 *
 * <tt>m_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)</tt><br>
 * <tt>m_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)</tt><br>
 *
 * For the solve operations, the lhs and rhs arguments may not be the same
 * in general so don't assume that you can alias the lhs with the rhs and
 * get correct results.
 *
 * Any nonsingular matrix abstraction that can be used to solve for nonlinear
 * systems should also be able to support the \c MatrixWithOp interface.
 * Therefore, this interface is more of an implementation artifact than
 * an a legitimate domain abstraction.  However, some linear solvers that
 * can implement this interface, can not easily implement the <tt>%MatrixWithOp</tt>
 * interface and therefore this interface is justified.  A general client should never
 * use this interface directly.  Instead, the combined interface \c MatrixWithOpNonsingular
 * should be used with fully formed matrix abstractions.
 *
 * All these Level-2 and Level-3 BLAS operations have default implementations
 * based on the Level-2 BLAS operations:
 *
 * <tt>v_lhs = inv(op(M_rhs1)) * vs_rhs2</tt><br>
 *
 * which allows for fast prototyping of new matrix subclasses.
 *
 * The member functions should not be called directly but instead through
 * the \ref MatrixNonsingular_funcs_grp "provided non-member functions".
 *
 * The multiple dispatch approach taken in <tt>MatrixWithOp</tt> is not taken
 * in this interface.  This is because it is considered here that the
 * nonsingular matrix takes procedence of a general matrix arguemnt and
 * we can not expect a general matrix to know how to solve for a linear
 * system with some other nonsigular matrix.
 */
class MatrixNonsingular : public virtual MatrixBase {
public:

	/** @name Public types */
	//@{

#ifndef DOXYGEN_COMPILE
	///
	typedef MemMngPack::ref_count_ptr<const MatrixNonsingular>    mat_mns_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<MatrixNonsingular>          mat_mns_mut_ptr_t;
#endif

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
	virtual mat_mns_mut_ptr_t clone_mns();

	///
	/** Clone the const matrix object (if supported).
	 *
	 * The behavior of this method is the same as for the non-const version
	 * above except it returns a smart pointer to a const matrix object.
	 */
	virtual mat_mns_ptr_t clone_mns() const;

	//@}

	/** @name Level-2 BLAS */
	//@{

	/// v_lhs	= inv(op(M_rhs1)) * vs_rhs2
	virtual void V_InvMtV(
		VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const VectorWithOp& v_rhs2) const = 0;
	/// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
	virtual void V_InvMtV(
		VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;
	/// result	= vs_rhs1' * inv(op(M_rhs2)) * vs_rhs3
	virtual value_type transVtInvMtV(
		const VectorWithOp& v_rhs1
		,BLAS_Cpp::Transp trans_rhs2, const VectorWithOp& v_rhs3) const;
	/// result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3
	virtual value_type transVtInvMtV(
		const SpVectorSlice& sv_rhs1
		,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;

	//		end Level-2 BLAS
	//@}

	/** @name Level-3 BLAS */
	//@{

	///
	/** m_lhs = alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right).
	 *
	 * The default implemention performs a <tt>dynamic_cast<MultiVectorSymMutable>(m_lhs)</tt>.
	 * If this \c dynamic_cast<> returns true, then this operation is implemented in terms of
	 * <tt>this->V_InvMtV()</tt> one row or column at a time.  If this \c dynamic_cast<> returns
	 * false, then this default implementation has no choice but to throw an exception.
	 */
	virtual void M_StInvMtM(
		MatrixWithOp* m_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		) const;
	///
	/** m_lhs = alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left).
	 *
	 * The default implemention performs a <tt>dynamic_cast<MultiVectorSymMutable>(m_lhs)</tt>.
	 * If this \c dynamic_cast<> returns true, then this operation is implemented in terms of
	 * <tt>this->V_InvMtV()</tt> one row or column at a time.  If this \c dynamic_cast<> returns
	 * false, then this default implementation has no choice but to throw an exception.
	 */
	virtual void M_StMtInvM(
		MatrixWithOp* m_lhs, value_type alpha
		,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2
		) const;

	//		end Level-3 BLAS
	//@}

};	// end class MatrixNonsingular

/** \defgroup MatrixNonsingular_funcs_grp MatrixNonsingular inline non-member operation functions
  * that call virtual functions.
  *
  * These allow nonmember functions to act like virtual functions
  * and thereby allow the same syntax as in LinAlgPack.
  */
//@{

/** @name Level-2 BLAS */
//@ {

/// v_lhs	= inv(op(M_rhs1)) * v_rhs2
inline void V_InvMtV(
	VectorWithOpMutable* v_lhs, const MatrixNonsingular& M_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const VectorWithOp& v_rhs2)
{
	M_rhs1.V_InvMtV(v_lhs,trans_rhs1,v_rhs2);
}

/// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
inline void V_InvMtV(
	VectorWithOpMutable* v_lhs, const MatrixNonsingular& M_rhs1
	,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2)
{
	M_rhs1.V_InvMtV(v_lhs,trans_rhs1,sv_rhs2);
}

/// result	= v_rhs1' * inv(op(M_rhs2)) * v_rhs3
inline value_type transVtInvMtV(
	const VectorWithOp& v_rhs1, const MatrixNonsingular& M_rhs2
	,BLAS_Cpp::Transp trans_rhs2, const VectorWithOp& v_rhs3)
{
	return M_rhs2.transVtInvMtV(v_rhs1,trans_rhs2,v_rhs3);
}

/// result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3
inline value_type transVtInvMtV(
	const SpVectorSlice& sv_rhs1, const MatrixNonsingular& M_rhs2
	,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3)
{
	return M_rhs2.transVtInvMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

//		end Level-2 BLAS
//@ }

/** @name Level-3 BLAS */
//@ {

/// m_lhs	= alpha * inv(op(mwo_rhs1)) * op(mwo_rhs2) (right)
inline void M_StInvMtM(
	MatrixWithOp* m_lhs, value_type alpha
	,const MatrixNonsingular&  M_rhs1,     BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp&       mwo_rhs2,   BLAS_Cpp::Transp trans_rhs2 )
{
	M_rhs1.M_StInvMtM(m_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2);
}

/// m_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
	MatrixWithOp* m_lhs, value_type alpha
	,const MatrixWithOp&      mwo_rhs1,  BLAS_Cpp::Transp trans_rhs1
	,const MatrixNonsingular& M_rhs2,    BLAS_Cpp::Transp trans_rhs2 )
{
	M_rhs2.M_StMtInvM(m_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2);
}

//		end Level-3 BLAS
//@ }

//		end Inline non-member operation functions
//@}

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_NONSINGULAR_H
