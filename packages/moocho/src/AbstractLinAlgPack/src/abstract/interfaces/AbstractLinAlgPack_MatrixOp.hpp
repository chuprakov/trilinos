// //////////////////////////////////////////////////////////////////////////////////
// MatrixWithOp.h
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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_H
#define ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_H

#include <iosfwd>

#include "MatrixBase.h"
#include "MatrixSpace.h"
#include "Range1D.h"

namespace AbstractLinAlgPack {

///
/** Base class for all matrices that support basic matrix operations.
  * 
  * These basic operations are:
  *
  * Level-1 BLAS
  *
  * #mwo_lhs += alpha * op(M_rhs)# (BLAS xAXPY)\\
  * #mwo_lhs += alpha * op(M_rhs)  * op(P_rhs)#\\
  * #mwo_lhs += alpha * op(P_rhs)  * op(M_rhs)#\\
  * #mwo_lhs += alpha * op(P1_rhs) * op(M_rhs) * op(P2_rhs)#\\
  *
  * #M_lhs += alpha * op(mwo_rhs)# (BLAS xAXPY)\\
  * #M_lhs += alpha * op(mwo_rhs) * op(P_rhs)#\\
  * #M_lhs += alpha * op(P_rhs)   * op(mwo_rhs)#\\
  * #M_lhs += alpha * op(P1_rhs)  * op(mwo_rhs) * op(P2_rhs)#\\
  *
  * Level-2 BLAS
  *
  * #vs_lhs = alpha * op(M_rhs1) * v_rhs2 + beta * vs_lhs# (BLAS xGEMV)\\
  * #vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs# (BLAS xGEMV)\\
  * #vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * v_rhs3 + beta * vs_lhs#\\
  * #vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * sv_rhs3 + beta * vs_lhs#\\
  * #result = v_rhs1' * op(M_rhs2) * v_rhs3#\\
  * #result = sv_rhs1' * op(M_rhs2) * sv_rhs3#\\
  *
  * Level-3 BLAS
  *
  * #mwo_lhs = alpha * op(M_rhs1)   * op(mwo_rhs2) + beta * mwo_lhs# (right) (xGEMM)\\
  * #mwo_lhs = alpha * op(mwo_rhs1) * op(M_rhs2)   + beta * mwo_lhs# (left)  (xGEMM)\\
  *
  * #M_lhs   = alpha * op(mwo_rhs1) * op(mwo_rhs2) + beta * M_lhs#           (xGEMM)\\
  *
  * All of the Level-1, Level-2 and Level-3 BLAS operations have default implementations
  * based on the Level-2 BLAS operation:\\
  *
  * #vs_lhs = alpha * op(M_rhs1) * v_rhs2 + beta * vs_lhs# (BLAS xGEMV) \\
  *
  * The only methods that have to be overridden are #space_rows()#, #space_cols()#
  * and the single #Vp_StMtV(...)# method shown above for dense vectors.
  *
  * This is to allow fast prototyping of the matrix and for postponement of writting
  * specialized methods of other time critical operations until later if they are needed.
  *
  * All of the Level-1 and Level-3 methods with #this# matrix acting as a rhs
  * argument are only implemented for the case where the lhs matrix is a mutable
  * matrix with row, column and diagonal access.  All of the default implementations
  * for the Level-1 and Level-3* methods with #this# acting as the lhs argument throw
  * exceptions.  The idea behind this set of methods and the default implementations
  * being the way they are is to try to deal with the multiple dispatch problem
  * (i.e. selecting the proper method at run time based on two or more abstract
  * arguments).  The idea is that the client calls a non-member function of the
  * for #Mp_StM(...)# or #Mp_StMtM(...)# (these are defined later).  These non-member
  * functions call the method on the first rhs abstract MatrixWithOp argument.  These
  * are the methods that have default implementations based on #MatrixWithOpMutable#.
  * Therefore, if there is a common mutable matrix associated with any particular
  * application area (i.e. a standard dense matrix class for a serial application)
  * then these default methods will do the trick based on the implementation of
  * #Vp_StMtV(...)#.  However, in more specialized applications, the lhs matrix
  * argument may not be able to support the #MatrixWithOpMutable# interface.
  * If this is the case, then the default implementations boot the problem
  * to the next matrix object in hope that it can implement the methods.
  * If this is the default implementation for #Mp_StM(MatrixWithOp*,value_type,Transp)#
  * for instance, then it calls #Mp_StM(value_type,const MatrixWithOp&,Transp)# in
  * hopes that the lhs matrix object has implemented this method.  If not, the
  * default implementation throw an excpetion.  Therefore, we have done all we
  * can do to provide a default implementation but in the end we may fail.  This
  * is the price we have to pay if we are to allow fully abstract lhs and rhs
  * matrix arguments in matrix operations.  The only alternative is to implement
  * a full blown multiple dispatch method calling infrastructure which would
  * be a burden even the simplest matrix classes.
  */
class MatrixWithOp : public virtual MatrixBase {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<const MatrixWithOp>    mat_ptr_t;

	/// Thrown if a method is not implemented
	class MethodNotImplemented : public std::runtime_error
	{public: MethodNotImplemented(const std::string& what_arg) : std::runtime_error(what_arg) {}};

	/// Thrown if matrices are not compatible
	class IncompatibleMatrices : public std::logic_error
	{public: IncompatibleMatrices(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/** @name Vector spaces for the rows and columns of the matrix.
	 *
	 * Note that the vectors space objects returned are specifically bound to this
	 * matrix object.  The vector space objects returned should only be considered
	 * to be transient and may become invalid if #this# is modified in some way
	 * (but not through the MatrixWithOp interface obviously).
	 */
	//@{

	/// Vector space for vectors that are compatible with the rows of the matrix.
	virtual const VectorSpace& space_rows() const = 0;

	/// Vector space for vectors that are compatible with the columns of the matrix.
	virtual const VectorSpace& space_cols() const = 0;

	//@}

	///
	/** Create a transient constant sub-matrix view of this matrix (if supported).
	 *
	 * This view is to be used immediatly and then discarded.
	 *
	 * This method can only be expected to return #return.get() != NULL# if
	 * #this->space_cols().sub_space(row_rng) != NULL# and
	 * #this->space_rows().sub_space(col_rng) != NULL#.  The default
	 * implementation returns #return.get() == NULL# unless #row_rng# and
	 * #col_rng# select the whole matrix. (i.e. #row_rng = [1,this->rows()]#
	 * or #row_rng.full_range() == true#.
	 */
	virtual mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
	
	///
	/** Inlined implementation calls #this->sub_view(Range1D(rl,ru),Range1D(cl,cu))#.
	 */
	virtual mat_ptr_t sub_view(
		const index_type& rl, const index_type& ru
		,const index_type& cl, const index_type& cu
		) const;

	///
	/** Zero out the matrix.
	 *
	 * The default implementation throws an exception.  This is not
	 * the best design but it meets some needs.  Any matrix implementation
	 * can implement this method and mimic the behavior.  However, only
	 * matrices that are going to be on the lhs (non-const) of a Level-1
	 * or Level-3 BLAS operation need every implement this method.
	 */
	virtual MatrixWithOp& zero_out();

	///
	/** Virtual assignment operator.
	  *
	  * The default implementation just throws a std::logic_error
	  * exception if it is not assignment to self.  A more specialized
	  * implementation could use this to copy the state to #this# matrix
	  * object from a compatible #M# matrix object.
	  */
	virtual MatrixWithOp& operator=(const MatrixWithOp& M);

	///
	/** Virtual output function.
	  *
	  * The default implementaion just extracts rows on at
	  * a time by calling #this->Vp_StMtV( ..., trans, ... )# and
	  * then prints the rows.
	  */
	virtual std::ostream& output(std::ostream& out) const;

	// /////////////////////////////////////////////////////
	/** @name Level-1 BLAS */
	//@{
	
	/** @name Implemented by rhs matrix object
	 *
	 * The default implementations of the methods require that
	 * #dynamic_cast<MatrixWithOpMutable*>(mwo_lhs) != NULL#.
	 * If not, then these default implementations call the
	 * versions of these methods on the mwo_lhs argument
	 * in hopes that it can implement the operation.
	 */
	//@{

	/// mwo_lhs += alpha * op(M_rhs) (BLAS xAXPY)
	virtual void Mp_StM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs) const;

	/// mwo_lhs += alpha * op(M_rhs) * op(P_rhs)
	virtual void Mp_StMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp M_trans
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		) const;

	/// mwo_lhs += alpha * op(P) * op(M_rhs)
	virtual void Mp_StPtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		, BLAS_Cpp::Transp M_trans
		) const;

	/// mwo_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
	virtual void Mp_StPtMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		) const;

	//@}

	/** @name Implemented by lhs matrix object
	 *
	 * The default implementation of these operations throw
	 * exceptions.
	 */
	//@{

	/// mwo_lhs += alpha * op(M_rhs) (BLAS xAXPY)
	virtual void Mp_StM(
		value_type alpha,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs);

	/// mwo_lhs += alpha * op(M_rhs) * op(P_rhs)
	virtual void Mp_StMtP(
		value_type alpha
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		);

	/// mwo_lhs += alpha * op(P) * op(M_rhs)
	virtual void Mp_StPtM(
		value_type alpha
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		);

	/// mwo_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
	virtual void Mp_StPtMtP(
		value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		);

	//@}

	//		end Level-1 BLAS
	//@}

	// ////////////////////////////////////////////////////
	/** @name Level-2 BLAS */
	//@{

	/// vs_lhs = alpha * op(M_rhs1) * v_rhs2 + beta * vs_lhs (BLAS xGEMV)
	virtual void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorWithOp& v_rhs2, value_type beta) const = 0;

	/// vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
	virtual void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;

	/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * v_rhs3 + beta * v_rhs
	virtual void Vp_StPtMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorWithOp& v_rhs3, value_type beta) const;

	/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * sv_rhs3 + beta * v_rhs
	virtual void Vp_StPtMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;

	/// result = v_rhs1' * op(M_rhs2) * v_rhs3
	virtual value_type transVtMtV(
		const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const VectorWithOp& v_rhs3) const;

	/// result = sv_rhs1' * op(M_rhs2) * sv_rhs3
	virtual value_type transVtMtV(
		const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;

	///
	/** Perform a specialized rank-2k update of a dense symmetric matrix of the form:
	 *
	 * #symwo_lhs += alpha*op(P1')*op(M)*op(P2) + alpha*op(P2')*op(M')*op(P1) + beta*symwo_lhs#
	 *
	 * The reason that this operation is being classified as a level-2 operation is that the
	 * total flops should be of O(n^2) and not O(n^2*k).
	 *
	 * The default implementation is based on #Mp_StMtP(...)# and #Mp_StPtM(...)#.  Of course
	 * in situations where this default implemention is inefficient the subclass should
	 * override this method.
	 */
	virtual void syr2k(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
		, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
		, value_type beta, MatrixSymWithOp* symwo_lhs ) const;

	// ////////////////////////////////////////////////////
	/** @name Level-3 BLAS */
	//@{

	///
	/** mwo_lhs = alpha * op(M_rhs1) * op(mwo_rhs2) + beta * mwo_lhs (left) (xGEMM)
	 *
	 * The default implementation of this method calls the method
	 * #mwo_rhs2.Mp_StMtM(...)#  in hopes that this object can implement
	 * the method.  If this method throws the exception #MethodNotImplemented#
	 * then the method #mwo_lhs.Mp_StMtM(...)#  is called.  If this method throws
	 * the exception #MethodNotImplemented# then
	 * #dynamic_cast<MatrixWithOpMutable*>(mwo_lhs)# is tried.  If this
	 * returns #!= NULL# then the operation is implemented in terms of
	 * #Vp_StMtV(...)#.  If it returns #== NULL# then the excpetion
	 * #IncompatibleMatrices# is thrown and the client is given the heads
	 * up that something is wrong.
	 */
	virtual void Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;

	///
	/** mwo_lhs = alpha * op(mwo_rhs1) * op(M_rhs2) + beta * mwo_lhs (right) (xGEMM)
	 *
	 * The default implementation just throws the exception #MethodNotImplemented#.
	 */
	virtual void Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;

	///
	/** M_lhs = alpha * op(mwo_rhs1) * op(mwo_rhs2) + beta * mwo_lhs (left) (xGEMM)
	 *
	 * The default implementation of this method throws the exception #MethodNotImplemented#.
	 * In some specialized applications, it is possible for the lhs matrix object to implment
	 * this method.
	 */
	virtual void Mp_StMtM(
		value_type alpha
		,const MatrixWithOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
		,value_type beta );
	
	///
	/** Perform a rank-k update of a dense symmetric matrix of the form:
	 *
	 * #symwo_lhs += op(M)*op(M') + beta*symwo_lhs#
	 *
	 * The default implementation is based on #Vp_StMtV(...)#.  Of course
	 * in situations where this default implemention is inefficient
	 * the subclass should override this method.
	 */
	virtual void syrk(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, value_type beta, MatrixSymWithOp* sym_lhs ) const;

	//		end Level-3 BLAS
	//@}

	/** @name Overridden from MatrixBase */
	//@{

	/// Returns #space_cols().dim()#
	size_type rows() const;

	/// Returns #space_rows().dim()#
	size_type cols() const;

	//@}

};	// end class MatrixWithOp

// ///////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////
/** @name MatrixWithOp inline non-member functions that call virtual functions.
  *
  * @memo These allow nonmember functions to act like virtual functions
  * and thereby allow the same syntax as in LinAlgPack.
  */
//@{

// /////////////////////////////////////////////////////
/** @name Level-1 BLAS */
//@{

/// mwo_lhs += alpha * op(M_rhs) (BLAS xAXPY)
inline void Mp_StM(
	MatrixWithOp* mwo_lhs, value_type alpha, const MatrixWithOp& M_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	M_rhs.Mp_StM(mwo_lhs,alpha,trans_rhs);
}

/// mwo_lhs += alpha * op(M_rhs) * op(P_rhs)
inline void Mp_StMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	)
{
	M_rhs.Mp_StMtP(mwo_lhs,alpha,M_trans,P_rhs,P_rhs_trans);
}

/// mwo_lhs += alpha * op(P) * op(M_rhs)
inline void Mp_StPtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	)
{
	M_rhs.Mp_StPtM(mwo_lhs,alpha,P_rhs,P_rhs_trans,M_trans);
}

/// mwo_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
inline void Mp_StPtMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs
	, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	)
{
	M_rhs.Mp_StPtMtP(mwo_lhs,alpha,P_rhs1,P_rhs1_trans,trans_rhs,P_rhs2,P_rhs2_trans);
}

//		end Level-1 BLAS
//@}

// ////////////////////////////////////////////////////
/** @name Level-2 BLAS */
//@{

/// vs_lhs = alpha * op(M_rhs1) * v_rhs2 + beta * vs_lhs (BLAS xGEMV)
inline void Vp_StMtV(
	VectorWithOpMutable* vs_lhs, value_type alpha, const MatrixWithOp& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorWithOp& v_rhs2, value_type beta = 1.0)
{
	M_rhs1.Vp_StMtV(vs_lhs,alpha,trans_rhs1,v_rhs2,beta);
}

/// vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
inline void Vp_StMtV(
	VectorWithOpMutable* vs_lhs, value_type alpha, const MatrixWithOp& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2, value_type beta = 1.0)
{
	M_rhs1.Vp_StMtV(vs_lhs,alpha,trans_rhs1,sv_rhs2,beta);
}

/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * v_rhs3 + beta * v_rhs
inline void Vp_StPtMtV(
	VectorWithOpMutable* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp M_rhs2_trans
	, const VectorWithOp& v_rhs3, value_type beta = 1.0) 
{
	M_rhs2.Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,v_rhs3,beta);
}

/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * sv_rhs3 + beta * v_rhs
inline void Vp_StPtMtV(
	VectorWithOpMutable* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta = 1.0)
{
	M_rhs2.Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,sv_rhs3,beta);
}

/// result = v_rhs1' * op(M_rhs2) * v_rhs3
inline value_type transVtMtV(
	const VectorWithOp& v_rhs1, const MatrixWithOp& M_rhs2
	, BLAS_Cpp::Transp trans_rhs2, const VectorWithOp& v_rhs3)
{
	return M_rhs2.transVtMtV(v_rhs1,trans_rhs2,v_rhs3);
}

/// result = sv_rhs1' * op(M_rhs2) * sv_rhs3
inline value_type transVtMtV(
	const SpVectorSlice& sv_rhs1, const MatrixWithOp& M_rhs2
	, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3)
{
	return M_rhs2.transVtMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

//		end Level-2 BLAS
//@}

// ////////////////////////////////////////////////////
/** @name Level-3 BLAS */
//@{

/// mwo_lhs = alpha * op(mwo_rhs1) * op(mwo_rhs2) + beta * mwo_lhs (right) (xGEMM)
inline void Mp_StMtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta = 1.0 )
{
	mwo_rhs1.Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta);
}

//		end Level-3 BLAS
//@}

//		end Inline non-member functions
//@}

// //////////////////////////////////////////////////
// Inlined methods for MatrixWithOp

inline
MatrixWithOp::mat_ptr_t
MatrixWithOp::sub_view(
	const index_type& rl, const index_type& ru
	,const index_type& cl, const index_type& cu
	) const
{
	return this->sub_view(Range1D(rl,ru),Range1D(cl,cu));
}

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_H
