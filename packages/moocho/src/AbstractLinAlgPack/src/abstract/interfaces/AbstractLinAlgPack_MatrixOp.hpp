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
 * <tt>mwo_lhs += alpha * op(M_rhs)</tt> (BLAS xAXPY)<br>
 * <tt>mwo_lhs += alpha * op(M_rhs)  * op(P_rhs)</tt><br>
 * <tt>mwo_lhs += alpha * op(P_rhs)  * op(M_rhs)</tt><br>
 * <tt>mwo_lhs += alpha * op(P1_rhs) * op(M_rhs) * op(P2_rhs)</tt><br>
 *
 * <tt>M_lhs += alpha * op(mwo_rhs)</tt> (BLAS xAXPY)<br>
 * <tt>M_lhs += alpha * op(mwo_rhs) * op(P_rhs)</tt><br>
 * <tt>M_lhs += alpha * op(P_rhs)   * op(mwo_rhs)</tt><br>
 * <tt>M_lhs += alpha * op(P1_rhs)  * op(mwo_rhs) * op(P2_rhs)</tt><br>
 *
 * Level-2 BLAS
 *
 * <tt>vs_lhs = alpha * op(M_rhs1) * v_rhs2 + beta * vs_lhs</tt> (BLAS xGEMV)<br>
 * <tt>vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs</tt> (BLAS xGEMV)<br>
 * <tt>vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * v_rhs3 + beta * vs_lhs</tt><br>
 * <tt>vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * sv_rhs3 + beta * vs_lhs</tt><br>
 * <tt>result = v_rhs1' * op(M_rhs2) * v_rhs3</tt><br>
 * <tt>result = sv_rhs1' * op(M_rhs2) * sv_rhs3</tt><br>
 *
 * Level-3 BLAS
 *
 * <tt>mwo_lhs = alpha * op(M_rhs1)   * op(mwo_rhs2) + beta * mwo_lhs</tt> (right) (xGEMM)<br>
 * <tt>mwo_lhs = alpha * op(mwo_rhs1) * op(M_rhs2)   + beta * mwo_lhs</tt> (left)  (xGEMM)<br>
 * <tt>M_lhs   = alpha * op(mwo_rhs1) * op(mwo_rhs2) + beta * M_lhs</tt>           (xGEMM)<br>
 *
 * All of the Level-1, Level-2 and Level-3 BLAS operations have default implementations
 * based on the Level-2 BLAS operation:<br>
 *
 * <tt>vs_lhs = alpha * op(M_rhs1) * v_rhs2 + beta * vs_lhs</tt> (BLAS xGEMV)<br>
 *
 * The only methods that have to be overridden are \c space_cols(), \c space_rows()
 * and the single \c Vp_StMtV() method shown above.  This is to allow fast prototyping of
 * matrix subclasses and for postponement of writting specialized methods of other time
 * critical operations until later if they are needed.
 *
 * The vector space objects returned by the methods \c space_cols() and \c space_rows() 
 * are specifically bound to this matrix object.  The vector space objects returned should
 * only be considered to be transient and may become invalid if \c this is modified in some
 * significant way (but not through this <tt>%MatrixWithOp</tt> interface obviously).
 *
 * Most of the Level-1 through Level-3 BLAS methods should not be called directly by clients,
 * but instead be called through the \ref MatrixWithOp_func_grp "provided non-member functions".
 * The Level-1 and Level-3 methods of this class have a special protocal in order to deal
 * with the multiple dispatch problem.  In essence, the "Chain of responsibility" design
 * pattern is used to allow each of the participating matrix objects a chance to implement
 * an operation.  All of the Level-1 and Level-3 methods with <tt>this</tt> matrix acting
 * as the first rhs argument are implemented for the case where the lhs matrix supports
 * the multi-vector interface.  All of the other default implementations return false.
 * The idea behind this set of methods and the default implementations
 * being the way they are is to try to deal with the multiple dispatch problem
 * (i.e. selecting the proper method at run time based on two or more abstract
 * arguments).  The idea is that the client calls a non-member function of the
 * for <tt>Mp_StM(...)</tt> or <tt>Mp_StMtM(...)</tt> (these are defined later).
 * These non-member functions call the method on the first rhs abstract \c MatrixWithOp
 * argument.  The first thing that tests methods do is to call the methods <tt>Mp_StM(...)</tt>
 * or <tt>Mp_StMtM(...)</tt> implemented by the lhs matrix object (or the other rhs matrix
 * object for <tt>Mp_StMtM(...)</tt>) and give them a chance to implement the operation.
 * If the other matrix objects return false from these method calls and refuse to
 * implement the method, that these default implementations try to cast the lhs matrix to
 * <tt>MultiVectorMutable</tt>.  Therefore, if there is a common multi-vector matrix
 * class associated with any particular application area (i.e. a standard dense matrix
 * class for a serial application) then these default methods will do the trick
 * based only on the implementation of <tt>Vp_StMtV(...)</tt>.  If the
 * <tt>MultiVectorMutable</tt> interface is not supported then \c false is returned to signal
 * the the client that the method is not implemented.  Therefore, we have done all we can do to
 * provide a default implementation (no matter how inefficient it may be) while allowing the other matrix
 * objects to handle the method if possible, but in the end we may fail (and hense the
 * exception is thrown).  This is the price we have to pay if we are to allow fully
 * abstract lhs and rhs matrix arguments in matrix operations.  The only alternative is to implement
 * a full blown multiple dispatch method calling infrastructure which would
 * be a burden even the simplest matrix classes.
 *
 * ToDo: Add more detailed documentation for the default Level-1 and Level-3 BLAS
 * methods.
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

	/** @name Vector spaces for the rows and columns of the matrix */
	//@{

	/// Vector space for vectors that are compatible with the columns of the matrix.
	virtual const VectorSpace& space_cols() const = 0;

	/// Vector space for vectors that are compatible with the rows of the matrix.
	virtual const VectorSpace& space_rows() const = 0;

	//@}

	/** @name Minimal modifying methods */
	//@{

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
	  * implementation could use this to copy the state to <tt>this</tt> matrix
	  * object from a compatible <tt>M</tt> matrix object.
	  */
	virtual MatrixWithOp& operator=(const MatrixWithOp& M);

	//@}

	/** @name Output */
	//@{

	///
	/** Virtual output function.
	  *
	  * The default implementaion just extracts rows one at
	  * a time by calling <tt>this->Vp_StMtV()</tt> with 
	  * <tt>EtaVector</tt> objects and then prints the rows.
	  */
	virtual std::ostream& output(std::ostream& out) const;

	//@}

	/** @name Sub-matrix views */
	//@{

	///
	/** Create a transient constant sub-matrix view of this matrix (if supported).
	 *
	 * This view is to be used immediatly and then discarded.
	 *
	 * This method can only be expected to return <tt>return.get() != NULL</tt> if
	 * <tt>this->space_cols().sub_space(row_rng) != NULL</tt> and
	 * <tt>this->space_rows().sub_space(col_rng) != NULL</tt>.
	 *
	 * It is allows for a matrix implementation to return <tt>return.get() == NULL</tt>
	 * for any arbitrary subview.
	 *
	 * The default implementation uses the matrix subclass \c MatrixWithOpSubView
	 * and therefore, can return any arbitrary subview.  More specialized implementations
	 * may want to restrict the subview that can be created somewhat.
	 */
	virtual mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
	
	///
	/** Inlined implementation calls <tt>this->sub_view(Range1D(rl,ru),Range1D(cl,cu))</tt>.
	 */
	mat_ptr_t sub_view(
		const index_type& rl, const index_type& ru
		,const index_type& cl, const index_type& cu
		) const;

	//@}

	/** @name Level-1 BLAS */
	//@{

	/// mwo_lhs += alpha * op(M_rhs) (BLAS xAXPY)
	virtual bool Mp_StM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs) const;

	/// mwo_lhs += alpha * op(M_rhs) * op(P_rhs)
	virtual bool Mp_StMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp M_trans
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		) const;

	/// mwo_lhs += alpha * op(P) * op(M_rhs)
	virtual bool Mp_StPtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		, BLAS_Cpp::Transp M_trans
		) const;

	/// mwo_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
	virtual bool Mp_StPtMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		) const;

	/// mwo_lhs += alpha * op(M_rhs) (BLAS xAXPY)
	virtual bool Mp_StM(
		value_type alpha,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs);

	/// mwo_lhs += alpha * op(M_rhs) * op(P_rhs)
	virtual bool Mp_StMtP(
		value_type alpha
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		);

	/// mwo_lhs += alpha * op(P) * op(M_rhs)
	virtual bool Mp_StPtM(
		value_type alpha
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		);

	/// mwo_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
	virtual bool Mp_StPtMtP(
		value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		);

	//		end Level-1 BLAS
	//@}

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
	 * <tt>symwo_lhs += alpha*op(P1')*op(M)*op(P2) + alpha*op(P2')*op(M')*op(P1) + beta*symwo_lhs</tt>
	 *
	 * The reason that this operation is being classified as a level-2 operation is that the
	 * total flops should be of <tt>O(n^2)</tt> and not <tt>O(n^2*k)</tt>.
	 *
	 * The default implementation is based on <tt>Mp_StMtP(...)</tt> and <tt>Mp_StPtM(...)</tt>.
	 * Of course in situations where this default implemention is inefficient the subclass should
	 * override this method.
	 */
	virtual void syr2k(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
		, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
		, value_type beta, MatrixSymWithOp* symwo_lhs ) const;

	//@}

	/** @name Level-3 BLAS */
	//@{

	///
	/** mwo_lhs = alpha * op(M_rhs1) * op(mwo_rhs2) + beta * mwo_lhs (left) (xGEMM)
	 *
	 * The default implementation of this method first calls the method
	 * <tt>mwo_rhs2.Mp_StMtM(...)</tt> in hopes that this object can implement
	 * the method.  If this method returns false then the method <tt>mwo_lhs.Mp_StMtM(...)</tt>
	 * is called.  If this method returns false then
	 * <tt>dynamic_cast<MultiVectorMutable*>(mwo_lhs)</tt> is tried.  If this
	 * returns <tt>!= NULL</tt> then the operation is implemented in terms of
	 * <tt>this->Vp_StMtV()</tt>.  If \c dynamic_cast it returns <tt>== NULL</tt> then
	 * this default implementation returns \c false.  This method should not be
	 * called by a direct client.  Instead, the non-member function
	 * <tt>AbstractLinAlgPack::Mp_StMtM()</tt> should be used since it will throw
	 * an exceptiion.
	 */
	virtual bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;

	///
	/** mwo_lhs = alpha * op(mwo_rhs1) * op(M_rhs2) + beta * mwo_lhs (right) (xGEMM)
	 *
	 * The default implementation just returns \c false.
	 */
	virtual bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;

	///
	/** M_lhs = alpha * op(mwo_rhs1) * op(mwo_rhs2) + beta * mwo_lhs (left) (xGEMM)
	 *
	 * The default implementation of this method returns \c false.
	 * In some specialized applications, it is possible for the lhs matrix object to implment
	 * this method.
	 */
	virtual bool Mp_StMtM(
		value_type alpha
		,const MatrixWithOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
		,value_type beta );
	
	///
	/** Perform a rank-k update of a dense symmetric matrix of the form:
	 *
	 * <tt>symwo_lhs += op(M)*op(M') + beta*symwo_lhs</tt>
	 *
	 * The default implementation relays on <tt>dynamic_cast<MultiVectorSymMutable*>(sym_lhs)</tt>
	 * is based on <tt>this->Vp_StMtV()</tt>.  If \c dynamic_cast returns \c NULL, the this
	 * default implementation throws \c MethodNotDefined.  Of course in situations where this
	 * default implemention is inefficient the subclass should override this method.
	 */
	virtual void syrk(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, value_type beta, MatrixSymWithOp* sym_lhs ) const;

	//		end Level-3 BLAS
	//@}

	/** @name Overridden from MatrixBase */
	//@{

	/// Returns <tt>space_cols().dim()</tt>
	size_type rows() const;

	/// Returns <tt>space_rows().dim()</tt>
	size_type cols() const;

	//@}

};	// end class MatrixWithOp

// ////////////////////////////////////////////////////////////////////////////////////////////////
/** \defgroup MatrixWithOp_func_grp MatrixWithOp non-member functions that call virtual functions.
  *
  * These allow nonmember functions to act like virtual functions.  If any of these methods
  * on the subclasses are not implemented for a particular set of matrix arguments, then the
  * exception <tt>AbstractLinAlgPack::MatrixWithOp::MethodNotImplemented</tt> is thrown.
  * This will not happen as long as a compatible (vector spaces are compatible) lhs matrix
  * argument is passed in and <tt>dynamic_cast<MultiVectorMatrix*>(lhs) != NULL</tt>.
  */
//@{

/** @name Level-1 BLAS */
//@{

/// mwo_lhs += alpha * op(M_rhs) (BLAS xAXPY)
void Mp_StM(
	MatrixWithOp* mwo_lhs, value_type alpha, const MatrixWithOp& M_rhs
	, BLAS_Cpp::Transp trans_rhs);

/// mwo_lhs += alpha * op(M_rhs) * op(P_rhs)
void Mp_StMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	);

/// mwo_lhs += alpha * op(P) * op(M_rhs)
void Mp_StPtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
	);

/// mwo_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
void Mp_StPtMtP(
	MatrixWithOp* mwo_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs
	, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	);

//		end Level-1 BLAS
//@}

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

/** @name Level-3 BLAS */
//@{

/// mwo_lhs = alpha * op(mwo_rhs1) * op(mwo_rhs2) + beta * mwo_lhs (right) (xGEMM)
void Mp_StMtM(
	MatrixWithOp* mwo_lhs, value_type alpha
	,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta = 1.0 );

//		end Level-3 BLAS
//@}

//		end non-member functions
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
