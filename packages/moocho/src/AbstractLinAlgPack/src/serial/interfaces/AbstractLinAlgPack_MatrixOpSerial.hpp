// //////////////////////////////////////////////////////////////////////////////////
// MatrixWithOpSerial.h
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

#ifndef SPARSE_LINALG_PACK_MATRIX_WITH_OP_SERIAL_H
#define SPARSE_LINALG_PACK_MATRIX_WITH_OP_SERIAL_H

#include <iosfwd>

#include "SparseLinAlgPackTypes.h"
#include "VectorSpaceSerial.h"
#include "AbstractLinAlgPack/include/MatrixWithOp.h"

namespace SparseLinAlgPack {

///
/** Base class for all matrices implemented in a shared memory address space.
 * 
 * This base class inherites from <tt>AbstractLinAlgPack::MatrixWithOp</tt> and
 * provides a means of transition from purely abstract (i.e. parallel, out-of-core)
 * linear algebra to serial linear algebra (with explicit dense vectors and
 * matrices.
 *
 * This matrix subclass simply uses <tt>dynamic_cast<></tt> to map from the abstract interfaces
 * <tt>AbstractLinAlgPack::VectorWithOp</tt>, <tt>AbstractLinAlgPack::VectorWithOpMutable</tt>,
 * <tt>AbstractLinAlgPack::MatrixWithOp</tt> to a concreate types using a few standarized serial
 * interfaces <tt>VectorWithOpGetSparse</tt>, <tt>VectorWithOpMutableDense</tt> ...
 * If these interfaces are not supported by the vector and matrix arguments, then exceptions may
 * be thrown.
 *
 * Note that every method from \c MatrixWithOp is not overridden here.  Specifically, the methods
 * <tt>MatrixWithOp::zero_out()</tt>  and <tt>MatrixWithOp::sub_view()</tt>
 * are not overridden.  These methods have default implementations that should be sufficient
 * for most uses but subclasses may want to provide specialized implementations anyway.
 *
 * ToDo: Finish documentation!
 */
class MatrixWithOpSerial : public virtual MatrixWithOp {
public:

	/** @name Level-1 BLAS */
	//@{

	/// gms_lhs += alpha * op(M_rhs) (BLAS xAXPY)
	virtual void Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs) const;

	/// gms_lhs += alpha * op(M_rhs) * op(P_rhs)
	virtual void Mp_StMtP(GenMatrixSlice* gms_lhs, value_type alpha
		, BLAS_Cpp::Transp M_trans
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		) const;

	/// gms_lhs += alpha * op(P) * op(M_rhs)
	virtual void Mp_StPtM(GenMatrixSlice* gms_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		, BLAS_Cpp::Transp M_trans
		) const;

	/// gms_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
	virtual void Mp_StPtMtP(GenMatrixSlice* gms_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_trans
		, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		) const;

	//@}

	/** @name Level-2 BLAS */
	//@{

	/// vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
	virtual void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const = 0;

	/// vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
	virtual void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;

	/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * vs_rhs3 + beta * vs_rhs
	virtual void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorSlice& vs_rhs3, value_type beta) const;

	/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * sv_rhs3 + beta * vs_rhs
	virtual void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;

	/// result = vs_rhs1' * op(M_rhs2) * vs_rhs3
	virtual value_type transVtMtV(const VectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const VectorSlice& vs_rhs3) const;

	///
	/** Perform a specialized rank-2k update of a dense symmetric matrix of the form:
	 *
	 * #sym_lhs += alpha*op(P1')*op(M)*op(P2) + alpha*op(P2')*op(M')*op(P1) + beta*sym_lhs#
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
		, value_type beta, sym_gms* sym_lhs ) const;

	//@}

	/** @name Level-3 BLAS */
	//@{

	/// gms_lhs = alpha * op(M_rhs1) * op(gms_rhs2) + beta * gms_lhs (right) (xGEMM)
	virtual void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	/// gms_lhs = alpha * op(gms_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xGEMM)
	virtual void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
		, BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	/// gms_lhs = alpha * op(M_rhs1) * op(mwo_rhs2) + beta * gms_lhs (right) (xGEMM)
	virtual void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOpSerial& mwo_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	/// gms_lhs = alpha * op(M_rhs1) * op(sym_rhs2) + beta * gms_lhs (right) (xSYMM)
	virtual void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const sym_gms& sym_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	/// gms_lhs = alpha * op(sym_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xSYMM)
	virtual void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const sym_gms& sym_rhs1
		, BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	/// gms_lhs = alpha * op(M_rhs1) * op(tri_rhs2) + beta * gms_lhs (right) (xTRMM)
	virtual void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const tri_gms& tri_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	/// gms_lhs = alpha * op(tri_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xTRMM)
	virtual void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const tri_gms& tri_rhs1
		, BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;


	///
	/** Perform a rank-k update of a dense symmetric matrix of the form:
	 *
	 * #sym_lhs += op(M)*op(M') + beta*sym_lhs#
	 *
	 * The default implementation is based on #Vp_StMtV(...)#.  Of course
	 * in situations where this default implemention is inefficient
	 * the subclass should override this method.
	 */
	virtual void syrk(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, value_type beta, sym_gms* sym_lhs ) const;

	//@}

	/** @name Overridden from MatrixWithOp */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;
	///
	std::ostream& output(std::ostream& out) const;
	///
	bool Mp_StM(
		MatrixWithOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs
		) const;
	///
	bool Mp_StMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp M_trans
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		) const;
	///
	bool Mp_StPtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		, BLAS_Cpp::Transp M_trans
		) const;
	///
	bool Mp_StPtMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		) const;
	///
	bool Mp_StM(
		value_type alpha,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs);
	///
	bool Mp_StMtP(
		value_type alpha
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		);
	///
	bool Mp_StPtM(
		value_type alpha
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		);
	///
	bool Mp_StPtMtP(
		value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		);
	///
	void Vp_StMtV(
		VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorWithOp& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* v_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorWithOp& v_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* v_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;
	///
	value_type transVtMtV(
		const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const VectorWithOp& v_rhs3) const;
	///
	value_type transVtMtV(
		const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;
	///
	void syr2k(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
		, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
		, value_type beta, MatrixSymWithOp* symwo_lhs ) const;
	///
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		,value_type beta ) const;
	///
	void syrk(
		BLAS_Cpp::Transp M_trans, value_type alpha
		,value_type beta, MatrixSymWithOp* sym_lhs ) const;

	//@}

private:

	// ////////////////////////////////////
	// Private data members

	mutable VectorSpaceSerial       space_cols_;
	mutable VectorSpaceSerial       space_rows_;

	// ////////////////////////////////////
	// Private member functions


};	// end class MatrixWithOpSerial

/** \group MatrixWithOpSerial_funcs MatrixWithOpSerial inline non-member functions that call virtual functions.
  *
  * @memo These allow nonmember functions to act like virtual functions
  * and thereby allow the same syntax as in <tt>LinAlgPack</tt>.
  */
//@{

/** @name Level-1 BLAS */
//@{

/// gms_lhs += alpha * op(M_rhs) (BLAS xAXPY)
inline void Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha, const MatrixWithOpSerial& M_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	M_rhs.Mp_StM(gms_lhs,alpha,trans_rhs);
}

/// gms_lhs += alpha * op(M_rhs) * op(P_rhs)
inline void Mp_StMtP(GenMatrixSlice* gms_lhs, value_type alpha
	, const MatrixWithOpSerial& M_rhs, BLAS_Cpp::Transp M_trans
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	)
{
	M_rhs.Mp_StMtP(gms_lhs,alpha,M_trans,P_rhs,P_rhs_trans);
}

/// gms_lhs += alpha * op(P) * op(M_rhs)
inline void Mp_StPtM(GenMatrixSlice* gms_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
	, const MatrixWithOpSerial& M_rhs, BLAS_Cpp::Transp M_trans
	)
{
	M_rhs.Mp_StPtM(gms_lhs,alpha,P_rhs,P_rhs_trans,M_trans);
}

/// gms_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
inline void Mp_StPtMtP(GenMatrixSlice* gms_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, const MatrixWithOpSerial& M_rhs, BLAS_Cpp::Transp trans_rhs
	, const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
	)
{
	M_rhs.Mp_StPtMtP(gms_lhs,alpha,P_rhs1,P_rhs1_trans,trans_rhs,P_rhs2,P_rhs2_trans);
}

//@}

/** @name Level-2 BLAS */
//@{

/// vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
inline void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, const MatrixWithOpSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2, value_type beta = 1.0)
{
	M_rhs1.Vp_StMtV(vs_lhs,alpha,trans_rhs1,vs_rhs2,beta);
}

/// vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
inline void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, const MatrixWithOpSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2, value_type beta = 1.0)
{
	M_rhs1.Vp_StMtV(vs_lhs,alpha,trans_rhs1,sv_rhs2,beta);
}

/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * vs_rhs3 + beta * vs_rhs
inline void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, const MatrixWithOpSerial& M_rhs2, BLAS_Cpp::Transp M_rhs2_trans
	, const VectorSlice& vs_rhs3, value_type beta = 1.0) 
{
	M_rhs2.Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,vs_rhs3,beta);
}

/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * sv_rhs3 + beta * vs_rhs
inline void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
	, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
	, const MatrixWithOpSerial& M_rhs2, BLAS_Cpp::Transp M_rhs2_trans
	, const SpVectorSlice& sv_rhs3, value_type beta = 1.0)
{
	M_rhs2.Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,sv_rhs3,beta);
}

/// result = vs_rhs1' * op(M_rhs2) * vs_rhs3
inline value_type transVtMtV(const VectorSlice& vs_rhs1, const MatrixWithOpSerial& M_rhs2
	, BLAS_Cpp::Transp trans_rhs2, const VectorSlice& vs_rhs3)
{
	return M_rhs2.transVtMtV(vs_rhs1,trans_rhs2,vs_rhs3);
}

/// result = sv_rhs1' * op(M_rhs2) * sv_rhs3
inline value_type transVtMtV(const SpVectorSlice& sv_rhs1, const MatrixWithOpSerial& M_rhs2
	, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3)
{
	return M_rhs2.transVtMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

//@}

/** @name Level-3 BLAS */
//@{

/// gms_lhs = alpha * op(M_rhs1) * op(gms_rhs2) + beta * gms_lhs (right) (xGEMM)
inline void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const MatrixWithOpSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta = 1.0)
{
	M_rhs1.Mp_StMtM(gms_lhs,alpha,trans_rhs1,gms_rhs2,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(gms_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xGEMM)
inline void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOpSerial& M_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta = 1.0)
{
	M_rhs2.Mp_StMtM(gms_lhs,alpha,gms_rhs1,trans_rhs1,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(mwo_rhs1) * op(mwo_rhs2) + beta * gms_lhs (right) (xGEMM)
inline void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const MatrixWithOpSerial& mwo_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOpSerial& gms_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta = 1.0)
{
	mwo_rhs1.Mp_StMtM(gms_lhs,alpha,trans_rhs1,gms_rhs2,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(M_rhs1) * op(sym_rhs2) + beta * gms_lhs (right) (xSYMM)
inline void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const MatrixWithOpSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const sym_gms& sym_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta = 1.0)
{
	M_rhs1.Mp_StMtM(gms_lhs,alpha,trans_rhs1,sym_rhs2,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(sym_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xSYMM)
inline void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const sym_gms& sym_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOpSerial& M_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta = 1.0)
{
	M_rhs2.Mp_StMtM(gms_lhs,alpha,sym_rhs1,trans_rhs1,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(M_rhs1) * op(tri_rhs2) + beta * gms_lhs (right) (xTRMM)
inline void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const MatrixWithOpSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const tri_gms& tri_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta = 1.0)
{
	M_rhs1.Mp_StMtM(gms_lhs,alpha,trans_rhs1,tri_rhs2,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(tri_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xTRMM)
inline void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOpSerial& M_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta = 1.0)
{
	M_rhs2.Mp_StMtM(gms_lhs,alpha,tri_rhs1,trans_rhs1,trans_rhs2,beta);
}

//@}

//@}

}	// end namespace SparseLinAlgPack 

#endif	// SPARSE_LINALG_PACK_MATRIX_WITH_OP_SERIAL_H
