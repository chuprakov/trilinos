// /////////////////////////////////////
// MatrixHessianRelaxed.h

#ifndef MATRIX_HESSIAN_RELAXED_H
#define MATRIX_HESSIAN_RELAXED_H

#include "ConstrainedOptimizationPackTypes.h"
#include "SparseLinAlgPack/include/MatrixSymWithOp.h"

namespace ConstrainedOptimizationPack {

///
/** Represents a symmetric Hessian matrix with a relaxation variable
  * added.
  *
  * This class is used to represent the matrix:
  \begin{verbatim}
	     [ H       ]
	G =  [    bigM ]
  \end{verbatim}
  *
  */
class MatrixHessianRelaxed : public MatrixSymWithOp {
public:

	/// Construct to uninitialized
	MatrixHessianRelaxed();

	///
	/** Initialize.
	  *
	  * ToDo: Finish documentation!
	  *
	  */
	void initialize(
		  const MatrixSymWithOp	&H
		, value_type			bigM
		);

	// ///////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	// //////////////////////////////
	// Overridden from MatrixWithOp

	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorSlice& vs_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;
	///
	value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const ;

private:
	size_type				n_;	// size of H
	const MatrixSymWithOp	*H_;
	value_type				bigM_;

};	// end class MatrixHessianRelaxed

}	// end namespace ConstrainedOptimizationPack

#endif 	// MATRIX_HESSIAN_RELAXED_H
