// //////////////////////////////////////////////////////////////////////////////////
// MatrixSymPosDefInvCholFactor.h

#ifndef MATRIX_SYM_POS_DEF_INV_CHOL_FACTOR_H
#define MATRIX_SYM_POS_DEF_INV_CHOL_FACTOR_H

#include "SymInvCholMatrixClass.h"
#include "MatrixSymSecantUpdateable.h"
#include "MatrixExtractInvCholFactor.h"
#include "SparseLinAlgPack/include/MatrixWithOpConcreteEncap.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpFactorized.h"

namespace ConstrainedOptimizationPack {

///
/** Implementation of MatrixWithOp abstract interface for SymInvCholMatrix
  */
class MatrixSymPosDefInvCholFactor
	: public virtual MatrixWithOpConcreteEncap<SymInvCholMatrix>
	, public virtual MatrixSymWithOpFactorized
	, public MatrixSymSecantUpdateable
	, public MatrixExtractInvCholFactor
{
public:

	///
	MatrixSymPosDefInvCholFactor()
	{}

	///
	MatrixSymPosDefInvCholFactor(const SymInvCholMatrix& m)
		: MatrixWithOpConcreteEncap<SymInvCholMatrix>(m)
	{}

	// /////////////////////////////////////////////////////////
	/** @name Overridden from Matrix */
	//@{

	/// 
	size_type cols() const;

	//@}

	// /////////////////////////////////////////////////////////
	/** @name Overridden from MatrixWithOp */
	//@{

	///
	MatrixWithOp& operator=(const MatrixWithOp& m);
	///
	std::ostream& output(std::ostream& out) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	value_type transVtMtV(const VectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const VectorSlice& vs_rhs3) const;
	///
	value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;

	//@}

	// ////////////////////////////////////////////////////////////
	/** @name Overridden from MatrixFactorized */
	//@{

	///
	void V_InvMtV(Vector* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;
	///
	void V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;
	///
	void V_InvMtV(Vector* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;
	///
	void V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;
	///
	value_type transVtInvMtV(const VectorSlice& vs_rhs1
		, BLAS_Cpp::Transp trans_rhs2, const VectorSlice& vs_rhs3) const;
	///
	value_type transVtInvMtV(const SpVectorSlice& sv_rhs1
		, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;

	//@}

	// ////////////////////////////////////////////////////////////
	/** @name Overridden from MatrixSymFactorized */
	//@{

	///
	void M_StMtInvMtM( sym_gms* sym_gms_lhs, value_type alpha
		, const MatrixWithOp& mwo, BLAS_Cpp::Transp mwo_trans, EMatrixDummyArg
		) const;

	//@}

	// ///////////////////////////////////////////////////////////
	/** @name Overridden from MatrixSymSecantUpdateable */
	//@{

	///
	void init_identity( size_type n, value_type alpha );
	///
	void init_diagonal( const VectorSlice& diag );
	///
	void secant_update(VectorSlice* s, VectorSlice* y, VectorSlice* Bs);

	//@}

	// ///////////////////////////////////////////////////////////
	/** @name Overridden from MatrixExtractInvCholFactor */
	//@{

	///
	void extract_inv_chol( tri_ele_gms* InvChol ) const;

	//@}

};	// end class MatrixSymPosDefInvCholFactor

}	// end namespace ConstrainedOptimizationPack 

#endif	// MATRIX_SYM_POS_DEF_INV_CHOL_FACTOR_H
