// //////////////////////////////////////////////////////////////////////////////////
// COOMatrixPartitionViewSubclass.h

#ifndef COO_MATRIX_PARTITION_VIEW_SUBCLASS_H
#define COO_MATRIX_PARTITION_VIEW_SUBCLASS_H

#include "MatrixWithOp.h"
#include "COOMatrixWithPartitionedView.h"

namespace SparseLinAlgPack {

// Could not derive this class form MatrixWithOpConcreteEncap because the assignment
// operator is not defined for the partition class.

///
/** Implementation of MatrixWithOp abstract interface for
  * COOMatrixWithPartitionedView::partition_type.
  */
class COOMatrixPartitionViewSubclass : public MatrixWithOp
{
public:

	///
	typedef COOMatrixWithPartitionedView::partition_type	M;

	///
	COOMatrixPartitionViewSubclass()
		: trans_(BLAS_Cpp::no_trans)
	{}

	///
	COOMatrixPartitionViewSubclass(BLAS_Cpp::Transp trans)
		: trans_(trans)
	{}

	///
	COOMatrixPartitionViewSubclass(const M& m)
		: m_(m), trans_(BLAS_Cpp::no_trans)
	{}

	///
	COOMatrixPartitionViewSubclass(const M& m, BLAS_Cpp::Transp trans)
		: m_(m), trans_(trans)
	{}

	///
	void set_trans(BLAS_Cpp::Transp trans) {
		trans_ = trans;
	}

	// /////////////////////////////////////////////////////
	/** @name Representation access */
	//@{

	/// Get the underlying M object
	M& m() {
		return m_;
	}

	///
	const M& m() const {
		return m_;
	}

	//		end Representation access
	//@}

	// /////////////////////////////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	///
	size_type cols() const;

	// /////////////////////////////////////////////////////
	// Overridden from MatrixWithOp

	///
	MatrixWithOp& operator=(const MatrixWithOp& m);

	// /////////////////////////////////////////////////////
	/** @name Level-1 BLAS */
	//@{

	/// (1) gms_lhs += alpha * op(M_rhs) (BLAS xAXPY)
	void Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const;

	//		end Level-1 BLAS
	//@}

	// ////////////////////////////////////////////////////
	/** @name Level-2 BLAS */
	//@{

	/// (2) vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;

	/// (3) vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;

	/// (4) result = vs_rhs1' * op(M_rhs2) * vs_rhs3
	value_type transVtMtV(const VectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const VectorSlice& vs_rhs3) const;

	/// (5) result = sv_rhs1' * op(M_rhs2) * sv_rhs3
	value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;

	//		end Level-2 BLAS
	//@}

	// ////////////////////////////////////////////////////
	/** @name Level-3 BLAS */
	//@{

	/// (6) gms_lhs = alpha * op(M_rhs1) * op(gms_rhs2) + beta * gms_lhs (right) (xGEMM)
	void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	/// (7) gms_lhs = alpha * op(gms_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xGEMM)
	void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
		, BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	//		end Level-3 BLAS
	//@}

private:
	M m_;
	BLAS_Cpp::Transp trans_;	// for how the matrix if viewed as.

	BLAS_Cpp::Transp op(BLAS_Cpp::Transp trans) const {
		using BLAS_Cpp::trans_not;
		return trans_ == BLAS_Cpp::no_trans ? trans : trans_not(trans);
	}

};	// end class COOMatrixPartitionViewSubclass

}	// end namespace SparseLinAlgPack 

#endif	// COO_MATRIX_PARTITION_VIEW_SUBCLASS_H
