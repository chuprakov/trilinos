// ////////////////////////////////////////////////////////////////////////////////////////////////////
// MatrixIdentConcat.h

#ifndef MATRIX_IDENT_CONCAT_H
#define MATRIX_IDENT_CONCAT_H

#include "ConstrainedOptimizationPackTypes.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"

namespace ConstrainedOptimizationPack {

///
/** Matrix class for a matrix vertically concatonated with an identity matrix {abstract}.
 *
 * Represents an interface for a matrix that represents:
 * #M = [ alpha*op(D); I ] (D_rng = [1,rows(op(D))], I_rng = [rows(op(D))+1,rows(op(D))+cols(op(D))])#
 * or
 * #M = [ I; alpha*op(D) ] (D_rng = [cols(op(D))+1,rows(op(D))+cols(op(D))], I_rng = [1,cols(op(D)))#
 * where #I# is #op(D).cols() x op(D).cols()# and where M is #op(D).rows() + op(D).cols() x op(D).cols()#.
 */
class MatrixIdentConcat : public MatrixWithOp {
public:

	/** @name Access to representation.
	 */
	//@{
	///
	virtual Range1D D_rng() const = 0;
	///
	virtual Range1D I_rng() const = 0;
	///
	virtual value_type alpha() const = 0;
	///
	virtual const MatrixWithOp& D() const = 0;
	///
	virtual BLAS_Cpp::Transp D_trans() const = 0;
	//@}

	/** @name Overridden from Matrix */
	//@{
	///
	size_type rows() const;
	///
	size_type cols() const;
	///
	size_type nz() const;
	//@}

	/** @name Overridden from MatrixWithOp */
	//@{
	///
	std::ostream& output(std::ostream& out) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	//@}

}; // end class MatrixIdentConcat

} // end namespace ConstrainedOptimizationPack

#endif // MATRIX_IDENT_CONCAT_H
