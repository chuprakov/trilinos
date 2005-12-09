// /////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymDiagSparse.hpp
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
#ifndef SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H
#define SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H

#include "AbstractLinAlgPack_MatrixSymOpSerial.hpp"
#include "AbstractLinAlgPack_MatrixConvertToSparse.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract base class for all serial symmetric diagonal matrices with
 * significant zeros along the diagonal.
 */
class MatrixSymDiagSparse
	: virtual public MatrixSymOpSerial
	, virtual public MatrixConvertToSparse
{
public:

	///
	/** <<std member comp>> members for how many updates to compute
	  * at once in the operation M_MtMtM(....).
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_updates_at_once )

	///
	/** The default value of num_updates_at_once == 0 is set to allow
	  * this class to determine the appropriate size internally.
	  */
	MatrixSymDiagSparse();

	/** @name To be overridden by subclass */
	//@{

	/// Give access to the sparse diagonal
	virtual const SpVectorSlice diag() const = 0;

	//@}

	/** @name Overridden from MatrixBase */
	//@{

	///
	size_type rows() const;

	//@}

	/** @name Overridden from MatrixOp */
	//@{

	///
	std::ostream& output(std::ostream& out) const;

	//@}

	/** @name Overridden from MatrixOpSerial */
	//@{

	///
	void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const DVectorSlice& vs_rhs2, value_type beta) const;

	//@}

	/** @name Overridden from MatrixSymOpSerial */
	//@{

	///
	/** Computes the dense symmetric matrix B += a*op(A')*M*op(A).
	 *
	 * This matrix is computed using a set of rank-1 updates.
	 *
	 * Runtime ~ O( (m^2)*nz )
	 *
	 * Storage ~ O( num_updates_at_once * m )
	 *
	 * Where:<ul>
	 * <li> <tt>n = A.rows() == this->rows()</tt>
	 * <li> <tt>m = A.cols()</tt>
	 * <li> <tt>nz = this->diag().nz()</tt>
	 * </ul>
	 *
	 * Note that a necessary condition for \c B to be full rank is for
	 * <tt>nz >= m</tt>.
	 *
	 * Also note that this default implementation is only for nonnegative
	 * diagonal entries.
	 */
	void Mp_StMtMtM( DMatrixSliceSym* sym_lhs, value_type alpha
		, EMatRhsPlaceHolder dummy_place_holder
		, const MatrixOpSerial& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
		, value_type beta ) const;

	//@}

	/** @name Overridden from MatrixConvertToSparse */
	//@{

	///
	index_type num_nonzeros(
		EExtractRegion        extract_region
		,EElementUniqueness   element_uniqueness
		) const;
	///
	void coor_extract_nonzeros(
		EExtractRegion                extract_region
		,EElementUniqueness           element_uniqueness
		,const index_type             len_Aval
		,value_type                   Aval[]
		,const index_type             len_Aij
		,index_type                   Arow[]
		,index_type                   Acol[]
		,const index_type             row_offset
		,const index_type             col_offset
		) const;

	//@}

};	// end class MatrixSymDiagSparse

}	// end namespace AbstractLinAlgPack

#endif	// SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_H
