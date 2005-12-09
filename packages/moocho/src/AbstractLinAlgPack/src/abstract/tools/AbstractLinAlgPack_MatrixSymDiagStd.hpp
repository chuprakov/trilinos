// ///////////////////////////////////////////
// AbstractLinAlgPack_MatrixSymDiagStd.hpp
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

#ifndef MATRIX_SYM_DIAGONAL_STD_H
#define MATRIX_SYM_DIAGONAL_STD_H

#include "AbstractLinAlgPack_MatrixSymInitDiag.hpp"
#include "AbstractLinAlgPack_MatrixSymDiag.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

///
/** Simple diagonal matrix class.
 *
 * ToDo: Implement clone_mswons() and deal with this->diag_ptr().count() > 1
 * by cloning vector if told to.  This allows lazy evaluation of the clone_mswons()
 * method.
 */
class MatrixSymDiagStd
	: public virtual MatrixSymInitDiag
	, public virtual MatrixSymDiag
{
public:

	///
	/** PostMod class to use with <tt>MemMngPack::AbstractFactorStd</tt>.
	 */
	class PostMod {
	public:
		PostMod(VectorSpace::space_ptr_t vectorSpace)
			: vectorSpace_(vectorSpace) {}

		void initialize(MatrixSymDiagStd* matrix) const
		    { matrix->initialize(vectorSpace_->create_member()); }
				 
	private:
		VectorSpace::space_ptr_t vectorSpace_;

	}; // end PostMod


	/** @name Constructors/initalizers */
	//@{

	/// Calls <tt>this->initialize()</tt>.
	MatrixSymDiagStd(
		const VectorSpace::vec_mut_ptr_t& diag   = Teuchos::null
		,bool                             unique = true
		);

	///
	/** Initialize given the diagonal vector (or no vector at all).
	 *
	 * @param  diag   [in] Vector to be used for the diagonal.  If <tt>diag.get() == NULL</tt>
	 *                then \c this will be uninitialized.
	 * @param  unique [in] Determines if the underlying \c diag vector is guaranteed to be
	 *                unique and not shared.
	 */
	void initialize(
		const VectorSpace::vec_mut_ptr_t& diag
		,bool                             unique = true
		);

	//@}

	/** @name Access */
	//@{

	///
	/** Give non-const access to the diagonal vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->diag_ptr().get() != NULL</tt> (throw <tt>???</tt>)
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	VectorMutable& diag();
	///
	const VectorSpace::vec_mut_ptr_t& diag_ptr() const;
	///
	bool unique() const;

	//@}

	/** @name Overridden from MatrixBase */
	//@{

	/// Returns 0 if not initalized (this->diag() == NULL).
	size_type rows() const;
	///
	size_type nz() const;

	//@}

	/** @name Overridden from MatrixOp */
	//@{

	///
	const VectorSpace& space_rows() const;
	///
	const VectorSpace& space_cols() const;
	///
	MatrixOp& operator=(const MatrixOp& mwo_rhs);
	///
	/** Add to a mutable matrix lhs.
	 *
	 * Preconditions:<ul>
	 * <li> #dynamic_cast<MultiVectorMutable*>(m_lhs) != NULL#.
	 * </ul>
	 */
	bool Mp_StM(MatrixOp* g_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const;
	///
	void Vp_StMtV(VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const Vector& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	/** Implements the symmetric rank-k update for all diagonal matrix lhs
	 *
	 * @return Returns <tt>true</tt> if <tt>dynamic_cast<MatrixSymDiagStd>(sym_lhs) != NULL</tt>.
	 * Otherwise, returns false.
	 */
	bool syrk(
		BLAS_Cpp::Transp   M_trans
		,value_type        alpha
		,value_type        beta
		,MatrixSymOp   *sym_lhs
		) const;

	//@}

	/** Overridden from MatrixOpNonsing */
	//@{

	///
	void V_InvMtV(VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const Vector& v_rhs2) const;
	///
	void V_InvMtV(VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;

	//@}

	/** @name Overridden from MatrixSymInitDiag */
	//@{

	///
	void init_identity( const VectorSpace& space_diag, value_type alpha );
	///
	void init_diagonal( const Vector& diag );

	//@}

	/** @name Overridden from MatrixSymDiag */
	//@{

	///
	const Vector& diag() const;

	//@}

private:

	VectorSpace::vec_mut_ptr_t     diag_;
	bool                           unique_;

	void copy_unique();

}; // end class MatrixSymDiagStd

// ////////////////////////////////////////
// Inline members

inline
bool MatrixSymDiagStd::unique() const
{
	return unique_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SYM_DIAGONAL_STD_H
