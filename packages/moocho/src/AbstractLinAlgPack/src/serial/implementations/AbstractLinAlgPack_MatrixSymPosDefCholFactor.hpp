// /////////////////////////////////////////////////////////////////
// MatrixSymPosDefCholFactor.hpp
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

#ifndef MATRIX_SYM_POS_DEF_CHOL_FACTOR_H
#define MATRIX_SYM_POS_DEF_CHOL_FACTOR_H

#include "SparseLinAlgPackTypes.hpp"
#include "MatrixExtractInvCholFactor.hpp"
#include "MatrixSymAddDelUpdateable.hpp"
#include "MatrixSymWithOpNonsingularSerial.hpp"
#include "MatrixSymDenseInitialize.hpp"
#include "MatrixSymWithOpGetGMSSymMutable.hpp"
#include "AbstractLinAlgPack/src/MatrixSymSecantUpdateable.hpp"
#include "LinAlgPack/src/GenMatrixClass.hpp"
#include "LinAlgPack/src/GenMatrixAsTriSym.hpp"
#include "ref_count_ptr.hpp"
#include "ReleaseResource.hpp"

namespace SparseLinAlgPack {
///
/** A do all class for dense symmetric positive definite matrices
 * that stores the original matrix and/or its upper cholesky factor.
 *
 * This matrix class supports a boat load of interfaces.  It is ment to be
 * a do all matrix class for all dense positive definite matrices stored as
 * the upper cholesky factor and/or the lower symmetric nonfactorized
 * portion.  It is designed to meet many different needs and apply in many different
 * contexts.  Objects from this class can represent matrices of the form:
 \verbatim
     M = scale * U' * U\endverbatim
 * where <tt>U</tt> is an upper triangular cholesky factor (i.e. the diagonal of
 * <tt>U</tt> is positive).  Also, the lower part of <tt>M</tt> may also be stored 
 * and manipulated.  The reason for maintaining the unfactorized matrix
 * also is to allow its use in contexts where the cholesky factor may not
 * be needed since linear systems are never solved for.  It is also useful
 * for (extended precision) iterative refinement.
 *
 * The purpose of <tt>scale</tt> in <tt>M = scale * U' * U</tt> is in order to allow
 * matrices of this type to represent positive definite (scale > 0)
 * or negative definite (scale < 0) matrices.
 *
 * This class allows you to do a bunch of stuff with these matrix objects:
 * \begin{itemize}
 * \item BLAS operations (\Ref{MatrixSymWithOp})
 * \item Solve for linear systems (\Ref{MatrixSymFactorized})
 * \item Initialize from a dense symmetric matrix, provided that the matrix
 *   is positive definite (\Ref{MatrixSymDenseInitialize}) (default implementation)
 * \item Extract a dense inverse cholesky factor (\Ref{MatrixExtractInvCholFactor})
 * \item Perform BFGS positive definite secant updates (\Ref{MatrixSymSecantUpdateable})
 *       (default implementation)
 * \item Add rows/cols (provided new matrix is p.d. (scale > 0) or n.d. (scale < 0)) and
 *   remove rows/cols (\Ref{MatrixSymAddDelUpdateable}) and update the factors.
 *   Note that these operations will change the size of <tt>M</tt> and/or <tt>U</tt>.
 * \end{itemize}
 *
 * The <tt>tri_gms</tt> view <tt>U</tt> is really a subview of another <tt>GenMatrixSlice</tt>
 * <tt>MU_store</tt>.  The matrix <tt>U</tt> is defined as:
 *
 * <tt>U = tri(MU_store(U_l_r:U_l_r+M_size-1,U_l_c+1:U_l_c+M_size),upper,nonunit)</tt>.
 *
 * The <tt>sym_gms</tt> view <tt>M</tt> is really another subview of <tt>MU_store</tt>:
 *
 * <tt>M = sym(MU_store(M_l_r+1:M_l_r+M_size,M_l_c:M_l_c+M_size-1),lower)</tt>.
 *
 * The reason for the offsets in the definition of <tt>M</tt> and <tt>U</tt> above is to keep the
 * diagonal of <tt>MU_store</tt> open for use as workspace.
 *
 * To allow for more flexible resourse management, this class maintains a
 * <tt>ref_count_ptr<ReleaseResource></tt> object so whatever memory needs to be released
 * will be released when <tt>MU_store</tt> is no longer needed.  This class will also
 * allocate its own memory if none is set through the constructor or the
 * init_setup(...) method.
 *
 * In short, this class should be able to handle most uses for a dense symmetric
 * positive definite matrix stored as a cholesky factor.
 * Also, direct access is given to <tt>U</tt>, <tt>M</tt>, <tt>MU_store</tt> and the ability
 * to change the setup.  I don't know what more do you want!
 *
 * More operations will be overridden as they are needed by various applications.
 */
class MatrixSymPosDefCholFactor
	: virtual public SparseLinAlgPack::MatrixSymWithOpNonsingularSerial  // doxygen needs full name
	, virtual public SparseLinAlgPack::MatrixSymDenseInitialize          // ""
	, virtual public SparseLinAlgPack::MatrixSymWithOpGetGMSSymMutable   // ""
	, virtual public MatrixExtractInvCholFactor
	, virtual public MatrixSymSecantUpdateable
	, virtual public MatrixSymAddDelUpdateable
{
public:
	
	/** @name Public types */
	//@{

	///
	typedef MemMngPack::ref_count_ptr<
		MemMngPack::ReleaseResource>  release_resource_ptr_t;

	///
	/** PostMod class to use with <tt>MemMngPack::AbstractFactorStd</tt>.
	 */
	class PostMod {
	public:
		///
		PostMod(
			bool   maintain_original = true
			,bool  maintain_factor   = false
			,bool  allow_factor      = true
			)
			:maintain_original_(maintain_original)
			,maintain_factor_(maintain_factor)
			,allow_factor_(allow_factor)
		{}
		///
		void initialize(MatrixSymPosDefCholFactor* p) const
		{
			p->init_setup(
				NULL,MemMngPack::null,0 // Allocate own storage
				,maintain_original_,maintain_factor_,allow_factor_
				);
		}
	private:
		bool  maintain_original_;
		bool  maintain_factor_;
		bool  allow_factor_;
	}; // end PostMod

	//@}

    /** @name Constructors/initalizers */
	//@{

	///
	/** Initialize to uninitialized.
	 *
	 * By default if init_setup(...) is not called then this object
	 * will allocate its own storage and maintain_original == true
	 * and maintain_factor = false.
	 */
	MatrixSymPosDefCholFactor();

	///
	/** Initialize with possible storage.
	 *
	 * This constructor just calls <tt>this->init_setup(...)</tt>.
	 */
	MatrixSymPosDefCholFactor(
		GenMatrixSlice                    *MU_store
		,const release_resource_ptr_t&    release_resource_ptr = MemMngPack::null
		,size_type                        max_size             = 0
		,bool                             maintain_original    = true
		,bool                             maintain_factor      = false
		,bool                             allow_factor         = true
		,bool                             set_full_view        = true
		,value_type                       scale                = 1.0
		);

	///
	/** Initial storage setup and possibly the view.
	 *
	 * @param  MU_store [in/state]
	 *                  If <tt>MU_store != NULL</tt> then this matrix is 
	 *                  used to store the original matrix and/or its
	 *                  cholesky factor.  This matrix may or may not be
	 *                  initialized and ready to go.  The maxinum size
	 *                  for <tt>M</tt> that can be stored is <tt>max_size =</tt>
	 *                  <tt>min(MU_store.rows(),MU_store.cols())-1</tt>.  The reason
	 *                  for this is that the diagonal of <tt>MU_store</tt> is used
	 *                  as workspace. If <tt>maintain_original == true</tt> then any portion
	 *                  of the lower triangular region below the center
	 *                  diagonal is open game to be accessed.  If
	 *                  <tt>allow_factor == true</tt> then any portion of the
	 *                  center diagonal or the upper triangular region
	 *                  above the diagonal is fair game to be accessed.
	 *                  If <tt>MU_store == NULL</tt> then any current
	 *                  storage is unset and <tt>this</tt> matrix object is then
	 *                  free to allocate its own storage as needed.
	 * @param  release_resource_ptr
	 *                  [in] Only significant if <tt>MU_store != NULL</tt>.
	 *                  Points to a resource to be released when
	 *                  <tt>MU_store</tt> is no longer needed.
	 * @param  max_size [in] Only significant if<tt>MU_store == NULL</tt>.
	 *                  If <tt>max_size > 0</tt> then this is the maximum size the matrix will be sized to.
	 *                  If <tt>max_size == 0</tt> then the max size will be determined by the other initialization
	 *                  functions.
	 * @param  maintain_original
	 *                  [in] Always significant.  If true then the original
	 *                  matrix will be maintained in the lower triangle of
	 *                  <tt>this->MU_store()</tt> (see intro).
	 * @param  maintain_factor
	 *                  [in] Always significant.  If true then the cholesky
	 *                  factor will be maintained in the upper triangle of
	 *                  <tt>this->MU_store()</tt> (see intro).
	 * @param  allow_factor
	 *                  [in] Only significant if <tt>maintain_factor == false</tt>.
	 *                  If true then the factorization can be
	 *                  computed in the upper triangular portion
	 *                  of <tt>MU_store</tt>.  Otherwise it can not.
	 * @param  set_full_view
	 *                  [in] Only significant if <tt>MU_store != NULL</tt>.
	 *                  If true then <tt>this</tt> will be setup as the
	 *                  full view of <tt>MU_store</tt> and <tt>MU_store</tt> should
	 *                  already be initialized properly.
	 * @param  scale    [in] Only significant if <tt>MU_store != NULL</tt> and <tt>set_full_view == true</tt>
	 *                  (see intro).
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>MU_store!=NULL && set_full_view==true</tt>]
	 *      <tt>this->rows() == min( MU_store->rows(), MU_store->cols() ) - 1</tt>
	 * </ul>
	 */
	void init_setup(
		GenMatrixSlice                    *MU_store
		,const release_resource_ptr_t&    release_resource_ptr = MemMngPack::null
		,size_type                        max_size             = 0
		,bool                             maintain_original    = true
		,bool                             maintain_factor      = false
		,bool                             allow_factor         = true
		,bool                             set_full_view        = true
		,value_type                       scale                = 1.0
		);
	///
	/** Resize the view <tt>M</tt> and <tt>U</tt> of <tt>MU_store</tt>.
	 *
	 * Preconditions:\begin{itemize}
	 * \item [<tt>!allocates_storage()</tt>] <tt>MU_store().rows() > 0</tt>
	 * \item <tt>M_size >= 0</tt>
	 * \item [<tt>M_size > 1</tt>] <tt>scale != 0.0</tt>
	 * \item <tt>maintain_original || maintain_factor</tt>
	 * \item [<tt>maintain_original</tt>] <tt>1 <= M_l_r <= M_l_c</tt>
	 * \item [<tt>maintain_original</tt>] <tt>M_l_r + M_size <= MU_store_.rows()</tt>
	 * \item <tt>U_l_r >= U_l_c</tt>
	 * \item [<tt>U_l_r > 0</tt>] <tt>U_l_c + M_size <= MU_store_.cols()</tt>
	 * \end{itemize}
	 *
	 * @param  M_size  [in] Size of <tt>M</tt> (see intro).  If <tt>M_size == 0</tt> then the
	 *                 matrix will be set to a size of zero.
	 * @param  scale   [in] Only significant if <tt>M_size > 0</tt>.  If <tt>scale > 0</tt> then <tt>M</tt> must
	 *                 be p.d. and if <tt>scale < 0</tt> then n.d.
	 * @param  maintain_original
	 *                 [in] Always significant.
	 *                 If true then original <tt>M</tt> is maintained and also
	 *                 elements in <tt>this->M()</tt> are expected to be initialized
	 *                 already.  If false then the lower triangular portion of
	 *                 <tt>MU_store()</tt> is strictly off limits and will never be
	 *                 touched.
	 * @param  M_l_r   [in] Only significant if <tt>M_size > 0</tt>.  Lower row index for <tt>M</tt> (see intro)
	 * @param  M_l_c   [in] Only significant if <tt>M_size > 0</tt>.  Lower column index for <tt>M</tt> (see intro)
	 * @param  maintain_factor
	 *                 [in] Always significant.
	 *                 If true then the factor <tt>U</tt> is maintained and also
	 *                 the elements in <tt>this->U()</tt> are expected to be initialized
	 *                 already.  If false then the factor may still be computed
	 *                 and stored in the upper triangular part of <tt>MU_store()</tt>
	 *                 if needed by members (such as <tt>V_InvMtV(...)</tt>) but only if
	 *                 <tt>U_l_r > 0</tt>.
	 * @param  U_l_r   [in] If <tt>M_size == 0</tt> then only <tt>U_l_r == 0</tt> and <tt>U_l_r > 0</tt> is significant.
	 *                 Otherwise <tt>U_l_r</tt> is the lower row index for <tt>U</tt> (see intro).
	 *                 If <tt>U_l_r == 0</tt> then the upper triangular portion of <tt>MU_store()</tt> is strictly
	 *                 off limits and the factorization will never be computed.
	 *                 Any functions that need this factorization will throw
	 *                 exceptions in these cases.
	 * @param  U_l_c   [in] Only significant if <tt>U_l_r > 0</tt>.
	 *                  Lower column index for <tt>U</tt> (see intro).
	 */
	void set_view(
		size_t            M_size
		,value_type       scale
		,bool             maintain_original
		,size_t           M_l_r
		,size_t           M_l_c
		,bool             maintain_factor
		,size_t           U_l_r
		,size_t           U_l_c
		);
	///
	/** Set the default pivot tolerance.
	 */
	void pivot_tols( PivotTolerances pivot_tols );
	///
	PivotTolerances pivot_tols() const;

	//@}

    /**  @name Access representation */
	//@{

	///
	/** Get the current setup.
	 */
	void get_view_setup(
		size_t            *M_size
		,value_type       *scale
		,bool             *maintain_original
		,size_t           *M_l_r
		,size_t           *M_l_c
		,bool             *maintain_factor
		,size_t           *U_l_r
		,size_t           *U_l_c
		) const;
	///
	/** Return if this object owns and allocates storage.
	 */
	bool allocates_storage() const;

	///
	/** Get access to MU_store.
	 */
	GenMatrixSlice& MU_store();
	///
	const GenMatrixSlice& MU_store() const;
	///
	/** Get view of U.
	 */
	tri_gms U();
	///
	const tri_gms U() const;
	///
	/** Get view of lower part of M.
	 */
	sym_gms M();
	///
	const sym_gms M() const;

	//@}

	/** @name Overridden from MatrixBase */
	//@{

	///
	size_type rows() const;

	//@}

	/** @name Overridden from MatrixWithOp */
	//@{

	///
	void zero_out();
	///
	std::ostream& output(std::ostream& out) const;
	///
	bool Mp_StM(
		MatrixWithOp* m_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs
		) const;
	///
	bool Mp_StM(
		value_type alpha,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs
		);

	//@}

	/** @name Overridden from MatrixWithOpSerial */
	//@{

	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& vs_rhs2, value_type beta) const;
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

	//@}

	/** @name Overridden form MatrixSymWithOpSerial */
	//@{

	void Mp_StPtMtP( sym_gms* sym_lhs, value_type alpha
		, EMatRhsPlaceHolder dummy_place_holder
		, const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
		, value_type beta ) const;

	//@}

	/** @name Overridden from MatrixNonsingularSerial */
	//@{

	/// With throw exception if factorization is not allowed.
	void V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;
	/// With throw exception if factorization is not allowed.
	void V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;

	//@}

	/** @name Overridden from MatrixSymNonsingularSerial */
	//@{

	/// Will throw exception if factorization is not allowed.
	void M_StMtInvMtM(
		sym_gms* sym_gms_lhs, value_type alpha
		, const MatrixWithOpSerial& mwo, BLAS_Cpp::Transp mwo_trans, EMatrixDummyArg
		) const;

	//@}

	/** @name Overridden from MatrixSymDenseInitialize */
	//@{

	/// Will resize view of matrices and reset scale
	void initialize( const sym_gms& M );

	//@}

	/** @name Overridden from MatrixSymWithOpGetGMSSym */
	//@{

	///
	const LinAlgPack::sym_gms get_sym_gms_view() const;
	///
	void free_sym_gms_view(const LinAlgPack::sym_gms* sym_gms_view) const;

	//@}

	/** @name Overridden from MatrixSymWithOpGetGMSSymMutable */
	//@{

	///
	LinAlgPack::sym_gms get_sym_gms_view();
	///
	void commit_sym_gms_view(LinAlgPack::sym_gms* sym_gms_view);

	//@}

	/** @name Overridden from MatrixExtractInvCholFactor */
	//@{

	///
	void extract_inv_chol( tri_ele_gms* InvChol ) const;
	
	//@}

	/** @name Overridden from MatrixSymSecantUpdateble */
	//@{

	/// Will reset view and set scale
	void init_identity( const VectorSpace& space_diag, value_type alpha );
	/// Will reset view and set scale
	void init_diagonal( const VectorWithOp& diag );
	/// Must agree with current scale
	void secant_update(VectorWithOpMutable* s, VectorWithOpMutable* y, VectorWithOpMutable* Bs);

	//@}

	/** @name Overridden from MatrixSymAddDelUpdateble */
	//@{

	/// Will reset view of U and M and reset scale
	void initialize(
		value_type         alpha
		,size_type         max_size
		);
	/// Will reset view of U and M and reset scale
	void initialize(
		const sym_gms      &A
		,size_type         max_size
		,bool              force_factorization
		,Inertia           inertia
		,PivotTolerances   pivot_tols
		);
	///
	size_type max_size() const;
	/// Will be (rows(),0,0) if scale < 0 or (0,0,rows()) if scale > 0.
	Inertia inertia() const;
	/// Will set rows() == 0
	void set_uninitialized();
	/// Will throw exceptions if not p.d. (scale > 0) or n.d. (scale < 0).
	void augment_update(
		const VectorSlice  *t
		,value_type        alpha
		,bool              force_refactorization
		,EEigenValType     add_eigen_val
		,PivotTolerances   pivot_tols
		);
	/// Should always succeed unless user gives wrong value for drop_eigen_val
	void delete_update(
		size_type          jd
		,bool              force_refactorization
		,EEigenValType     drop_eigen_val
		,PivotTolerances   pivot_tols
		);

	//@}

private:
	
	// /////////////////////////////
	// Private data members

	bool                    maintain_original_;
	bool                    maintain_factor_;
	bool                    factor_is_updated_;    // only used if maintain_factor_ == false, otherwise assumed true
	bool                    allocates_storage_;    // If then then this object allocates the storage
	release_resource_ptr_t  release_resource_ptr_;
	GenMatrixSlice          MU_store_;
	size_t                  max_size_;
	size_t                  M_size_,               // M_size == 0 is flag that we are completely uninitialized
		                    M_l_r_,
		                    M_l_c_,
	                        U_l_r_,
	                        U_l_c_;
	value_type              scale_;
	bool                    is_diagonal_;
	PivotTolerances         pivot_tols_;
	Vector                  work_; // workspace.

	// /////////////////////////////
	// Private member functions

	void assert_storage() const;
	void allocate_storage(size_type max_size) const;
	void assert_initialized() const;
	void resize_and_zero_off_diagonal(size_type n, value_type scale);
	void update_factorization() const;

}; // end class MatrixSymPosDefCholFactor

// ///////////////////////////////////////////////////////
// Inline members for MatrixSymPosDefCholFactor

inline
bool MatrixSymPosDefCholFactor::allocates_storage() const
{
	return allocates_storage_;
}

inline
GenMatrixSlice& MatrixSymPosDefCholFactor::MU_store()
{
	return MU_store_;
}

inline
const GenMatrixSlice& MatrixSymPosDefCholFactor::MU_store() const
{
	return MU_store_;
}

inline
void MatrixSymPosDefCholFactor::get_view_setup(
	size_t            *M_size
	,value_type       *scale
	,bool             *maintain_original
	,size_t           *M_l_r
	,size_t           *M_l_c
	,bool             *maintain_factor
	,size_t           *U_l_r
	,size_t           *U_l_c
	) const
{
	*M_size               = M_size_;
	*scale                = scale_;
	*maintain_original    = maintain_original_;
	*M_l_r                = maintain_original_ ? M_l_r_ : 0;
	*M_l_c                = maintain_original_ ? M_l_c_ : 0;
	*maintain_factor      = maintain_factor_;
	*U_l_r                = maintain_factor_ ? U_l_r_ : 0;
	*U_l_c                = maintain_factor_ ? U_l_c_ : 0;
}

inline
tri_gms MatrixSymPosDefCholFactor::U()
{
	return LinAlgPack::nonconst_tri(
		MU_store_(U_l_r_,U_l_r_+M_size_-1,U_l_c_+1,U_l_c_+M_size_)
		, BLAS_Cpp::upper, BLAS_Cpp::nonunit
		);
}

inline
const tri_gms MatrixSymPosDefCholFactor::U() const
{
	return LinAlgPack::tri(
		MU_store_(U_l_r_,U_l_r_+M_size_-1,U_l_c_+1,U_l_c_+M_size_)
		, BLAS_Cpp::upper, BLAS_Cpp::nonunit
		);
}

inline
sym_gms MatrixSymPosDefCholFactor::M()
{
	return LinAlgPack::nonconst_sym(
		MU_store_(M_l_r_+1,M_l_r_+M_size_,M_l_c_,M_l_c_+M_size_-1)
		, BLAS_Cpp::lower
		);
}

inline
const sym_gms MatrixSymPosDefCholFactor::M() const
{
	return LinAlgPack::sym(
		MU_store_(M_l_r_+1,M_l_r_+M_size_,M_l_c_,M_l_c_+M_size_-1)
		, BLAS_Cpp::lower
		);
}

} // end namespace SparseLinAlgPack

#endif // MATRIX_SYM_POS_DEF_CHOL_FACTOR_H
