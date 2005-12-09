// //////////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_MatrixSymPosDefLBFGS.hpp
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

#ifndef MATRIX_SYM_POS_DEF_LBFGS_H
#define MATRIX_SYM_POS_DEF_LBFGS_H

#include <vector>

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixSymSecant.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

///
/** Implementation of limited Memory BFGS matrix for arbitrary vector spaces.
 *
 * The function <tt>set_num_updates_stored()</tt> must be called first to set the maximum number of
 * the most recent updates that can be stored.  The storage requirements for this class are
 * <tt>O( n*m + m*m )</tt> which is <tt>O(n*m)</tt> when <tt>n >> m</tt> which is expected
 * (where \c n is the dimension of the vector space and \c m is the maximum number of updates
 * stored).
 *
 * This implementation is based on:
 *
 * Byrd, Nocedal, and Schnabel, "Representations of quasi-Newton matrices
 * and their use in limited memory methods", Mathematical Programming, 63 (1994)
 * 
 * Consider BFGS updates of the form:
 \verbatim

 ( B^{k-1}, s^{k-1}, y^{k-1} ) -> B^{k}
 
 where:
 
 B^{k} = B^{k-1} - ( (B*s)*(B*s)' / (s'*B*s) )^{k-1} + ( (y*y') / (s'*y) )^{k-1}
 
 B <: R^(n x n)
 s <: R^(n)
 y <: R^(n)
 \endverbatim
 * Now let us consider limited memory updating.  For this implementation we set:
 \verbatim

 Bo = ( 1 / gamma_k ) * I
 
 where:
               / (s^{k-1}'*y^{k-1})/(y^{k-1}'*y^{k-1})           :  if auto_rescaling() == true
     gamma_k = |
	           \  1/alpha from last call to init_identity(n,alpha) :  otherwise
 
 \endverbatim
 * Now let us define the matrices \c S and \c Y that store the update vectors
 * <tt>s^{i}</tt> and <tt>y^{i}</tt> for <tt>i = 1 ... m_bar</tt>:
 \verbatim

 S = [ s^{1}, s^{2},...,s^{m_bar} ] <: R^(n x m)
 Y = [ y^{1}, y^{2},...,y^{m_bar} ] <: R^(n x m)
  \endverbatim
 * Here we are only storing the <tt>m_bar <= m</tt> most recent update vectors
 * and their ordering is significant.  The columns \c S(:,m_bar) and \c Y(:,m_bar)
 * contain the most recent update vectors.  This is all client needs to know
 * in order to reconstruct the updates themselves.
 *
 * This class allows matrix-vector products <tt>x = B*y</tt> and the inverse
 * matrix-vector products <tt>x = inv(B)*y</tt> to be performed at a cost of
 * about <tt>O(n*m_bar^2)</tt>.
 */
class MatrixSymPosDefLBFGS
	: public AbstractLinAlgPack::MatrixSymOpNonsing // doxygen needs full path
	, public MatrixSymSecant
{
public:

	/** @name Public types */
	//@{

	///
	typedef Teuchos::RefCountPtr<const MultiVector>   multi_vec_ptr_t;

	///
	/** PostMod class to use with <tt>MemMngPack::AbstractFactorStd</tt>.
	 */
	class PostMod {
	public:
		///
		PostMod(
			size_type     m                   = 10
			,bool         maintain_original   = true
			,bool         maintain_inverse    = true
			,bool         auto_rescaling      = false
			)
			:m_(m)
			,maintain_original_(maintain_original)
			,maintain_inverse_(maintain_inverse)
			,auto_rescaling_(auto_rescaling)
		{}
		///
		void initialize(MatrixSymPosDefLBFGS* p) const
		{
			p->initial_setup(m_,maintain_original_,maintain_inverse_,auto_rescaling_);
		}
	private:
		size_type   m_;
		bool        maintain_original_;
		bool        maintain_inverse_;
		bool        auto_rescaling_;
	}; // end PostMod

	//@}

	/** @name Constructors and initializers */
	//@{

	/// Calls <tt>this->initial_setup()</tt>
	MatrixSymPosDefLBFGS(
	    size_type     m                   = 10
		,bool         maintain_original   = true
		,bool         maintain_inverse    = true
		,bool         auto_rescaling      = false
		);

	///
	/** Set whether automatic rescaling is used or not.
	  *
	  * This function must be called before a BFGS update is performed
	  * in order for it to take effect for that update.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, auto_rescaling )

	///
	/** Initial setup for the matrix.
	  *
	  * This function must be called before init_identity(n)
	  * is called in order for these setting to have affect.
	  * When this function is called all current updates are
	  * lost and the matrix becomes uninitialized.
	  *
	  * @param  m  [in] Max number of recent update vectors \a s
	  *            and \c y stored.
	  * @param  maintain_original
	  *            [in] If \c true then quantities needed to compute
	  *            <tt>x = Bk*y</tt> will be maintained, otherwise they
	  *            will not be unless needed.  This is to save
	  *            computational costs in case matrix-vector
	  *            products will never be needed.  However,
	  *            if a matrix vector product is needed then
	  *            these quantities will be computed on the fly
	  *            in order to satisfy the request.
	  * @param  maintain_inverse
	  *            [in] If \c true then quantities needed to compute
	  *            <tt>x = inv(Bk)*y = x = Hk*y</tt> will be maintained,
	  *            otherwise they will not be unless needed.
	  *            This is to save computational costs in case
	  *            inverse matrix-vector products will never be needed.
	  *            However, if the inverse product is ever needed
	  *            then the needed quantities will be computed
	  *            on the fly in order to satisfiy the request.
	  *            Because it takes so little extra work to maintain
	  *            the quantities needed for \c Hk it is recommended
	  *            to always set this to true.
	  * @param  auto_rescaling
	  *            [in] See intro.
	  */
	 void initial_setup(
	    size_type     m                   = 10
		,bool         maintain_original   = true
		,bool         maintain_inverse    = true
		,bool         auto_rescaling      = false
		 );

	// //////////////////////////////////
	// Representation access

	///
	size_type m() const;
	///
	size_type m_bar() const;
	///
	value_type gamma_k() const;
	///
	const multi_vec_ptr_t S() const;
	///
	const multi_vec_ptr_t Y() const;
	///
	bool maintain_original() const;
	///
	bool maintain_inverse() const;
	///
	/** Returns the total number of successful secant updates performed
	 * since <tt>this->init_identity()</tt> was called.
	 */
	size_type num_secant_updates() const;

	//@}

	/** @name Overridden from MatrixOp */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	std::ostream& output(std::ostream& out) const;
	///
	MatrixOp& operator=(const MatrixOp& mwo);
	///
	void Vp_StMtV(
		VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const Vector& v_rhs2, value_type beta ) const;

	//@}

	/** @name Overridden from MatrixOpNonsing */
	//@{

	///
	void V_InvMtV(
		VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const Vector& v_rhs2 ) const;

	//@}

	/** @name Overridden from MatrixSymSecant */
	//@{

	///
	void init_identity( const VectorSpace& space_diag, value_type alpha );
	///
	/** Actually this calls init_identity( diag.space(), diag.norm_inf() ).
	  *
	  * This initialization is not convienent for this implementation.
	  * Besides, when we are using automatric rescaling (auto_rescaling == true)
	  * then this will really not matter much anyway.
	  */
	void init_diagonal( const Vector& diag );
	///
	void secant_update(
		VectorMutable     *s
		,VectorMutable    *y
		,VectorMutable    *Bs
		);

	//@}

private:

	// //////////////////////////////////
	// Private types

	typedef VectorSpace::multi_vec_mut_ptr_t   multi_vec_mut_ptr_t;

	// //////////////////////////////////
	// Private data members

	bool        maintain_original_;  // If true, qualities needed for Bk will be maintained
	bool        original_is_updated_;// If true, qualities needed for Bk are already updated
	bool        maintain_inverse_;   // If true, quantities needed for Hk will be maintained
	bool        inverse_is_updated_; // If true, quantities needed for Hk are already updated

	VectorSpace::space_ptr_t
		        vec_spc_; // The vector space that everything is based on.

	size_type   n_,		// Size of the matrix.  If 0 then is uninitialized
	            m_,		// Maximum number of update vectors that can be stored.
                m_bar_,	// Current number of update vectors being stored.
		                // 0 <= m_bar <= m
                num_secant_updates_; // Records the number of secant updates performed
	value_type  gamma_k_;// Scaling factor for Bo = (1/gamma_k) * I.

	multi_vec_mut_ptr_t
                S_,     // (n x m) Matrix of stored update vectors = [ s1, ..., sm ]
		                // S(:,m_bar) is the most recently stored s update vector
	            Y_;     // (n_max x m) Matrix of stored update vectors = [ y1, ..., ym ]
						// Y(:,k_bar) is the most recently stored y update vector
	DMatrix   STY_,   // (m x m) The matrix S'Y
	            STSYTY_;// ((m+1) x (m+1)) The strictly upper triangular part stores the
	                    // upper triangular part Y'Y and the strictly lower triangular
	                    // part stores the lower triangular part of S'S.  The diagonal
	                    // can be used for workspace.

	mutable bool        Q_updated_; // True if Q has been updated for the most current update.
	mutable DMatrix   QJ_;        // Used to store factorization of the schur complement of Q.

	// //////////////////////////////////
	// Private member functions

	// Access to important matrices.

	///
	const DMatrixSliceTri R() const;
	/// Strictly lower triangular part of L
	const DMatrixSliceTri Lb() const;
	///
	DMatrixSlice STY();
	///
	const DMatrixSlice STY() const;
	///
	DMatrixSliceSym STS();
	///
	const DMatrixSliceSym STS() const;
	///
	DMatrixSliceSym YTY();
	///
	const DMatrixSliceSym YTY() const;
	/// y = inv(Q) * x
	void V_invQtV( DVectorSlice* y, const DVectorSlice& x ) const;
	/// y += D * x
	void Vp_DtV( DVectorSlice* y, const DVectorSlice& x ) const;

	// Updates

	/// Update Q
	void update_Q() const;

	///
	void assert_initialized() const;

};	// end class MatrixSymPosDefLBFGS

// //////////////////////////////////////////////
// Inline member functions

inline
size_type MatrixSymPosDefLBFGS::m() const
{
	return m_;
}

inline
size_type MatrixSymPosDefLBFGS::m_bar() const
{
	return m_bar_;
}

inline
value_type MatrixSymPosDefLBFGS::gamma_k() const
{
	return gamma_k_;
}

inline
const MatrixSymPosDefLBFGS::multi_vec_ptr_t
MatrixSymPosDefLBFGS::S() const
{
	return S_->mv_sub_view(1,n_,1,m_bar_);
}

inline
const MatrixSymPosDefLBFGS::multi_vec_ptr_t
MatrixSymPosDefLBFGS::Y() const
{
	return Y_->mv_sub_view(1,n_,1,m_bar_);
}

inline
bool MatrixSymPosDefLBFGS::maintain_original() const
{
	return maintain_original_;
}

inline
bool MatrixSymPosDefLBFGS::maintain_inverse() const
{
	return maintain_inverse_;
}

inline
size_type MatrixSymPosDefLBFGS::num_secant_updates() const
{
	return num_secant_updates_;
}

}	// end namespace ConstrainedOptPack 

#endif	// MATRIX_SYM_POS_DEF_LBFGS_H
