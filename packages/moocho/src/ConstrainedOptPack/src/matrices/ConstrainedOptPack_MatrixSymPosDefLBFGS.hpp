// //////////////////////////////////////////////////////////////////////////////////
// MatrixSymPosDefLBFGS.h

#ifndef MATRIX_SYM_POS_DEF_LBFGS_H
#define MATRIX_SYM_POS_DEF_LBFGS_H

#include <vector>

#include "MatrixSymSecantUpdateable.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpFactorized.h"
#include "LinAlgPack/include/GenMatrixAsTriSym.h"

namespace ConstrainedOptimizationPack {

///
/** Implementation of limited Memory BFGS matrix.
  *
  * The function set_num_updates_stored(l) must be called first
  * to set the number of most current updates stored.  The storage
  * requirements for this class are O( n*l + l*l ) which is
  * O(n*l) when n >> l which is expected.
  *
  *   
  */
class MatrixSymPosDefLBFGS
	: public MatrixSymWithOpFactorized
		, public MatrixSymSecantUpdateable
{
public:

	///
	MatrixSymPosDefLBFGS( int num_updates_stored = 10 );

	///
	/** Set the maximum number of update vectors
	  *
	  * Note: set_num_updates_stored(m) must be called before identity(n)
	  * is called.  When set_num_updates_stored(m) is called all current
	  * updates are lost and the matrix becomes uninitialized.
	  */
	void set_num_updates_stored(int);

	///
	/** Set whether automatic rescaling is used or not.
	  *
	  * This function must be called before a BFGS update is performed
	  * in order for it to take effect for that update.
	  */
	void auto_rescaling(bool);
	///
	bool auto_rescaling() const;

	// /////////////////////////////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	// /////////////////////////////////////////////////////////
	/** @name Overridden from MatrixWithOp */
	//@{

	///
	std::ostream& output(std::ostream& out) const;
	///
	MatrixWithOp& operator=(const MatrixWithOp& m);
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;

	//@}

	// ////////////////////////////////////////////////////////////
	/** @name Overridden from MatrixWithOpFactorized */
	//@{

	///
	void V_InvMtV( VectorSlice* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;

	//@}

	// ///////////////////////////////////////////////////////////
	/** @name Overridden from MatrixSymSecantUpdateable */
	//@{

	///
	void init_identity( size_type n, value_type alpha );

	///
	/** Actually this calls init_identity( (&diag)->size(), norm_inf(diag) ).
	  *
	  * This initialization is not convienent for this implementation.
	  * Besides, when we are using automatric rescaling (auto_rescaling == true)
	  * then this will really not matter much anyway.
	  */
	void init_diagonal( const VectorSlice& diag );

	///
	void secant_update(VectorSlice* s, VectorSlice* y, VectorSlice* Bs);

	//		end Overridden from MatrixSymSecantUpdateable
	//@}

private:

	// //////////////////////////////////
	// Private types

	// //////////////////////////////////
	// Private data members

	size_type	n_,		// Size of the matrix.
				m_,		// Maximum number of update vectors that can be stored.
				m_bar_,	// Current number of update vectors being stored.
						// 0 <= m_bar <= m
				k_bar_;	// Position of the most recently stored update vector in S & Y
						// 1 <= k_bar <= m_bar

	value_type	gamma_k_;// Scaling factor for Bo = (1/gamma_k) * I.
	bool		auto_rescaling_; // If true then gamma_k will be set to 1.0 always

	GenMatrix	S_,		// (n x m) Matrix of stored update vectors = [ s1, ..., sm ]
						// S(:,k_bar) is the most recently stored s update vector
				Y_,		// (n x m) Matrix of stored update vectors = [ y1, ..., ym ]
						// Y(:,k_bar) is the most recently stored y update vector
				STY_,	// (m x m) The matrix S'Y
				STSYTY_;// ((m+1) x (m+1)) The strictly upper triangular part stores the
						// upper triangular part Y'Y and the strictly lower triangular
						// part stores the lower triangular part of S'S.  The diagonal
						// can be used for work space.

	mutable bool		Q_updated_;	// True if Q has been updated for the most current update.
	mutable GenMatrix	QJ_;		// Used to store factorization of the schur complement of Q.

	mutable Vector		work_;	// workspace for performing operations.

	// //////////////////////////////////
	// Private member functions

	// Access to important matrices.

	///
	const GenMatrixSlice S() const;

	///
	const GenMatrixSlice Y() const;

	///
	const tri_gms R() const;

	/// Strictly lower triangular part of L
	const tri_gms Lb() const;

	///
	const sym_gms STS() const;

	///
	const sym_gms YTY() const;

	/// y = inv(Q) * x
	void V_invQtV( VectorSlice* y, const VectorSlice& x ) const;

	/// y += D * x
	void Vp_DtV( VectorSlice* y, const VectorSlice& x ) const;

	// Updates

	/// Update Q
	void update_Q() const;

	///
	void assert_initialized() const;

};	// end class MatrixSymPosDefLBFGS

// //////////////////////////////////////////////
// Inline member functions

inline
void MatrixSymPosDefLBFGS::auto_rescaling(bool auto_rescaling)
{
	auto_rescaling_ = auto_rescaling;
}

inline
bool MatrixSymPosDefLBFGS::auto_rescaling() const
{
	return auto_rescaling_;
}

}	// end namespace ConstrainedOptimizationPack 

#endif	// MATRIX_SYM_POS_DEF_LBFGS_H
