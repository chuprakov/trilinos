// ////////////////////////////////////////////////////
// MatrixSymAddDelBunchKaufman.h

#ifndef MATRIX_SYM_POS_DEF_BUNCH_KAUFMAN_H
#define MATRIX_SYM_POS_DEF_BUNCH_KAUFMAN_H

#include <vector>

#include "MatrixSymAddDelUpdateableWithOpFactorized.h"
#include "MatrixSymAddDelUpdateable.h"
#include "MatrixSymPosDefCholFactor.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpFactorized.h"
#include "LinAlgPack/include/GenMatrixAsTriSym.h"

namespace ConstrainedOptimizationPack {

///
/** This class maintains the factorization of symmetric indefinite matrix
 * using a Bunch & Kaufman factorization.
 *
 * When the matix in question is positive definite or negative definite
 * then a cholesky factorization will be used for greater efficiency.
 * As much as possible the class trys to keep runtime costs down to
 * a minimum while still maintaining a stable factorization.
 *
 * ToDo: This matrix class could also handle the case of a single negative
 * or positive eigen value in an efficient manner as well.
 */
class MatrixSymAddDelBunchKaufman
	:public virtual MatrixSymWithOpFactorized
	,public virtual MatrixSymAddDelUpdateable
	,public virtual MatrixSymAddDelUpdateableWithOpFactorized
{
public:

	/// Initializes with 0x0 and pivot_tols == (0.0,0.0,0.0).
	MatrixSymAddDelBunchKaufman();

	/// Pivot tolerance used durring the cholesky factorization (it may be zero).
	void pivot_tols( PivotTolerances pivot_tols );
	///
	PivotTolerances	pivot_tols() const;

	// /////////////////////////////////////////////////////////////
	// Overridden from MatrixSymAddDelUpdateableWithOpFactorized

	///
	const MatrixSymWithOpFactorized& op_interface() const;
	///
	MatrixSymAddDelUpdateable& update_interface();
	///
	const MatrixSymAddDelUpdateable& update_interface() const;

	// /////////////////////////////////////////////////////////////
	// Overridden from MatrixSymAddDelUpdateable

	///
	void initialize(
		value_type         alpha
		,size_type         max_size
		);
	///
	void initialize(
		const sym_gms      &A
		,size_type         max_size
		,bool              force_factorization
		,Inertia           inertia
		,PivotTolerances   pivot_tols
		);
	///
	size_type max_size() const;
	///
	Inertia inertia() const;
	///
	void set_uninitialized();
	///
	void augment_update(
		const VectorSlice  *t
		,value_type        alpha
		,bool              force_refactorization
		,EEigenValType     add_eigen_val
		,PivotTolerances   pivot_tols
		);
	///
	void delete_update(
		size_type          jd
		,bool              force_refactorization
		,EEigenValType     drop_eigen_val
		,PivotTolerances   pivot_tols
		);

	// ///////////////////////////////////////////////////////
	// Overridden from MatrixSymWithOpFactorized

	///
	size_type rows() const;
	///
	std::ostream& output(std::ostream& out) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;

private:

	// /////////////////////////////////////////////////////
	// Private types

	typedef std::vector<FortranTypes::f_int> IPIV_t;

	// /////////////////////////////////////////////////////
	// Private data members

	size_type S_size_;      // The size of the current symmetric matrix.  Size == 0 if flag for uninitialized.
	bool      S_indef_;     // True if the current matrix is indefinite.
	bool      fact_updated_;// True if the factorization for the current matrix is updated.  Only meaningful
	                        // if S_indef_==true.
	bool      fact_in1_;    // If true then the current factorization is in S_store1_
	                        // otherwise it is in S_store2_.  Only meaningful if S_indef_==true and
	                        // fact_updated_==true
	MatrixSymAddDelUpdateable::Inertia
	          inertia_;     // Inertial for the indefinite L*D*L' factorization.  If fact_updated_ == false
	                        // then this will be UNKNOWN.  IF S_indef_==false then this is meaningless.
	GenMatrix S_store1_;    // Storage for the factored matrix in the
	                        // upper triangle as well as the original matrix
	                        // in the lower triangle.  This uses the same storage scheme as
	                        // in MatrixSymPosDefCholFactor.
	GenMatrix S_store2_;    // Storage for the factorization also.  We need
	                        // two storage locations for the L*D*L factorization
	                        // in case an update is singular.  This will not
	                        // be initialized for a p.d. or n.d. matrix.
	IPIV_t    IPIV_;        // Stores permutations computed by LAPACK
	mutable Vector
              WORK_;        // workspace
	MatrixSymPosDefCholFactor
		      S_chol_;      // Used to factorize the matrix
	                        // when it is p.d. or n.d.

	// /////////////////////////////////////////////////////
	// Private member funcitons.

	///
	/** Get view of DU.
	 */
	tri_ele_gms DU(size_type S_size, bool fact_in1);
	///
	const tri_ele_gms DU(size_type S_size, bool fact_in1) const;
	///
	/** Get view of lower part of S.
	 */
	sym_gms S(size_type S_size);
	///
	const sym_gms S(size_type S_size) const;
	///
	void assert_initialized() const;
	///
	void resize_DU_store( bool in_store1 );
	///
	/** Copy the original matrix into the new storage location and factorize it.
	 *
	 * Will throw LinAlgLAPack::FactorizationException if singular.
	 */
	void copy_and_factor_matrix( size_type S_size, bool fact_in1 );
	///
	/** Factor the current set matrix in-place (do not copy the original). 
	 *
	 * Will throw LinAlgLAPack::FactorizationException if singular.
	 */
	void factor_matrix( size_type S_size, bool fact_in1 );
	///
	/** Compute the new inertia and validate that it is what the client says it was.
	 *
	 * Will throw exceptions if the matrix is singular or has the wrong inertia.  If
	 * the matrix is near singular then true will be returned, the update should
	 * succeed but a warning exception should be thrown
	 */
	bool compute_assert_inertia(
		size_type S_size, bool fact_in1
		,const Inertia& expected_inertia, const char func_name[]
		,PivotTolerances pivot_tols, Inertia* comp_inertia, std::ostringstream* err_msg, value_type* gamma );

	/// Not defined and not to be called.
	MatrixSymAddDelBunchKaufman( const MatrixSymAddDelBunchKaufman& );
	MatrixSymAddDelBunchKaufman& operator=( const MatrixSymAddDelBunchKaufman& );

};	// end class MatrixSymAddDelBunchKaufman

// ///////////////////////////
// Inline member functions

inline
tri_ele_gms MatrixSymAddDelBunchKaufman::DU(size_type S_size, bool fact_in1)
{
	resize_DU_store(fact_in1);
	return LinAlgPack::nonconst_tri_ele(
		( fact_in1 ? S_store1_ : S_store2_ )(1,S_size,2,S_size+1)
		,BLAS_Cpp::upper );
}

inline
const tri_ele_gms MatrixSymAddDelBunchKaufman::DU(size_type S_size, bool fact_in1) const
{
	return LinAlgPack::tri_ele(
		( fact_in1 ? S_store1_ : S_store2_ )(1,S_size,2,S_size+1)
		,BLAS_Cpp::upper);
}

inline
sym_gms MatrixSymAddDelBunchKaufman::S(size_type S_size)
{
	return LinAlgPack::nonconst_sym(
		S_store1_(2,S_size+1,1,S_size)
		, BLAS_Cpp::lower );
}

inline
const sym_gms MatrixSymAddDelBunchKaufman::S(size_type S_size) const
{
	return LinAlgPack::sym(
		S_store1_(2,S_size+1,1,S_size)
		, BLAS_Cpp::lower );
}

}	// namespace ConstrainedOptimizationPack 

#endif	// MATRIX_SYM_POS_DEF_BUNCH_KAUFMAN_H
