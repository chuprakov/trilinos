// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystemVarReduct.h

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_H

#include "DecompositionSystem.h"
#include "SparseSolverPack/include/BasisSystem.h"

namespace ConstrainedOptimizationPack {

///
/** Specialization of DecompositionSystem for Variable Reduction basises {abstract}.
  *
  * This interface abstracts a variable reduction decomposition where:
  *
  \begin{verbatim}
  
  P*A'*Q = [C  N] 
           [E  F]

  Z = [ D ]
      [ I ]

      where:   D = -inv(C) * N

  A(indep)'*Y = C

  \end{verbatim}
  *
  * Here #C# is a r x r nonsingular matrix.
  */
class DecompositionSystemVarReduct : public DecompositionSystem {
public:

	/** @name types */
	//@{

	///
	typedef DecompositionSystem	inherited;

	//@}

	/** @name Exceptions */
	//@{

	/// Thrown if the basis given in set_basis(..) is structurally singular
	typedef BasisSystem::StructuralSingularityException	StructuralSingularityException;

	/// Thrown if the basis given in set_basi(...) is numerically singular
	typedef BasisSystem::NumericalSingularityException	NumericalSingularityException;

	/// Thrown if an invalid basis is passed to #this->set_basis()#
	typedef BasisSystem::InvalidBasis					InvalidBasis;								

	/// Thrown if #this->get_basis()# is called before #this->set_basis()#, or #this->update_decomp()#
	typedef BasisSystem::NoBasisSelected				NoBasisSelected;								
	
	/// Thrown if #this->solve_C()#, or #this->N()# is called before #this->update_decomp()#
	typedef inherited::NoDecompositionSelected			NoDecompositionSelected;

	///
	typedef inherited::InvalidMatrixType				InvalidMatrixType;

	//@}	

	/// Return the current basis being used A.
	virtual void get_basis(IVector* row_perm, IVector* col_perm
		, size_type* rank) const = 0;

	///
	/** Solves the systems #op(C) * x = b#.
	  *
	  * This operation solves the linear system:\\
	  *
	  * #C * x = b#,  for #trans == BLAS_Cpp::no_trans#\\
	  * #C' * x = b#,  for #trans == BLAS_Cpp::trans#\\
	  *
	  * Preconditions (rank as returned from the last call to get_basis(....,&rank)):
	  * \begin{itemize}
	  * \item #b.size() == rank# (throw std::length_error)
	  * \item #x.size() == rank# (throw std::length_error)
	  * \end{itemize} 
	  */
	virtual void solve_C(const VectorSlice& b, BLAS_Cpp::Transp trans
		, VectorSlice* x) const = 0;

	///
	/** Solves the set of linear systems #op(C) * X = B#.
	  *
	  * This operation solves the set of linear systems:\\
	  *
	  * #C * X = B#,  for #trans == BLAS_Cpp::no_trans#\\
	  * #C' * X = B#,  for #trans == BLAS_Cpp::trans#\\
	  *
	  * Preconditions (rank as returned from the last call to get_basis(....,&rank)):
	  * \begin{itemize}
	  * \item #B.rows() == rank# (throw std::length_error)
	  * \item #X.rows() == rank# (throw std::length_error)
	  * \item #X.cols() == B.cols()# (throw std::length_error)
	  * \end{itemize} 
	  */
	virtual void solve_C(const GenMatrixSlice& B, BLAS_Cpp::Transp trans
		, GenMatrixSlice* X) const = 0;

	/// Return a reference to the C matrix
	virtual const MatrixWithOp& C() const = 0;

	/// Return a reference to the N matrix
	virtual const MatrixWithOp& N() const = 0;

	/// Return a reference to the E matrix if it exists
	virtual const MatrixWithOp& E() const = 0;

	/// Return a reference to the F matrix if it exists
	virtual const MatrixWithOp& F() const = 0;

};	// end class DecompositionSystemVarReduct

}	// end namespace ConstrainedOptimizationPack

#endif	// DECOMPOSITION_SYSTEM_VAR_REDUCT_H