// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystem.h

#ifndef DECOMPOSITION_SYSTEM_H
#define DECOMPOSITION_SYSTEM_H

#include <stdexcept>
#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** This class abstracts a decomposition choose for the
  * range space #Y#, and null space #Z#, matrices for the linearly
  * independent set of columns of #A# {abstract}.
  *
  *	#A = [ A(indep),  A(dep)  ]#
  *
  * where A is n x m, A(indep) is n x r and A(dep) is n x (m - r)
  *
  * Note that the columns in #A(dep)# may be linearly dependent with
  * the columns in #A(indep)# but they may just be undecomposed
  * linearly independent equality constraints.
  *
  * The decomposition formed by subclasses must have the properties:
  *
  * #A(indep)' * Z = 0# (orthoganol property)\\
  * #[Y , Z]# is nonsingular\\
  *
  * This object forms the decomposition and solves the linear systems
  * #op(A(indep)' * Y) * x = b#.
  *
  * It is assumed that the matrix #A# is not modified between calls to
  * #update_decomp(...)#.
  */
class DecompositionSystem {
public:

	/** @name Exceptions */
	//@{

	/// Thrown if #this->solve_transAtY()# is called before #this->update_decomp()#
	class NoDecompositionSelected : public std::logic_error
	{public: NoDecompositionSelected(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if the matrix type for A, Z, Y, U, or V is not of the proper type
	class InvalidMatrixType : public std::logic_error
	{public: InvalidMatrixType(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}	

	// Defined to satisfy VC++ 5.0
	DecompositionSystem() {}
		
	///
	/** Creates the decomposition for #A'#.
	  *
	  * The decomposition is based on the linearly independent
	  * columns of A:\\
	  *
	  * #A = [ A(indep)  A(dep) ]#\\
	  *
	  * Specifically this operation finds the matrices:
	  *
	  * #Z# s.t. #A(indep)' * Z = 0#\\
	  * #Y# s.t. #[Z  Y]# is nonsingular\\
	  * #U = A(dep)' * Y#\\
	  * #V = A(dep)' * Z#\\
	  *
	  * After this operations is called then the MatrixWithOp view of A will
	  * be sorted so that #A = [ A(indep)  A(dep) ]#.
	  *
	  * If there is some problem creating the decomposition then exceptions
	  * with the base class std::exception may be thrown.  The meaning
	  * of these exceptions are more associated with the subclasses
	  * that implement this operation.
	  *
	  * The concrete types for #A#, #Z#, #Y#, #U#, and #V# must be compatable with
	  * the concrete subclass or an #InvalidMatrixType# exeption will be thrown..
	  *
	  * If m() == r() then on output #U# and #V# should not be accessed.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #A.rows() > A.cols()# (throw std::length_error)
	  * \end{itemize} 
	  *
	  * Postconditions:\begin{itemize}
	  * \item #Z.rows() == r()#
	  * \item #Z.cols() == n() - this->r()#
	  * \item #Y.rows() == r()#
	  * \item #Y.cols() == n()#
	  * \item [#r() < m()#] #U.rows() == m() - r()#
	  * \item [#r() < m()#] #U.cols() == r()#
	  * \item [#r() < m()#] #V.rows() == m() - r()#
	  * \item [#r() < m()#] #V.cols() == n() - r()#
	  * \item Z is orthoganal to A(indep)'
	  * \item [Y  Z] is nonsingular
	  * \end{itemize}
	  *
	  * @param	A	[I/O]	On input represents the matrix to have its
	  *						decomposition formed.  On output its columns
	  *						are conceptually partitioned into independent
	  *						and dependent parts where A(:,indep) is n x r
	  *						and A(:,dep) is n x (m-r).  The matrix
	  *						MatrixWithOp interface for A supports this.
	  *	@param	Z	[O]		On output represents the n x (n-r) null space
	  *						matrix such that (A' * Z)(indep,:) == 0.
	  *	@param	Y	[O]		On output represents the n x r Range space matrix
	  *						such that [Y Z] is nonsingular.
	  *	@param	U	[O]		On output represents the (m-r) x r matrix
	  *						A(:,dep) * Y.
	  *	@param	V	[O]		On output represents the (m-r) x (n-r) matrix
	  *						A(:,dep) * Z
	  */
	virtual void update_decomp(MatrixWithOp* A, MatrixWithOp* Z, MatrixWithOp* Y
		, MatrixWithOp* U, MatrixWithOp* V) = 0;

	///
	/** Solves the linear systems #op(A' * Y) * x = b#.
	  *
	  * This operation solves the linear system:\\
	  *
	  * (A' * Y) * x = b,  for #trans == BLAS_Cpp::no_trans#\\
	  * (Y' * A) * x = b,  for #trans == BLAS_Cpp::trans#\\
	  *
	  * Here #x# and #b# can be the same.
	  *
	  * Preconditions (A is matrix returned in the last call to this->update_decomp(..)):
	  * \begin{itemize}
	  * \item #b.size() == A.cols()# (throw std::length_error)
	  * \item #x.size() == A.cols()# (throw std::length_error)
	  * \end{itemize} 
	  */
	virtual void solve_transAtY(const VectorSlice& b, BLAS_Cpp::Transp trans
		, VectorSlice* x) const = 0;

	///
	/** Solves the set of linear systems #op(A' * Y) * X = B#.
	  *
	  * This operation solves the set of linear systems:\\
	  *
	  * (A' * Y) * X = B,  for #trans == BLAS_Cpp::no_trans#\\
	  * (Y' * A) * X = B,  for #trans == BLAS_Cpp::trans#\\
	  *
	  * Here #X# and #B# can be the same.
	  *
	  * Preconditions (A is matrix returned in the last call to this->update_decomp(..)):
	  * \begin{itemize}
	  * \item #B.rows() == A.cols()# (throw std::length_error)
	  * \item #X.rows() == A.cols()# (throw std::length_error)
	  * \item #X.cols() == B.cols()# (throw std::length_error)
	  * \end{itemize} 
	  */
	virtual void solve_transAtY(const GenMatrixSlice& B, BLAS_Cpp::Transp trans
		, GenMatrixSlice* X) const = 0;

	/// Return the number of columns in #A'#
	virtual size_type n() const = 0;

	/// Return the number of rows in #A'#
	virtual size_type m() const = 0;

	/// Returns the rank of #A(indep)#
	virtual size_type r() const = 0;

	///
	/** Setup of the trase of calculations and checking of errors.
	  *
	  * @param	out		[O]	If out == 0 then no output is produced/
	  *	@param	trase	[I]	If true then calculations and tests
	  *						will be output to out.
	  *	@param	check_results	[I]	If true then results will be checked and
	  *								if any fails then an std:runtime_error
	  *								exception will be thrown
	  */
	virtual void setup_trase( std::ostream* out, bool trase = true
		, bool check_results = false ) = 0;

};	// end class DecompositionSystem

}	// end namespace ConstrainedOptimizationPack

#endif	// DECOMPOSITION_SYSTEM_H
