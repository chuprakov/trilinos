// /////////////////////////////////////////////////////////////
// DirectSparseSolver.h
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

#ifndef	DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_H
#define DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_H

#include "SparseSolverPackTypes.h"
#include "SparseLinAlgPack/include/MatrixNonsingularSerial.h"
#include "SparseLinAlgPack/include/MatrixConvertToSparse.h"
#include "ref_count_ptr.h"
#include "AbstractFactory.h"

namespace SparseSolverPack {

///
/** Abstract interface to serial direct sparse linear solvers.
 *
 * This interface is designed to accommodate both unsymmetic and symmetric
 * direct solvers.  The motivation for using one common interface is that
 * we would like to be able to easily use an unsymmetic solver on a
 * symmetric system.  Of couse this would also allow a client to attempt
 * to solve an unsymmetric system with a symmetric solver which of course
 * will not work.  A symmetric implementation can of course check that
 * <tt>A.rows() == A.cols()</tt> but it would be much more difficult to determine
 * if the matrix was indeed symmetric.  It is likely that the symmetric
 * solver would only extract the upper or lower triangular region of an
 * unsymmetric matrix and would never check the other triangular region
 * for compatibility.
 *
 * The interface is designed to allow for maximum memory sharing
 * and recycling.  For example, we would like for the client to be able
 * to maintain more than one set of nonzero factors using the same sparsity
 * and fill-in structure determined by the analize_and_factorize(...)
 * operation.
 *
 * The basic idea is that given an <tt>m x n</tt> rectangular matrix \c A,
 * where <tt>m <= n</tt>, objects of this type will find
 * row and column permutations \c P and \c Q and a rank \c r such that
 * <tt>C = (P'*A*Q)(1:r,1:r)</tt> is the largest square nonsingular
 * (well conditioned) matrix that can be found.  This basis matrix
 * \c C is represented as a \c BasisMatrix object that supports the
 * \c MatrixNonsingularSerial interface.
 *
 * *** ToDo: Change this!
 * A smart pointer to this \c BasisMatrix
 * object is returned from the analyze_and_factor(...) and the
 * factor(...) methods.
 *
 * All the information needed to solve for a linear
 * system and to factorize other matrices with the same structure
 * is contained in the \c BasisMatrix object.  Therefore,
 * a \c DirectSparseSolver (\c DSS) object can
 * be considered to be stateless if needed.  The first point to make
 * clear is that once the client has a reference to a \c BasisMatrix object
 * in a smart pointer like \c C1_ptr, that \c BasisMatrix object will not be altered
 * by any future operations on the \c DSS object that created it unless the
 * client specifically requests an alteration.  This behavior allows \c BasisMatrix
 * objects to be completely decoupled from the \c DSS object that created them.
 * This behavior is explained in more detail later.
 *
 * It is through the \c MatrixNonsingularSerial interface of the \c BasisMatrix object
 * that clients can solve for linear systems.
 *
 * The usage of this interface will be illustrated by considering
 * a few different scenarios:
 *
 * 1) Let the DSS object manage the BasisMatrix object as we only
 * need to maintian the factorization of one unsymmetric matrix at a time.
 *
 \code

	// Storage for the permutations and rank
	IVector row_perm, col_perm;
	size_type rank;
	
	// Analyze and factor A1
	direct_solver.analyze_and_factor( A1, &row_perm, &col_perm, &rank );

	// Solve for x = inv(C1)*b1
	Vector x;
	V_InvMtV( &x, *direct_solver.get_basis_matrix(), BLAS_Cpp::no_trans, b1 );

	// Factor another matrix with the same structure
	direct_solver.factor( A2 );

	// Solve for x = inv(C2')*b2
	V_InvMtV( &x, *direct_solver.get_basis_matrix(), BLAS_Cpp::trans, b2 );

 \endcode
 *
 * 2) Maintain the factorization of two unsymmetic matrices with the same
 *	structure at the same time.
 *
 \code
	
	// Storage for the permutations and rank
	IVector row_perm, col_perm;
	size_type rank;

	// Factorize A1
	typedef DirectSparseSolver::basis_matrix_ptr_t C_ptr_t;
	C_ptr_t
		C1_ptr = direct_solver.analize_and_factor( A1, &row_perm, &col_perm, &rank );

	// Solve for x = inv(C1)*b1
	Vector x;
	V_InvMtV( &x, *C1_ptr, BLAS_Cpp::no_trans, b1 );

	// Factor another matrix with the same structure but new nonzeros.
	// The fact that C1_ptr still points to the previous BasisMatrix
	// ensures that this operation will not overwrite its nonzeros.
	C_ptr_t
		C2_ptr = direct_solver.factor( A2 );

	// Solve for x = inv(C2)*b2
	V_InvMtV( &x, *C2_ptr, BLAS_Cpp::no_trans, b2 );

	// Recycle the storage for the factorization of A1 and factorize A3
	C1_ptr = direct_solver.factor( A3, C1_ptr->get_fact_struc(), C1_ptr );

	// Solve for x = inv(C3)*b3
	V_InvMtV( &x, *C1_ptr, BLAS_Cpp::no_trans, b3 );

 \endcode
 *
 * 3) Maintain the factorization of three unsymmetic matrices, two with the
 * same structure and one with a different structure.
 *
 \code
	
	// Storage for the permutations and rank
	IVector row_perm1, col_perm1;
	size_type rank1;
	IVector row_perm2, col_perm2;
	size_type rank2;

	// Factorize A1
	typedef DirectSparseSolver::basis_matrix_ptr_t C_ptr_t;
	C_ptr_t
		C1_ptr = direct_solver.analize_and_factor( A1, &row_perm1, &col_perm1, &rank1 );

	// Solve for x = inv(C1)*b1
	Vector x;
	V_InvMtV( &x, *C1_ptr, BLAS_Cpp::no_trans, b1 );

	// Factor another matrix with a different structure.
	C_ptr_t
		C2_ptr = direct_solver.analize_and_factor( A2, &row_perm2, &col_perm2, &rank2 );

	// Solve for x = inv(C2)*b2
	V_InvMtV( &x, *C2_ptr, BLAS_Cpp::no_trans, b2 );

	// Factor another matrix A3 with the same structure of as A1 but preserve
	// the factorization nonzeros of A1.
	C_ptr_t
		C3_ptr = direct_solver.factor( A3, C1_ptr->get_fact_struc() );

	// Solve for x = inv(C3)*b3
	V_InvMtV( &x, *C3_ptr, BLAS_Cpp::no_trans, b3 );

 \endcode
 *
 * 4) Factor of two unsymmetic matrices with different structure and then
 * recycle the storage of the first matrix to factor a third matrix with
 * an all new structure.
 *
 \code
	
	// Storage for the permutations and rank
	IVector row_perm1, col_perm1;
	size_type rank1;
	IVector row_perm2, col_perm2;
	size_type rank2;

	typedef DirectSparseSolver::basis_matrix_ptr_t
		C_ptr_t;

	// Factorize A1
	C_ptr_t
		C1_ptr = direct_solver.analize_and_factor( A1, &row_perm1, &col_perm1, &rank1 );

	// Solve for x = inv(C1)*b1
	Vector x;
	V_InvMtV( &x, *C1_ptr, BLAS_Cpp::no_trans, b1 );

	// Factor another matrix with a different structure.
	C_ptr_t
		C2_ptr = direct_solver.analize_and_factor( A2, &row_perm2, &col_perm2, &rank2 );

	// Solve for x = inv(C2)*b2
	V_InvMtV( &x, *C2_ptr, BLAS_Cpp::no_trans, b2 );

	// Analyze and factor another matrix A3 with an all new structure
	// and recycle storage of C1.
	C1_ptr = direct_solver.analyze_and_factor( A3, &row_perm1, &col_perm1, &rank1
		, C1_ptr );

	// Solve for x = inv(C3)*b3
	V_InvMtV( &x, *13_ptr, BLAS_Cpp::no_trans, b3 );

	// Factor another matrix A4 with the same structure of as A2 but preserve
	// the factorization of A2.
	C_ptr_t
		C4_ptr = direct_solver.factor( A4, C2_ptr->get_fact_struc() );

	// Solve for x = inv(C4)*b4
	V_InvMtV( &x, *C4_ptr, BLAS_Cpp::no_trans, b4 );

	// Analyze and factor another matrix A5 with an all new structure
	// and recycle storage of C4.  Note that the behavoir of C2_ptr must
	// not change by this operation.
	C4_ptr = direct_solver.analyze_and_factor( A5, &row_perm5, &col_perm5, &rank5
		, C4_ptr );

 \endcode
 *
 * In general the client should not try to manually modify the smart pointers
 * returned by any of the objects involved in this interface or try to modify
 * the objects they point to except through this interface.  To help clients
 * avoid problems associated with accidentlly modifying these smart pointers,
 * all of the references to the smart points are declared constant.
 * The client should not cast this constantness away and fiddle with these
 * const references.  The objects being pointed to themselves are not declared
 * const since there are no modifying operations on these objects anyway.
 *
 * If by some major mistake the client were to pass in basis_matrix
 * or fact_struc objects of the incorrect type from non-compatible
 * DSS objects then the InvalidObjectType exception would be thrown.
 * If the client follows the protocal defined in this interface, this will
 * never happen.
 */
class DirectSparseSolver {
public:

	/** @name Public Types */
	//@{

	class FactorizationStructure;

	///
	/** Abstract class for objects that represent the factorized
	  * matrix and can be used to solve for different right-hand-sides.
	  *
	  * This object encapsulates the factorzation structure and the nonzero
	  * values of the factorized basis matrix.
	  */
	class BasisMatrix : public MatrixNonsingularSerial {
	public:
		///
		typedef ReferenceCountingPack::ref_count_ptr<FactorizationStructure>    fact_struc_ptr_t;
		/// Return a smart pointer to the object that represents the factorization structure.
		virtual const fact_struc_ptr_t&		get_fact_struc() const = 0;
	};

	///
	/** Abstract class for objects that represent the factorization
	  * structure of a particular matrix.
	  *
	  * This structure can be reused over and over again for factorizing matrices
	  * with the same structure but different nonzero elements.
	  */
	class FactorizationStructure {
	public:
		///
		virtual ~FactorizationStructure() {}
	};

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractFactoryPack::AbstractFactory<BasisMatrix> >   basis_matrix_factory_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<BasisMatrix>       basis_matrix_ptr_t;

	///
	class UnsymmetricRankDeficientException : public std::logic_error
	{public: UnsymmetricRankDeficientException (const std::string& what_arg)
	: std::logic_error(what_arg) {}};

	///
	class IncompatibleMatrixStructureException : public std::logic_error
	{public: IncompatibleMatrixStructureException (const std::string& what_arg)
	: std::logic_error(what_arg) {}};

	///
	class InvalidObjectType : public std::logic_error
	{public: InvalidObjectType (const std::string& what_arg)
	: std::logic_error(what_arg) {}};

	//@}
	
	///
	virtual ~DirectSparseSolver() {}

	///
	/** Return a factory object that can create the basis matrix.
	 *
	 *
	 */
	virtual const basis_matrix_factory_ptr_t basis_matrix_factory() const = 0;

	///
	/** Set the estimate of the fill-in ratio of the factorization.
	 *
	 * This value is used to set the initial guess at the amount
	 * of storage needed for the factorization of the matrix
	 * plus some "elbow" room for the factorization algorithm.
	 *
	 * The estimated amount of nonzeros for the factorization is:
	 *
	 * num_factor_nonzeros = estimated_fillin_ratio * A.nz()
	 */
	virtual void estimated_fillin_ratio( value_type estimated_fillin_ratio ) = 0;

	///
	/** Analyze and factor a matrix.
	 *
	 * This function extracts the nonzero elements from a matrix and
	 * then finds a nonsingular basis matrix and factors it (see intro).
	 *
	 * Given the <tt>m x n</tt> matrix \c A (<tt>m <= n</tt>), this method finds the
	 * row and column permutations \c P and \c Q respectively and the rank \c r such that
	 * the matrix <tt>C = <tt>(P'*A*Q)(1:r,1:r)</tt> is nonsingular.
	 *
	 * If a symmetric system is being solved and it is feared that \c A is not full rank
	 * and the client wants to use symmetric pivoting (<tt>P = Q</tt>) to find the basis
	 * \c C, then the client can pass in NULL \c for \c col_perm and force the solver to
	 * return a symmetric permutation.  Note that the row and column permutations of
	 * \c row_perm(k) and \c col_perm(k) for <tt>k = 1...r</tt> does not neccesarly indicate
	 * the exact permutations used within the sparse solver.  Instead, what is significant is
	 * partitioning of rows and columns between <tt>k = 1..r</tt> and <tt>k = r+1..m</tt>
	 * for row permutations and <tt>k = r+1..n</tt> for column permutations.  Therefore,
	 * if <tt>col_perm==NULL</tt> on input and an unsymmetric solver is used and the matrix
	 * is not full rank, then since in reality internally unsymmetric pivoting is being used,
	 * the implementation must throw an \c UnsymmetricRankDeficientException exception.
	 * This convention is needed to allow an unsymmetric solver to be used for a symmetric
	 * system but still meet the needs of the user for symmetric systems.  In summary, as long
	 * as the symmetric system is full rank, then an unsymmetic solver can always be used and
	 * the client need not care how the system is solved.  On the other hand, if the client does
	 * not care if unsymmetric pivoting is used on a symmetric rank deficient system, then
	 * the client can set \c col_perm to point to a valid \c IVector object and deal
	 * with the unsymmetic permutations that are returned.  Who knows what the use of this
	 * would be but there is no reason not to allow it.
	 *
	 * This method allows for the carefull sharing and recycling of the factorization structure
	 * and nonzeros of a factorized matrix
	 * while also maintaining a "current" factorization.  To explain
	 * how this works let C1_ptr be set to the smart point to the basis
	 * matrix returned by analyze_and_factor(...), factor(...)
	 * or get_basis_matrix(...).
	 * 
	 * If the client whats to recycle the factorization structure and
	 * nonzeros that are contained in C1_ptr to anylize and
	 * factor A3 then the client would perform:
	 *
	 \code
	 C1_ptr = direct_solver.analyze_and_factor(A3,&row_perm,&col_perm,&rank,C1_ptr);
	 \endcode
	 * 
	 * It is perfectly legal in the above example for C1_ptr==this->get_basis_matrix()
	 * before this->analyze_and_factor(...,C1_ptr) is called.
	 *
	 * Note that once C1_ptr is passed to this->analyze_and_factor(...,C1_ptr)
	 * that it will be transformed such that after this->analyze_and_factor(...,C1_ptr)
	 * returns that this->get_basis_matrix() == C1_ptr.  In this case,
	 * it is not really needed to set C1_ptr = analyze_and_factor(...,C1_ptr)
	 * like in the above statement but is makes explicit that C1_ptr changes.
	 * If C1_ptr->get_fact_struc()->count() > 1 before this method is called,
	 * then the client has a reference to another basis matrix that also uses
	 * this factorization structure.  In this case a new factorization structure
	 * must be allocated as not to disturb the other basis matrices.
	 *
	 * If basis_matrix.get()==NULL and this->get_basis_matrix().get()!=NULL and
	 * this->get_basis_matrix().count()==1 then no new storage is computed
	 * and the internally referenced storage is recycled automatically.  However if
	 * basis_matrix.get()==NULL and this->get_basis_matrix().get()==NULL
	 * then no basis has be factored yet and new storage for structure and nonzeros
	 * must be allocated.  If basis_matrix.get()==NULL, this->get_basis_matrix().get()!=NULL
	 * and this->get_basis_matrix().count() > 1, then the client has maintained a reference
	 * to the internal basis matrix.  In this case #this# must not alter that
	 *	basis matrix object so new storage for the factoraztion structure and nonzeros
	 * must be computed for the current matrix being analyzed and factored.
	 *
	 * Preconditions: \begin{itemize}
	 *		\item #A.rows() <= A.cols()# (throw #std::invalid_argument#)
	 *		\end{itemize}
	 *
	 * Postconditions: \begin{itemize}
	 *		\item #get_basis_matrix() == C_ptr# where #C_ptr# was just returned
	 *			from #analyze_and_factor(...)#
	 *		\end{itemize}
	 *
	 *	@param	A			[in] Matrix who's elements are extracted in order
	 *						to form the factorization.  Here A.rows() must
	 *						be <= A.cols().
	 * @param	A_trans		[in] If A_trans == no_trans then op(A) = A.  If
	 *						A_trans == trans then op(A) = A'.
	 *	@param	row_perm	[out] On output holds an array (size A.rows()) of row
	 *						permutations (P) used to select the basis matrix C.
	 *						The ith row of P*A*Q is row (*row_perm)(i) of A.
	 *	@param	col_perm	[out] On output holds an array (size A.cols()) of column
	 *						permutations (Q) used to select the basis matrix C.
	 *						The jth column of P*A*Q is column (*col_perm)(j) of A.
	 *						If the client is solving a symmetric system then
	 *						col_perm may be set to NULL on input which forces
	 *						the implementation to use symmetric permutations
	 *						(see above).
	 *	@param	rank		[out] On output holds the rank (r) (r <= A.rows()) of the
	 *						basis matrix C.
	 *	@param	basis_matrix
	 *						[in/out] This is used to recycle storage for the structure
	 *						of the factorization and its nonzeros.
	 *						The default value is basis_matrix.get() == 0
	 *						in which no recycling will be performed and
	 *						references to previously computed basis matrices
	 *						will not be altered.  If basis_matrix.get() != 0 then
	 *						the referenced basis matrix will be transformed
	 *						such that basis_matrix == this->get_basis_matrix()
	 *						after this method successfully returns.
	 *	@param	out			[in/out] If out != NULL, then some information about the
	 *						operations performed internally may be printed
	 *						to *out.  The amount of this output should be
	 *						very minimal and should not significantly increase
	 *						with the size of the problem being solved.
	 *
	 *	@return
	 *		Returns the factorized basis matrix C for the input matrix A.
	 */
	virtual void analyze_and_factor(
		const MatrixConvertToSparse     &A
		,IVector                        *row_perm
		,IVector                        *col_perm
		,size_type                      *rank
		,BasisMatrix                    *basis_matrix   = NULL
		,std::ostream                   *out            = NULL
		) = 0;
	
	///
	/** Factor a matrix given its factorization structure.
	 *
	 * This method allows clients to factor a matrix given it predetermined
	 * factorization structure.  The client can use the factorization
	 * structure stored internally from the last call to this->analyze_and_factor(...)
	 * (fact_struc.get() == NULL) or use the factorization structure from another
	 * precomputed basis matrix (fact_struc.get() != NULL).  If the passed in matrix
	 * A is obviously not compatible with the precomputed factorization structure
	 * being used then an IncompatibleMatrixStructureException exception will be
	 * thrown.  If A1 was the matrix originally passed to this->analyze_and_factor(A1,...)
	 * and A2 is the matrix being passed to this->factor(A2,...) then if
	 * A2.rows() != A1.rows() or A2.cols() != A1.cols() or A2.nz() != A1.nz()
	 * then the matrices are obviously not compatible and the exception will
	 * be thrown.
	 *
	 * The argument basis_matrix is included to allow for the recycling of
	 * the nonzero values of the factorization of an existing basis matrix.
	 * For example, suppose that C1_ptr and C2_ptr point to two different
	 * basis matrices with different structure (i.e. C1_ptr->get_fact_struc() !=
	 * C2_ptr->get_fact_struc()).  Now suppose we would like to factorize
	 * another matrix A3 that has the same structure of A2 by reusing the nonzero
	 * values referenced in C1_ptr.  To do this the client would call:
	 *
	 \code
	 C1_ptr = direct_solver.factor(A3,C2_ptr->get_fact_struc(),C1_ptr);
	 \endcode
	 *
	 * Note that in the above example that setting C1_ptr = direct_solver.factor(...,C1_ptr)
	 * is not really necessary since C1_ptr will be modified anyway but it is done
	 * to explicity show this modification.
	 *
	 * Now consider the case where we would like to recycle the storage
	 * for the nonzeros in C1_ptr and also use the same structure as C1.
	 * To do this the client would call:
	 *
	 \code
	 C1_ptr = direct_solver.factor(A3,C1_ptr->get_fact_struc(),C1_ptr);
	 \endcode
	 *
	 * Above, the client must explicitly pass C1_ptr->get_fact_struc() to use
	 * this structure or the current interally stored structure given by
	 * this->get_basis_matrix()->get_fact_struc() would be used instead.
	 *
	 *	@param	A			[in] Matrix to be factored given a previously determined
	 *						factorization structure (see fact_struc).
	 * @param	A_trans		[in] If A_trans == no_trans then op(A) = A.  If
	 *						A_trans == trans then op(A) = A'.
	 *	@param	fact_struc	[in] The factorization structure to use.  If fact_struc.get()==NULL
	 *						then this->get_basis_matrix()->get_fact_struc() is
	 *						used in its place.
	 *	@param	basis_matrix
	 *						[in/out] The basis matrix used to recycle the nonzero
	 *						values of a previously computed basis matrix.
	 *						If basis_matrix.get()==NULL then new storage
	 *						for nonzeros will be allocated if
	 *						this->get_basis_matrix().count() > 1.  On output
	 *						if basis_matrix.get()!=NULL then the underlying
	 *						BasisMatrix object will be transformed to the
	 *						newly factorized basis matrix.
	 *	@param	out			[in/out] If out != NULL, then some information about the
	 *						operations performed internally may be printed
	 *						to *out.  The amount of this output should be
	 *						very minimal and should not significantly increase
	 *						with the size of the problem being solved.
	 *
	 *	@return
	 *		Returns the factorized basis matrix C for the input matrix A.
	 */
	virtual void factor(
		const MatrixConvertToSparse              &A
		,const BasisMatrix::fact_struc_ptr_t     &fact_struc    = ReferenceCountingPack::null
		,BasisMatrix                             *basis_matrix  = NULL
		,std::ostream                            *out           = NULL
		) = 0;

	///
	/** Get the current basis matrix.
	  *
	  * This returns a smart reference counted pointer to the basis matrix
	  * object returned by the last call to this->analize_and_factor(...)
	  * or this->factor(...).
	  */
	virtual const basis_matrix_ptr_t& get_basis_matrix() const = 0;

	///
	/** Release all allocated memory.
	  *
	  * After this function returns then get_basis_matrix().get() == 0.
	  */
	virtual void release_memory() = 0;

};	// end class DirectSparseSolver 

}	// end namespace SparseSolverPack 

#endif	// DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_H
