// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxed.h

#ifndef QP_SOLVER_RELAXED_H
#define QP_SOLVER_RELAXED_H

#include "QPSolverStats.h"

namespace ConstrainedOptimizationPack {

///
/** Solves Quadratic Programs (QPs) of several different forms while
  * relaxing the constraints.
  *
  * The formulation for the QP being solve is:
  \begin{verbatim}
  
  (1.1)  min          g'*d + 1/2*d'*G*d + M(eta)
         d <: R^n
         
         s.t.
  (1.2)               etaL <=  eta
  (1.3)               dL   <=  d                       <= dU
  (1.4)               eL   <=  op(E)*d - b*eta         <= eU
  (1.5)                        op(F)*d + (1 - eta) * f  = 0
  
  \end{verbatim}
  * The relaxation is used to ensure that the QP will have a solution
  * (eta = 1, d = 0 guaranteed if dL <= 0 <= dU and eL <= b <= eU).
  * If the function M(eta) in the objective is large enough, then
  * the constraint etaL <= eta will be active if a feasible region
  * exists.  The form of the function M(eta) is not specified by this
  * interface defined as appropriate for each individual QP solver
  * method and implementation. 
  *
  * The Lagrangian for the QP in (1) is:
  \begin{verbatim}
  
	L = g' * d + 1/2 * d' * G * d + M(eta)
	     + kappa * (etaL - eta)
	     + nuL' * (dL - d)
	     + nuU' * (d - dU)
	     + muL' * (eL - op(E)*d - b*eta)
	     + muU' * (op(E)*d - b*eta - eU)
		 + lambda' * (op(F)*d + (1 - eta) * f)
	
  \end{verbatim}
  * The optimality conditions for this QP are:
  \begin{verbatim}

	Linear dependence of gradients:
	  
	(2)  d(L)/d(d) = g + G*d - nuL + nuU + op(E)'*(- muL + muU) + op(F)'*lambda
	               = g + G*d + nu + op(E)'*mu + op(F)'*lambda = 0
	  
	     where: nu = nuU - nuL, mu = muU - muL

	(3)  d(L)/d(eta) = d(M)/d(eta) - kappa - b'*mu - f'*lambda = 0
	  
	Feasibility:
	  
	(4.1)  etaL <=  eta
	(4.2)  dL   <=  d                       <= dU
	(4.3)  eL   <=  op(E)*d - b*eta         <= eU
	(4.4)  op(F)*d + (1 - eta) * f  = 0

	Complementarity:

	(5.1)  nu(i) * (dL - d)(i) \approx 0, if nu(i) <= 0, i = 1...n
	(5.2)  nu(i) * (d - dU)(i) \approx 0, if nu(i) >= 0, i = 1...n
	(5.3)  mu(j) * (eL - op(E)*d + b*eta)(j) \approx 0, if mu(j) <= 0, j = 1...m_in
	(5.4)  mu(j) * (op(E)*d - b*eta - eU)(j) \approx 0, if mu(j) >= 0, j = 1...m_in

	Nonnegativity of Lagrange Multipliers for Inequality Constraints:

	(6.1)  nu(i) <= 0 if (dL - d)(i) \approx 0, i = 1...n
	(6.2)  nu(i) >= 0 if (d - dU)(i) \approx 0, i = 1...n
	(6.3)  mu(j) <= 0 if (eL - op(E)*d - b*eta)(j) \approx 0, j = 1...m_in
	(6.4)  mu(j) >= 0 if (op(E)*d - b*eta - eU)(j) \approx 0, j = 1...m_in

  \end{verbatim}
  * The optimal #d# and #eta# are determined as well as the lagrange multipliers
  * for the constriants:
  \begin{verbatim}

	nu :       dL <= d <= dU
 
	mu :       eL <= op(E)*d - b*eta <= eU

	lambda :   op(F)*d + (1 - eta) * f  = 0

  \end{verbatim}
  * The lagrange multiper #kappa# for the constraint #etaL <= eta# is not returned
  * since if this constraint is not active, then #kappa == 0# and all of the multiplier
  * estimates will be off because of the arbitrarily large value of #d(M)/d(eta)# in
  * the optimality condition (3).
  *
  * The bounds #dL#, #dU#, #eL# and #eU# are input as sparse vectors (SpVectorSlice).
  *
  * Note that this interface is geared toward active-set solvers but could also
  * be used for iterative solvers with care.
  */
class QPSolverRelaxed {
public:

	/** @name Public Types */
	//@{

	/// Thrown if the QP is unbounded.
	class Unbounded : public std::logic_error
	{public: Unbounded(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if the QP is infeasible.
	class Infeasible : public std::logic_error
	{public: Infeasible(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if there is invalid input
	class InvalidInput : public std::logic_error
	{public: InvalidInput(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if a test failed
	class TestFailed : public std::logic_error
	{public: TestFailed(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Enumeration for the amount of output to create from solve_qp(...).
	enum EOutputLevel {
		PRINT_NONE			= 0,
		PRINT_BASIC_INFO	= 1,
		PRINT_ITER_SUMMARY	= 2,
		PRINT_ITER_STEPS	= 3,
		PRINT_ITER_ACT_SET	= 4,
		PRINT_ITER_VECTORS	= 5,
		PRINT_EVERY_THING	= 6
		};

	/// Enumeration for if to run internal tests or not.
	enum ERunTests { RUN_TESTS, NO_TESTS };

	//@}

	///
	virtual ~QPSolverRelaxed() {}

	///
	/** Solve the QP.
	  *
	  * This function may throw many exceptions.  If there is some problem with
	  * the QP definition the exceptions #Unbounded#, #Infeasible# or #InvalidInput#
	  * may be thrown.  If the QP is illconditioned other exeptions may be thrown
	  * by this function or perhaps warning messages may be printed to #*out# and
	  * a value other than #OPTIMAL_SOLUTION# may be returned.
	  *
	  * After this function returns, #this->get_qp_stats()# can be called to return
	  * some statistics for the QP just solved (or attempted to be solved).
	  *
	  * Note, the variable bounds can be removed by passing in #dL.nz() == dU.nz() == 0#.
	  * 
	  * By default this function calls the function #this->solve_qp(...)# which accepts
	  * various sets of constraints.
	  *
	  *	@param	out			[out] If out != NULL then output is printed to this stream
	  *							depending on the value of #olevel#.
	  *	@param	olevel		[in] Determines the amount of output to print to *out.
	  *							The exact type of output is determined by the implementing
	  *							subclass but here is the sugguested behavior:
	  *							PRINT_NONE : Don't print anything (same as out == NULL).
	  *							PRINT_BASIC_INFO : Only print basic information about
	  *								the solution of the QP.  Amount of output = O(1).
	  *							PRINT_ITER_SUMMARY : Prints a summary of each iteration
	  *								in the algorithm.  Amount of output = O(num_iter).
	  *							PRINT_ITER_STEPS : Prints output about each iteration
	  *								in more detail than PRINT_ITER_SUMMARY but is still
	  *								O(num_iter).
	  *							PRINT_ITER_VECTORS : Prints out important vectors computed
	  *								in each QP iteration.  Mainly useful for debugging.
	  *								Amount of output = O((num_iter)(n)).
	  *							PRINT_EVERY_THING : Print out nearly every important
	  *								quantity within each QP iteration including
	  *								vectors and matrices.  Mainly useful for debugging.
	  *								Amount of output = O((num_iter)(n-m)^2)+O((num_iter)(n))
	  *	@param	test_what	[in] Determines if internal validation tests are performed.
	  *							The optimality conditions for the QP are not checked
	  *							internally, this is something that client can (and should)
	  *							do independently (see QPSolverRelaxedTester).
	  *							RUN_TESTS : As many validation/consistency tests
	  *								are performed internally as possible.  If a test
	  *								fails then a TestFailed execption will be thrown.
	  *								The subclasses determine what the tests are and
	  *								what failing a test means.
	  *							NO_TEST : No tests are performed internally.  This is
	  *								to allow the fastest possible execution.
	  *	@param	g			[in] vector (size n): objective 1st order
	  *	@param	G			[in] matrix (size n x n): objective, second order Hessian
	  *	@param	etaL		[in] scalar: Lower bound for relaxation variable (0 usually)
	  * @param	dL			[in] sparse vector (size n) (dL->is_sorted() == true): lower variable bounds
	  * @param	dU			[in] sparse vector (size n) (dU->is_sorted() == true): upper variable bounds
	  *	@param	E			[in] matrix (op(E)) (size m_in x n): inequality constraint Jacobian matrix 
	  * @param	trans_E		[in] E is transposed? 
	  * @param	b			[in] vector (size m_in): relaxation vector for inequalities
	  * @param	eL			[in] sparse vector (size m_in) (eL->is_sorted() == true): lower inequality bounds
	  * @param	eU			[in] sparse vector (size m_in) (eU->is_sorted() == true): upper inequality bounds
	  *	@param	F			[in] matrix (op(F) size m_eq x n): equality constraint Jacobian matrix
	  * @param	trans_F		[in] F is transposed? 
	  * @param	f			[in] vector (size m_eq): equality constraint right hand side
	  *	@param	obj_d		[out] If obj_d != NULL on input, then obj_d will be set with the
	  *							value of obj_d = g'*d + 1/2*d'*G*d for the value of d
	  *							returned.
	  *	@param	eta			[out] scalar:  Relaxation variable
	  *	@param	d			[in/out] vector (size n):  On input, it contains an intial estimate
	  *                         of the solution.  On output it is the estimate of the solution
	  *                         (see the return value).
	  *	@param	nu			[in/out] sparse vector (size n):  Lagrange multipilers
	  *							for variable bounds.  On input it contains the
	  *							estimate of the active set and multiplier values.  If
	  *							nu->nz() == 0 on input then there is no estimate for the
	  *							active-set.  Note that having nu->nz() > 0 on input is not
	  *							a commandment to perform a warm start.  Ultimatly this
	  *							decision is up to the subclass and lower level
	  *							subclass/client interactions.
	  *							On output nu contains the active-set for the
	  *							returned solution.  If nu->nz() > 0 on input then
	  *							nu->is_sorted() must be true.  On output nu->is_sorted()
	  *							will be true.
	  *	@param	mu			[in/out] sparse vector (size m_in):  Lagrange multipilers
	  *							for the general inequality constriants.
	  *							On input it contains the
	  *							estimate of the active set and multiplier values.  If
	  *							mu->nz() == 0 on input then there is no estimate for the
	  *							active-set.  Note that having mu->nz() > 0 on input is not
	  *							a commandment to perform a warm start.  Ultimatly this
	  *							decision is up to the subclass and lower level
	  *							subclass/client interactions.
	  *							On output mu contains the active-set for the
	  *							returned solution.  If mu->nz() > 0 on input then
	  *							mu->is_sorted() must be true.  On output mu->is_sorted()
	  *							will be true.
	  *	@param	Ed			[in/out] vector (size m_in) If Ed!=NULL on input, then on output
	  *							Ed will contain the product opt(E)*d for the value of
	  *							d on output.  This is included to save user from having
	  *							to perform this computation again if it has already been
	  *							done internally.
	  *	@param	lambda		[out] vector (size m_eq):  Lagrange multipilers for equality
	  *							constraints.
	  *	@param	Fd			[in/out] vector (size m_eq) If Fd!=NULL on input, then on output
	  *							Fd will contain the product opt(F)*d for the value of
	  *							d on output.  This is included to save user from having
	  *							to perform this computation again if it has already been
	  *							done internally.
	  *
	  *	@return
	  *		#OPTIMAL_SOLUTION# : Returned point satisfies the optimality conditions
	  *			in (2)-(6) above.  This will generally be the case if the maximum
	  *			number of QP iterations is not exceeded and none of the possible
	  *			exeptions are thrown.
	  *		#PRIMAL_FEASIBLE_POINT# : Returned point satisfies the feasibility
	  *         and complementarity conditions in (4)-(5) above.
	  *			For example, a primal, active-set QP algorithm may return this if
	  *			the maxinum number of iterations has been exceeded durring the
	  *			optimality phase (phase 2).
	  *		#DUAL_FEASIBLE_POINT# : Returned point satisfies the optimality conditions
	  *			in (2)-(4.1),(4.4),(5) and (6) but not the inequality constraints in (4.2),(4.3).
	  *			For example, a dual, active-set QP algorithm might return this
	  *			value if the maximum number of iterations has been exceeded.
	  *		#SUBOPTIMAL_POINT# : Returned point does not accurately enough satisfy any of the
	  *			optimality conditions in (2)-(6) above but the solution may still be
	  *			of some use.  For example, an active-set (primal, or dual) QP algorithm
	  *			might return this if there is some serious illconditioning in the QP
	  *			and a solution satisfying the desired tolerance could not be found.
	  *			Also, a primal-dual interior point QP solver might return this if the
	  *         maxinum number of iterations is exceeded.  The returned solution may
	  *         still be of some use to the client though.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
			, const SpVectorSlice& eL, const SpVectorSlice& eU
		, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
		, value_type* obj_d
		, value_type* eta, VectorSlice* d
		, SpVector* nu
		, SpVector* mu, VectorSlice* Ed
		, VectorSlice* lambda, VectorSlice* Fd
		);

	///
	/** Solve the QP without general equality constrants.
	  *
	  * By default this function calls #solve_qp(...)# which accepts
	  * various sets of constraints.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
			, const SpVectorSlice& eL, const SpVectorSlice& eU
		, value_type* obj_d
		, value_type* eta, VectorSlice* d
		, SpVector* nu
		, SpVector* mu, VectorSlice* Ed
		);

	///
	/** Solve the QP without general inequality constrants.
	  *
	  * By default this function calls #solve_qp(...)# which accepts
	  * various sets of constraints.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
		, value_type* obj_d
		, value_type* eta, VectorSlice* d
		, SpVector* nu
		, VectorSlice* lambda, VectorSlice* Fd
		);


	///
	/** Solve the QP without general equality or inequality constrants (no relaxation
	  * needed).
	  *
	  * By default this function calls #solve_qp(...)# which accepts
	  * various sets of constraints.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, const VectorSlice& g, const MatrixWithOp& G
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, value_type* obj_d
		, VectorSlice* d
		, SpVector* nu
		);

	///
	/** This is a more flexible function where the client can
	  * set different constraints to be included.
	  *
	  * The defualt implementation of this function validates that input
	  * sizes are correct and that the proper sets of constraints are set
	  * by calling validate_input(...) first.  Refere to the method
	  * \Ref{validate_input}(...) to see how arguments are set.  After validation
	  * print_qp_input(...) is called to print the QP intput arguments.  Then
	  * imp_solve_qp(...) is called which must be implemented by the subclass.
	  * Finally, print_qp_output(...) is called to print the QP output.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, value_type* obj_d
		, value_type* eta, VectorSlice* d
		, SpVector* nu
		, SpVector* mu, VectorSlice* Ed
		, VectorSlice* lambda, VectorSlice* Fd
		);

	///
	/** Get the statistics of the last QP solved.
	  *
	  *	solution:	Returns the type of solution found.  See solve_qp(...)\\
	  *	num_iter:	Gives the number of QP iterations on output.\\
	  *	num_adds:	Gives the number of iterations where a variable
	  *				was added to the active-set.  This does not include
	  *				variables that where part of the initial estimate in nu
	  *				for a warm start.\\
	  *	num_drops:	Gives the number QP iterations where a variable was
	  *				dropped.\\
	  *	warm_start:	Returns if a warm start was performed.\\
	  *	infeas_qp:	Returns if eta > 0.0.
	  */
	virtual QPSolverStats get_qp_stats() const = 0;

	///
	/** Release any memory that is being used.
	  */
	virtual void release_memory() = 0;

	///
	/** This is (static) function that is used as a utility to
	  * validate the input arguments to solve_qp(...).
	  * 
	  * The input arguments are validated as follows.
	  * 
	  * If the inequality constraints are excluded then:
	  * void(E) == void(b) == void(eL) == void(eU) == void(mu) == NULL
	  *
	  * If the equality constraints are excluded then:
	  * void(F) == void(f) == void(lambda) == NULL
	  *
	  * If the equality and inequality constraints are excluded then:
	  * eta == NULL
	  * void(E) == void(b) == void(eL) == void(eU) == void(mu) == NULL
	  * void(F) == void(f) == void(lambda) == NULL
	  * 
	  * If any errors are found an std::invalid_argument exception
	  * will be thrown.
	  */
	static void validate_input(
		  const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, const value_type* obj_d
		, const value_type* eta, const VectorSlice* d
		, const SpVector* nu
		, const SpVector* mu, const VectorSlice* Ed
		, const VectorSlice* lambda, const VectorSlice* Fd
		);

	/** @name Utility functions for dumping input and output argements.
	 *
	 * The amount of output is determined by #out# and #olevel#.
	 *
	 * @param out [out] stream printed to if #out != NULL#
	 * @param olevel [in] Determines what is printed.
	 *   \begin{description}
	 *   \item[(int)olevel >= (int)PRINT_ITER_STEPS] Prints O(1) information about the arguments.
	 *   \item[(int)olevel >= (int)PRINT_ITER_ACT_SET] Prints the contents of #nu#, #mu#, and #lambda#.
	 *      Output is proportional to the number of active constraints O(nu->nz() + mu->nz() + lambda->size()). 
	 *   \item[(int)olevel >= (int)PRINT_ITER_VECTORS] Prints the contents of all the vectors.
	 *      Output is proportional to O(d->size()).
	 *   \item[(int)olevel >= (int)PRINT_EVERY_THING] Prints the contents of all the vectors and matrices.
	 *      Output could be as large as O(d->size() * d->size()) or larger.
	 *   \end{description}
	 *  
	 */
	//@{

	///
	/** Utility (static) function for printing the input input/output arguments before
	  * the QP solver is run.  The QP solver subclasses can call this function.
	  */
	static void print_qp_input( 
		std::ostream* out, EOutputLevel olevel
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, value_type* eta, VectorSlice* d
		, SpVector* nu
		, SpVector* mu
		, VectorSlice* lambda
		);

	///
	/** Utility (static) function for printing the output input/output arguments after
	  * the QP solver is run.  The QP solver subclasses can call this function.
	  */
	static void print_qp_output(
		std::ostream* out, EOutputLevel olevel
		, const value_type* obj_d
		, const value_type* eta, const VectorSlice* d
		, const SpVector* nu
		, const SpVector* mu, const VectorSlice* Ed
		, const VectorSlice* lambda, const VectorSlice* Fd
		);

	//@}

protected:

	///
	/** Subclasses are to override this to implement the QP algorithm.
	  *
	  * Called by default implementations of solve_qp(...) methods.
	  */
	virtual QPSolverStats::ESolutionType imp_solve_qp(
		  std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, value_type* obj_d
		, value_type* eta, VectorSlice* d
		, SpVector* nu
		, SpVector* mu, VectorSlice* Ed
		, VectorSlice* lambda, VectorSlice* Fd
		) = 0;

};	// end class QPSovlerRelaxed

}	// end namespace ConstrainedOptimizationPack

#endif	// QP_SOLVER_RELAXED_H
