// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxed.h

#ifndef QP_SOLVER_RELAXED_H
#define QP_SOLVER_RELAXED_H

#include "QPSolverStats.h"

namespace ConstrainedOptimizationPack {

///
/** Solves Quadratic Programs (QPs) of the form:
  *
  * (1.1)  min          g'*x + 1/2*x'*G*x + M(eta)
  *        x <: R^n
  *
  * s.t.
  * (1.2)               etaL <=  eta
  * (1.3)               xL   <=  x                       <= xU
  *	(1.4)               eL   <=  op(E)*x + b*eta         <= eU
  *	(1.5)                        op(F)*x + (1 - eta) * f  = 0
  *
  * The relaxation is used to ensure that the QP will have a solution
  * (eta = 1, x = 0 guarrentted if xL <= 0 <= xU and eL <= b <= eU).
  * If the function of M(eta) in the objective is large enough, then
  * the constraint etaL <= eta will be active if a feasible region
  * exists.
  *
  * The Lagrangian for the QP is:
  *
  * L = g' * x + 1/2 * x' * G * x + M(eta)
  *		+ kappa * (etaL - eta)
  *		+ nuL' * (xL - x)
  *		+ nuU' * (x - xU)
  *		+ muL' * (eL - op(E)*x - b*eta)
  *		+ muU' * (op(E)*x + b*eta - eU)
  *		+ lambda' * (op(F)*x + (1 - eta) * f)
  *
  * The optimality conditions for this QP are:
  *
  * Linear dependence of gradients:
  *
  * (2)  d(L)/d(x) = g + G*x - nuL + nuU + op(E)'*(- muL + muU) + op(F)'*lambda
  *                = g + G*x + nu + op(E)'*mu - op(F)'*lambda = 0
  *
  *		where: nu = nuU - nuL, mu = muU - muL
  *
  * (3)  d(L)/d(eta) = d(M)/d(eta) - kappa + b'*mu - f'*lambda = 0
  *
  * Feasibility:
  *
  * (4.1)  etaL <=  eta
  * (4.2)  xL   <=  x                       <= xU
  *	(4.3)  eL   <=  op(E)*x + b*eta         <= eU
  *	(4.4)  op(F)*x + (1 - eta) * f  = 0
  *
  * The optimal x and eta are determined as well as the lagrange multipliers
  * for the constriants:
  *
  * nu :			xL <= x <= xU
  * 
  * mu :			eL <= op(E)*x + b*eta <= eU
  *
  * lambda :		op(F)*x + (1 - eta) * f  = 0
  *
  * The lagrange multiper for the constraint etaL <= eta kappa is not given
  * since if this constraint is not active, then kappa = 0 and all of the multiplier
  * estimates will be off because of the arbitrarily large value of M(eta) and
  * the optimality condition d(L)/d(eta) = 0.
  *
  * The bounds xL, xU, eL and eU are input as sparse vectors (SpVectorSlice).
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

	/// Enumeration for if the run internal tests or not.
	enum ERunTests { RUN_TESTS, NO_TESTS };

	//@}

	///
	virtual ~QPSolverRelaxed() {}

	///
	/** Solve the QP.
	  *
	  * This function may throw many exceptions.  If there is some problem with
	  * the QP definition the exceptions Unbounded, Infeasible or InvalidInput
	  * may be thrown.  If the QP is illconditioned other exeptions may be thrown
	  * by this function or perhaps a warning message may be printed to *out and
	  * a value of the return other than OPTIMAL_SOLUTION may be returned.
	  *
	  * After this function returns, get_qp_stats(...) can be called to return
	  * some statistics for the QP just solved (or attempted to be solved).
	  *
	  * Note skip the variable bounds can be removed by passing in
	  * pass in xL.nz() == xU.nz() == 0.
	  *
	  *	@param	g			[I]	vector (size n)
	  *	@param	G			[I]	matrix (size n x n).
	  *	@param	etaL		[I] scalar.
	  * @param	xL			[I] sparse vector (size n) (xL->is_sorted() == true)
	  * @param	xU			[I] sparse vector (size n) (xU->is_sorted() == true)
	  *	@param	E			[I] matrix (op(E) size m_in x n)
	  * @param	trans_E		[I]	E is transposed? 
	  * @param	b			[I] vector (size m_in)
	  * @param	eL			[I] sparse vector (size m_in) (eL->is_sorted() == true)
	  * @param	eU			[I] sparse vector (size m_in) (eU->is_sorted() == true)
	  *	@param	F			[I] matrix (op(F) size m_eq x n)
	  * @param	trans_F		[I]	F is transposed? 
	  * @param	f			[I] vector (size m_eq)
	  *	@param	out			[O]	If out != NULL then output is printed to this stream
	  *							depending on the value of olevel.
	  *	@param	olevel		[I]	Determines the amount of output to print to *out.
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
	  *	@param	test_what	[I]	Determines if internal validation tests are performed.
	  *							The optimality conditions for the QP are not checked
	  *							internally, this is something that client can (and should)
	  *							do independently (see TestQPSolverRelaxed).
	  *							RUN_TESTS : As many validation / consistency tests
	  *								are performed internally as possible.  If a test
	  *								fails then a TestFailed execption will be thrown.
	  *							NO_TEST : No tests are performed internally.  This is
	  *								to allow the fastest possilbe execution as possible.
	  *	@param	eta			[O]	scalar.  Relaxation variable.
	  *	@param	x			[O]	vector (size n).  Solution vector.
	  *	@param	nu			[I/O]	sparse vector (size n).  Lagrange multipilers
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
	  *							will be true..
	  *	@param	mu			[I/O]	sparse vector (size m_in).  Lagrange multipilers
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
	  *							will be true..
	  *	@param	lambda		[O]	vector (size m_eq).  Lagrange multipilers for equality
	  *							constraints.
	  *
	  *	@return
	  *		OPTIMAL_SOLUTION : Returned point satisfies the optimality conditions
	  *			in (2)-(4) above.  This will generally be the case if the maximum
	  *			number of QP iterations is not exceeded and none of the possible
	  *			exeptions are thrown.
	  *		PRIMAL_FEASIBLE_POINT : Returned point satisfies the feasibility conditions
	  *			in (4) above.
	  *			For example, a primal, active-set QP algorithm may return this if
	  *			the maxinum number of iterations has been exceeded durring the
	  *			optimality phase (phase 2).
	  *		DUAL_FEASIBLE_POINT : Returned point satisfies the optimality conditions
	  *			in (2)-(4.1) but not the variable bounds in (4.2),(4.3).
	  *			For example, a primal-dual, active-set QP algorithm might return this
	  *			value if the maximum number of iterations has been exceeded.
	  *		SUBOPTIMAL_POINT : Returned point does not accurately enough satisfy any of the
	  *			optimality conditions in (2)-(4) above but the solution may still be
	  *			of some use.  For example, an active-set (primal, or dual) QP algorithm
	  *			might return this if there is some serious illconditioning in the QP
	  *			and a solution satisfying the desired tolerance could not be found.
	  *			Also, an interior point QP solver might return this if the maxinum
	  *			number if iterations is exceeded.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& xL, const SpVectorSlice& xU
		, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
			, const SpVectorSlice& eL, const SpVectorSlice& eU
		, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
		, std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, value_type* eta, Vector* x, SpVector* nu, SpVector* mu, Vector* lambda
		) = 0;

	///
	/** Solve the QP without general equality constrants.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& xL, const SpVectorSlice& xU
		, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
			, const SpVectorSlice& eL, const SpVectorSlice& eU
		, std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, value_type* eta, Vector* x, SpVector* nu, SpVector* mu
		) = 0;

	///
	/** Solve the QP without general inequality constrants.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& xL, const SpVectorSlice& xU
		, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
		, std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, value_type* eta, Vector* x, SpVector* nu, Vector* lambda
		) = 0;


	///
	/** Solve the QP without general equality or inequality constrants (no relaxation
	  * needed).
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  const VectorSlice& g, const MatrixWithOp& G
		, const SpVectorSlice& xL, const SpVectorSlice& xU
		, std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, Vector* x, SpVector* nu
		) = 0;

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


protected:

	///
	/** Validates that input sizes are correct and the proper set of
	  * constraints is set.
	  *
	  * After checking the input arguments the pure virtual function
	  * imp_solve_qp(...) is called which must be implemented by the subclass.
	  */
	virtual QPSolverStats::ESolutionType solve_qp(
		  const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& xL, const SpVectorSlice& xU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, value_type* eta, Vector* x, SpVector* nu, SpVector* mu, Vector* lambda
		);
	
	///
	/** Subclasses are to override this to implement the QP algorithm.
	  *
	  * The sizes of the input arguments is validated as well as the proper
	  * sets of constraints.
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
	  */
	virtual QPSolverStats::ESolutionType imp_solve_qp(
		  const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& xL, const SpVectorSlice& xU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, value_type* eta, Vector* x, SpVector* nu, SpVector* mu, Vector* lambda
		) = 0;

};	// end class QPSovlerWithBounds

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// QP_SOLVER_RELAXED_H