// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxedTester.h

#ifndef QP_SOLVER_RELAXED_TESTER_H
#define QP_SOLVER_RELAXED_TESTER_H

#include "QPSolverRelaxed.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Tests the optimality conditions of the output from a \Ref{QPSolverRelaxed}
  * object.
  *
  * For the given QP and is solution (if solved) this class tests
  * the optimality conditions.
  * 
  * The optimality conditions checked are:
  \begin{verbatim}

	Linear dependence of gradients:
	  
	(2)  d(L)/d(d) = g + G*d - nuL + nuU + op(E)'*(- muL + muU) + op(F)'*lambda
	               = g + G*d + nu + op(E)'*mu + op(F)'*lambda = 0
	  
	     where: nu = nuU - nuL, mu = muU - muL

	Feasibility:
	  
	(4.1)  etaL <=  eta
	(4.2)  dL   <=  d                       <= dU
	(4.3)  eL   <=  op(E)*d - b*eta         <= eU
	(4.4)  op(F)*d + (1 - eta) * f  = 0

	Complementarity:

	(5.1)  nu(i) * (dL - d)(i), if nu(i) <= 0, i = 1...n
	(5.2)  nu(i) * (d - dU)(i), if nu(i) >= 0, i = 1...n
	(5.3)  mu(j) * (eL - op(E)*d + b*eta)(j), if mu(j) <= 0, j = 1...m_in
	(5.4)  mu(j) * (op(E)*d - b*eta - eU)(j), if mu(j) >= 0, j = 1...m_in
 
  \end{verbatim}
  *
  * The realtive error of each of these conditions is checked.  Specifically,
  * here is how the errors are computed which are compared to the error and warning
  * tolerances:
  \begin{verbatim}
  
    opt_err = | g + G*d + nu + op(E)'*mu + op(F)'*lambda | / (1 + opt_scale)
                  
    feas_err = ( b - op(A)*x ) / ( 1 + ||op(A)*x||inf )

    comp_err(i) = gamma(i) * ( op(A)*x - b )(i) / ( 1 + opt_scale + ||op(A).row(i)'*x||inf )
        ,for gamma(i) != 0

   	where:
	    op(A)*x <= b
	    opt_scale = ||g||inf + ||G*d||inf + ||nu||inf + ||op(E)'*mu||inf + ||op(F)'*lambda||inf

  \end{verbatim}
  *
  * Above, op(A)*x <= b can represent any of the constraints in (4.1)-(4.4).
  *
  * Any elements of opt_err(i) >= opt_warning_tol will result in an error
  * message printed to *out.  Any elements of opt_err(i) >= opt_error_tol
  * will cause the checks to stop and false to be returned from the function
  * check_optimiality_conditions(...).
  *
  * Any elements of feas_err(i) >= feas_warning_tol will result in an error
  * message printed to *out.  Any elements of feas_err(i) >= feas_error_tol
  * will cause the checks to stop and false to be returned from the function
  * check_optimiality_conditions(...).
  *
  * Any elements of comp_err(i) >= comp_warning_tol will result in an error
  * message printed to *out.  Any elements of comp_err(i) >= comp_error_tol
  * will cause the checks to stop and false to be returned from the function
  * check_optimiality_conditions(...).
  *
  * The goal of these tests is to first and foremost to catch gross programming
  * errors.  These tests can also be used to help flag and catch illconditioning
  * in the problem or instability in the QP solver.  The importance of such tests
  * can not be overestimated.  The scalings above are done to try to adjust for
  * the scaling of the problem.  Note that we are accounting for very big numbers
  * but not for very small numbers very well and therefore tests may be conserative
  * in some cases.  At the very least we account for loss of precision due to 
  * catastrophic cancelation that occurs when adding large numbers and expecting
  * to get zero. 
  *
  * As shown above, the complementarity conditions (5.1)-(5.4) are specifically checked.
  * These should be satisfied for any solution type other than a SUBOPTIMAL_POINT
  * return value from solve_qp(...).  Checking the complementarity conditions for
  * an active set QP solver is just checking that the active constraints are satisfied.
  * Checking them for an Iterior Point solver is critical to ensure that the system was
  * solved to satisfactory tolerance.
  * By scaling the active constraint violation by the Langrange multiplier
  * we emphasis the feasibility of those constraints that have the greatest
  * impact on the objective function.  In other words, all things being equal, we are
  * more concerned with a tight feasibility tolerance for constraints with
  * larger lagrange multipliers than for those with smaller multipliers.
  * The complementarity error is also scaled by the inverse of the sums
  * of the optimality scaling opt_scale and the size of the constraint residual.
  * By scaling by the max term opt_scale in the linear dependence of gradients we are
  * trying to adjust the effect of the lagrange multiplier.  Therefore if the gradient
  * of the objective g+G*d is large then opt_scale will account for this.
  */
class QPSolverRelaxedTester {
public:

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_error_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_error_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, comp_warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, comp_error_tol )

	///
	QPSolverRelaxedTester(
		  value_type	opt_warning_tol		= 1e-10
		, value_type	opt_error_tol		= 1e-5
		, value_type	feas_warning_tol	= 1e-10
		, value_type	feas_error_tol		= 1e-5
		, value_type	comp_warning_tol	= 1e-10
		, value_type	comp_error_tol		= 1e-5
		);

	///
	virtual ~QPSolverRelaxedTester() {}

	///
	/** Check the optimality conditions for the solved (or partially solved) QP.
	  *
	  * The default implementation calls the function check_optimality_conditions(...)
	  * which accepts various sets of constraints.
	  *
	  *	@param	solution_type
	  *						[in] Value returned from QPSolverRelaxed::solve_qp(...).
	  *							The exact optimality conditions enforced is determined
	  *							by this argument.  All of the optimality conditions are
	  *							checked.
	  *								OPTIMAL_SOLUTION : All of the optimality conditions
	  *									are enforced.
	  *								PRIMAL_FEASIBLE_POINT : Only the optimality conditions
	  *								 	(4.1)-(4.4) are enforced.
	  *								DUAL_FEASIBLE_POINT: Only the optimality condtions
	  *								 	(2) and (6.1)-(6.4) are enforced.
	  *								SUBOPTIMAL_POINT : None of the optimality conditions
	  *								 	are enforced.
	  *	@param	out			[out] If != NULL, the output is sent to this stream.
	  *	@param	print_all_warnings
	  *						[in] If true, then all errors greater than warning_tol will
	  *							be printed.
	  *	@param	g			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	G			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	etaL		[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	dL			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	dU			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	E			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	trans_E		[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	b			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	eL			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	eU			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	F			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	trans_F		[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	f			[in] Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	obj_d		[in] Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	eta			[in] Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	d 			[in] Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	nu			[in] Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	mu 			[in] Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	Ed			[in] Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	Fd			[in] Output from QPSolverRelaxed::solve_qp(...).
	  *
	  * @return #true# if all of the errors are greater than the error tolerances
	  * 	, otherwise it returns #false#
	  */
	virtual bool check_optimality_conditions(
		  QPSolverStats::ESolutionType solution_type
		, std::ostream* out, bool print_all_warnings, bool print_vectors
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
			, const SpVectorSlice& eL, const SpVectorSlice& eU
		, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
		, const value_type* obj_d
		, const value_type* eta, const VectorSlice* d
		, const SpVector* nu
		, const SpVector* mu, const VectorSlice* Ed
		, const VectorSlice* lambda, const VectorSlice* Fd
		);

	///
	/** Check the optimality conditions without general equality constrants.
	  */
	virtual bool check_optimality_conditions(
		  QPSolverStats::ESolutionType solution_type
		, std::ostream* out, bool print_all_warnings, bool print_vectors
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp& E, BLAS_Cpp::Transp trans_E, const VectorSlice& b
			, const SpVectorSlice& eL, const SpVectorSlice& eU
		, const value_type* obj_d
		, const value_type* eta, const VectorSlice* d
		, const SpVector* nu
		, const SpVector* mu, const VectorSlice* Ed
		);

	///
	/** Check the optimality conditions general inequality constrants.
	  */
	virtual bool check_optimality_conditions(
		  QPSolverStats::ESolutionType solution_type
		, std::ostream* out, bool print_all_warnings, bool print_vectors
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp& F, BLAS_Cpp::Transp trans_F, const VectorSlice& f
		, const value_type* obj_d
		, const value_type* eta, const VectorSlice* d
		, const SpVector* nu
		, const VectorSlice* lambda, const VectorSlice* Fd
		);


	///
	/** Check the optimality conditions without general equality or inequality
	  * constrants (no relaxation needed).
	  */
	virtual bool check_optimality_conditions(
		  QPSolverStats::ESolutionType solution_type
		, std::ostream* out, bool print_all_warnings, bool print_vectors
		, const VectorSlice& g, const MatrixWithOp& G
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const value_type* obj_d
		, const VectorSlice* d
		, const SpVector* nu
		);

	///
	/** This is a more flexible function where the client can
	  * set different constraints to be included.
	  *
	  */
	virtual bool check_optimality_conditions(
		  QPSolverStats::ESolutionType solution_type
		, std::ostream* out, bool print_all_warnings, bool print_vectors
		, const VectorSlice& g, const MatrixWithOp& G
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

protected:

	///
	/** Subclasses are to override this to implement the
	  * testing code.
	  *
	  * There is a default implementation that is very general and
	  * should be considered good enough for most applications.
	  */
	virtual bool imp_check_optimality_conditions(
		  QPSolverStats::ESolutionType solution_type
		, std::ostream* out, bool print_all_warnings, bool print_vectors
		, const VectorSlice& g, const MatrixWithOp& G
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

};	// end class QPSolverRelaxedTester

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// QP_SOLVER_RELAXED_TESTER_H
