// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxedTester.h

#ifndef QP_SOLVER_RELAXED_TESTER_H
#define QP_SOLVER_RELAXED_TESTER_H

#include "QPSolverRelaxed.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Tests the optimality conditions of the output from a QPSolverRelaxed
  * object.
  *
  * For the given QP and is solution (if solved) this class tests
  * the optimality conditions.
  * 
  * The optimality conditions checked are:
  \begin{verbatim}

	Linear dependence of gradients:
	  
	(2)  d(L)/d(d) = g + G*d - nuL + nuU + op(E)'*(- muL + muU) + op(F)'*lambda
	               = g + G*d + nu + op(E)'*mu - op(F)'*lambda = 0
	  
	     where: nu = nuU - nuL, mu = muU - muL

	Feasibility:
	  
	(4.1)  etaL <=  eta
	(4.2)  dL   <=  d                       <= dU
	(4.3)  eL   <=  op(E)*d + b*eta         <= eU
	(4.4)  op(F)*d + (1 - eta) * f  = 0

	Complementarity:

	(5.1)  nu(i) * (dL - d)(i), if nu(i) <= 0, i = 1...n
	(5.2)  nu(i) * (d - dU)(i), if nu(i) >= 0, i = 1...n
	(5.3)  mu(j) * (eL - op(E)*d - b*eta)(j), if mu(j) <= 0, j = 1...m_in
	(5.4)  mu(j) * (op(E)*d + b*eta - eU)(j), if mu(j) >= 0, j = 1...m_in
 
  \end{verbatim}
  *
  * The realtive error of each of these conditions is checked.  Specifically,
  * here is how the errors are computed which are compared to the error and warning
  * tolerances:
  \begin{verbatim}
  
    opt_err = | g + G*d + nu + op(E)'*mu - op(F)'*lambda |
               / ( 1 + ||g||inf + ||G*d||inf + ||nu||inf + ||op(E)'*mu||inf
                  + ||op(F)'*lambda||inf )
                  
    feas_err = ( b - op(A)*x ) / ( 1 + ||op(A)*x||inf )

    comp_err(i) = gamma(i) * ( op(A)*x - b )(i), for gamma(i) != 0

    	where: op(A)*x <= b

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
  * The complementarity conditions (5.1)-(5.4) are also checked.  These will should
  * be satisfied for any solution type other than a SUBOPTIMAL_POINT return value from
  * solve_qp(...).  Checking the complementarity conditions for an active set
  * QP solver is just checking that the active constraints are satisfied.
  * By scaling this violation by the Langrange multiplier you emphasis the
  * feasibility of those constraints that have the greatest impact on the
  * objective function.  In other words, all things being equal, we are
  * more concerned with a tight feasibility tolerance for constraints with
  * larger lagrange multipliers than for those with smaller multipliers.
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
	  *						[I]	Value returned from QPSolverRelaxed::solve_qp(...).
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
	  *	@param	out			[O]	If != NULL, the output is sent to this stream.
	  *	@param	print_all_warnings
	  *						[I]	If true, then all errors greater than warning_tol will
	  *							be printed.
	  *	@param	g			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	G			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	etaL		[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	dL			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	dU			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	E			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	trans_E		[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	b			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	eL			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	eU			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	F			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	trans_F		[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	f			[I]	Input to QPSolverRelaxed::solve_qp(...).
	  *	@param	obj_d		[I]	Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	eta			[I]	Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	d 			[I]	Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	nu			[I]	Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	mu 			[I]	Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	Ed			[I]	Output from QPSolverRelaxed::solve_qp(...).
	  *	@param	Fd			[I]	Output from QPSolverRelaxed::solve_qp(...).
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
