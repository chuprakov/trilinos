// ////////////////////////////////////////////////////////////////////////////
// rSQPSolverClientInterface.h

#ifndef RSQP_ALGO_SOLVER_CLIENT_INTERFACE_H
#define RSQP_ALGO_SOLVER_CLIENT_INTERFACE_H

#include <stdexcept>

#include "ReducedSpaceSQPPackTypes.h"
#include "rSQPTrack.h"
#include "NLPInterfacePack/include/NLP.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** This is the interface that dumb clients use to set the NLPReduced and Track objects
  * , solver the NLPReduced and get the solution {abstract}.
  */
class rSQPSolverClientInterface {
public:

	/** @name Public Types */
	//@{

	///
	enum EFindMinReturn { SOLUTION_FOUND, MAX_ITER_EXCEEDED, MAX_RUN_TIME_EXCEEDED
		, ALGORITHMIC_ERROR };

	/// Thrown if the setup is not valid
	class InvalidSetup : public std::logic_error
	{public: InvalidSetup(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}
	
	/// «std comp» members for nlp
	STANDARD_COMPOSITION_MEMBERS( NLP, nlp )

	/// «std comp» members for track
	STANDARD_COMPOSITION_MEMBERS( rSQPTrack, track )

	/** @name Solver Parameters*/
	//@{

	/// Set the maximum number of iterations the rSQP algorithm can perform (default = 100)
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, max_iter )

	///
	/** Set the maximum run_time (in min, default = infinity).
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, max_run_time )

	///
	/** Set the termination tolerance for the relative (scaled) linear dependence of the
	  * gradients part of the first order necessary optimality conditions.
	  * Default = 1e-6.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_tol )

	///
	/** Set the termination tolerance for the (scaled) equality constraints ||c(x*)||inf
	  * which is part of the first order necessary optimality conditions.
	  * Default = 1e-6.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_tol )

	///
	/** Set the termination tolerance for the change in the estimate of the solution.
	  * The test is: |d(i)|/(1+|x(i)|) < step_tol. 
	  * Default = 1e-6.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, step_tol )

	///
	/** Set the maximum absolute value for which the variable bounds may be violated
	  * by when computing function and gradient values.
	  *
	  * In other words the algorithm will never call on the NLP to compute
	  * a function and gradient evaluation outside of:
	  *
	  * xl - delta <= x <= xu + delta
	  *
	  * where delta = max_var_bounds_viol.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_var_bounds_viol )

	///
	/** Determine the amount of output to a journal file.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EJournalOutputLevel, journal_output_level )

	///
	/** Set the precesion of the journal output.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, journal_print_digits )

	///
	/** Set whether computations will be double checked or not.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_results )

	//@}

	///
	/** Construct with no nlp reference set, track set to rSQPTrack
	  * , max_iter = 100, opt_tol = 1e-6, feas_tol = 1e-6, and step_tol = 1e-6
	  * , max_var_bounds_viol = 1e-6.
	  */
	rSQPSolverClientInterface(
		  int					max_iter				= 100
		, double				max_run_time			= 1e+10 // run forever
		, value_type			opt_tol					= 1e-6
		, value_type			feas_tol				= 1e-6
		, value_type			step_tol				= 1e-2
		, value_type			max_var_bounds_viol		= 1e-5
		, EJournalOutputLevel 	journal_output_level 	= PRINT_NOTHING
		, int					journal_print_digits	= 6
		, bool					check_results			= false
		);

	///
	virtual ~rSQPSolverClientInterface() {}

	/** @name Solve */
	//@{

	///
	/** Find the minimun of the set NLP.
	  *
	  * This function returns #SOLUTION_FOUND# if the NLP has been solved to the desired
	  * tolerances.  In this case #this->track().output_final(...,TERMINATE_TRUE)#
	  * and #this->nlp().report_final_solution(...,true)# is called.
	  *
	  * If the solution is not found, then #this->nlp().report_final_solution(...,false)# is called.
	  * and one of the following occurs:
	  *
	  * If the maximum number of iterations has been exceeded then #MAX_ITER_EXCEEDED# will
	  * be returned.  In this case #this->track().output_final(...,MAX_ITER_EXCEEDED)# is called.
	  *
	  * If the maximum runtime has been exceeded then #MAX_RUN_TIME_EXCEEDED# will be
	  * returned.  In this case #this->track().output_final(...,MAX_RUN_TIME_EXCEEDED)# is called.
	  *
	  * The client should be prepaired to catch any exceptions thrown from this function.
	  * All of the purposefully thrown exceptions are derived from std::exception so the
	  * client can check the what() function to see and description of the error.  If an
	  * exception is thrown then #this->track().output_final(...,TERMINATE_FALSE)# will
	  * be called before this exception is rethrown out of this functiion.  If the constraints
	  * are found to be infeasible, then the exception #InfeasibleConstraints# will be thrown.
	  * If a line search failure occurs then the exception #LineSearchFailure# will be thrown.
	  * If some test failed then the exception #TestFailed# will be thrown.  Many other exceptions
	  * may be thrown but these are the main ones that the SQP algorithms known about and will
	  * purposefully generate.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #this->nlp() != 0# (throw InvalidSetup)
	  * \end{itemize}
	  *
	  * Postcondtions:\begin{itemized}
	  * \item	Minimum of NLPReduced is found to opt_tol, max_iter was reached
	  *			or max_run_time reached (throw std::exection)
	  * \end{itemize}
	  */
	virtual EFindMinReturn find_min() = 0;

	///
	/** Returns the state object containing the iteration quantities at the solution.
	  */
	virtual const rSQPState& state() const = 0;

	//@}

	///
	/** Prints a description of the algorithm.
	  */
	virtual void print_algorithm(std::ostream& out) const = 0;

	/** @name Algorithm Timing.
	  */
	//@{

	///
	/** Causes algorithm to be timed.
	  *
	  * Call with #algo_timing == true# before #find_min(...)# to have 
	  * the algorithm timed.
	  */
	virtual void set_algo_timing( bool algo_timing ) = 0;

	///
	virtual bool algo_timing() const = 0;

	///
	/** Outputs table of times for each step and the cummulative times.
	  *
	  * Call after #find_min(...)# has executed to get a table
	  * of times.
	  */
	virtual void print_algorithm_times( std::ostream& out ) const = 0;

	//@}

};	// end class rSQPSolverClientInterface

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_SOLVER_CLIENT_INTERFACE_H
