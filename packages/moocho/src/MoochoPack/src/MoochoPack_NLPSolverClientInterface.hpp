// ////////////////////////////////////////////////////////////////////////////
// rSQPSolverClientInterface.h

#ifndef RSQP_ALGO_SOLVER_CLIENT_INTERFACE_H
#define RSQP_ALGO_SOLVER_CLIENT_INTERFACE_H

#include <stdexcept>

#include "ReducedSpaceSQPPackTypes.h"
#include "rSQPTrack.h"
#include "NLPInterfacePack/include/NLPReduced.h"
#include "Misc/include/ref_count_ptr.h"

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
	enum EFindMinReturn { SOLUTION_FOUND, MAX_ITER_EXCEEDED, MAX_RUN_TIME_EXCEEDED };

	///
	typedef ReferenceCountingPack::ref_count_ptr<rSQPTrack>		track_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<NLPReduced>	nlp_ptr_t;

	/// Thrown if the setup is not valid
	class InvalidSetup : public std::logic_error
	{public: InvalidSetup(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	///
	/** Construct with no nlp reference set, track set to rSQPTrack
	  * , max_iter = 100, kkt_tol = 1e-6, feas_tol = 1e-6, and step_tol = 1e-6
	  * , max_var_bounds_viol = 1e-6.
	  */
	rSQPSolverClientInterface()
		: max_iter_(100), kkt_tol_(1e-6), feas_tol_(1e-6), step_tol_(1e-6)
			, max_var_bounds_viol_(1e-8)
	{}

	///
	virtual ~rSQPSolverClientInterface() {}

	/** @name «std comp» members for nlp. */
	//@{

	///
	virtual void set_nlp(const nlp_ptr_t& nlp)
	{	nlp_ = nlp; }
	///
	virtual nlp_ptr_t& get_nlp()
	{	return nlp_ ; }
	///
	virtual const nlp_ptr_t& get_nlp() const
	{	return nlp_; }
	///
	virtual NLPReduced& nlp()
	{	return *nlp_; }
	///
	virtual const NLPReduced& nlp() const
	{	return *nlp_; }

	//@}

	/** @name «std comp» members for track. */
	//@{

	///
	virtual void set_track(const track_ptr_t& track)
	{	track_ = track; }
	///
	virtual track_ptr_t& get_track()
	{	return track_ ; }
	///
	virtual const track_ptr_t& get_track() const
	{	return track_; }
	///
	virtual rSQPTrack& track()
	{	return *track_; }
	///
	virtual const rSQPTrack& track() const
	{	return *track_; }

	//@}

	/** @name Solver Parameters*/
	//@{

	/// Set the maximum number of iterations the rSQP algorithm can perform (default = 100)
	void set_max_iter(int max_iter)
	{	max_iter_ = max_iter; }
	///
	int max_iter() const
	{	return max_iter_; }

	/// Set the maximum run_time (in min, default = infinity)
	virtual void set_max_run_time(double max_run_time) = 0;
	///
	virtual double max_run_time() const = 0;

	///
	/** Set the termination tolerance for the relative KKT first order necessary conditions.
	  * Default = 1e-6.
	  */
	void kkt_tol(value_type kkt_tol)
	{	kkt_tol_ = kkt_tol; }
	///
	value_type kkt_tol() const
	{	return kkt_tol_; }

	///
	/** Set the termination tolerance for the equality constraints ||c(x*)||inf.
	  * Default = 1e-6.
	  */
	void feas_tol(value_type feas_tol)
	{	feas_tol_ = feas_tol; }
	///
	value_type feas_tol() const
	{	return feas_tol_; }

	///
	/** Set the termination tolerance for the change in the estimate of the solution.
	  * The test is: |d(i)|/(1+|x(i)|) < step_tol. 
	  * Default = 1e-6.
	  */
	void step_tol(value_type step_tol)
	{	step_tol_ = step_tol; }
	///
	value_type step_tol() const
	{	return step_tol_; }

	///
	/** Set the maximum absolute value for which the variable bounds may be violated
	  * by when computing function and gradient values.
	  *
	  * In other words the algorithm will never call an the NLP to compute
	  * a function and gradient evaluation outside of:
	  *
	  * xl - delta <= x <= xu + delta
	  *
	  * where delta = max_var_bounds_viol.
	  */
	void max_var_bounds_viol(value_type max_var_bounds_viol)
	{	max_var_bounds_viol_ = max_var_bounds_viol; }
	///
	value_type max_var_bounds_viol() const
	{	return max_var_bounds_viol_; }

	//@}

	/** @name Solve */
	//@{

	///
	/** Find the minimun of the set NLPReduced.
	  *
	  * If this function returns SOLUTION_FOUND then the solution vector to the NLPReduced can be obtained as
	  * #this->state().x().get_k(1)# and the value of the objective function is given as
	  * #this->state().f().get_k(1)#.
	  *
	  * If the maximum number of iterations has been exceeded then MAX_ITER_EXCEEDED will
	  * be returned.
	  *
	  * If the maximum runtime has been exceeded then MAX_RUN_TIME_EXCEEDED will be
	  * returned.
	  *
	  * Preconditions:\begin{itemize}
	  * \item #this->nlp() != 0# (throw InvalidSetup)
	  * \end{itemize}
	  *
	  * Postcondtions:\begin{itemized}
	  * \item	Minimum of NLPReduced is found to kkt_tol, max_iter was reached
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

private:
	int			max_iter_;	// The maximum number of iterations to perform
	value_type	kkt_tol_;
	value_type	feas_tol_;
	value_type	step_tol_;
	value_type	max_var_bounds_viol_;
	nlp_ptr_t	nlp_;		// Pointer to the nlp to solve.  Must be subclass of NLPReduced
	track_ptr_t	track_;		// Pointer to the track object that rSQPAlgo uses to output
							// itermediate results

};	// end class rSQPSolverClientInterface

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_SOLVER_CLIENT_INTERFACE_H